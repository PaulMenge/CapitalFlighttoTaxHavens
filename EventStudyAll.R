suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(lubridate); library(zoo)
  library(ggplot2)
})

# -------------------------
# 0) Parameters
# -------------------------
W <- 10
shock_qidx <- 2011L*4 + 1L - 1L   # 2011Q1

# -------------------------
# 1) Load & normalize
# -------------------------
dat <- read_excel("C:/Users/paulm/Desktop/Public Projekt/Long.xlsx", sheet = "Sheet1") |>
  rename_with(tolower) |>
  mutate(
    date  = as.Date(date),
    qidx  = year(date)*4 + quarter(date) - 1L
  )

# -------------------------
# 2) Collapse to one obs per (country, quarter), ignoring TH inside country
# -------------------------
country_q <- dat |>
  group_by(country, qidx) |>
  summarise(
    value = sum(value, na.rm = TRUE),      # sum within country×quarter across TH cells
    treat = as.integer(any(treat == 1L)),  # country treated that quarter?
    .groups = "drop"
  )

# Ever-treated flag per country
ever_tab <- country_q |>
  group_by(country) |>
  summarise(ever = as.integer(any(treat == 1L)), .groups = "drop")

country_q <- country_q |>
  left_join(ever_tab, by = "country") |>
  mutate(k = qidx - shock_qidx)

# -------------------------
# 3) AGGREGATE FIRST: sum across countries within each group×quarter
# -------------------------
gq_sum <- country_q |>
  group_by(ever, qidx) |>
  summarise(total = sum(value, na.rm = TRUE), .groups = "drop") |>
  mutate(k = qidx - shock_qidx) |>
  filter(k >= -W, k <= W)

# Ensure both groups have a baseline at k = -1
base <- gq_sum |> filter(k == -1) |> select(ever, base_total = total)
stopifnot(nrow(base) == 2, all(base$base_total > 0))

# Relative-to-baseline series and percent change
gq_rel <- gq_sum |>
  left_join(base, by = "ever") |>
  mutate(
    ln_rel = log(total) - log(base_total),   # log-point change vs group's own 2010Q4
    pct    = 100*(exp(ln_rel) - 1),
    group  = ifelse(ever == 1, "Treated", "Control")
  )

# -------------------------
# 4) Two-line plot: % change from each group's baseline (aggregate-first)
# -------------------------
ggplot(gq_rel, aes(k, pct, color = group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_line(size = 1) + geom_point(size = 1.6) +
  scale_color_manual(values = c("Arab countries"="#1f9ac9","Other"="#fb8161")) +
  labs(
    title = "Aggregate totals % change vs group’s 2010Q4",
    x = "Quarters relative 2011/Q1",
    y = "% change from 2010/Q4",
    color = ""
  ) +
  theme_minimal(base_size = 12)

# -------------------------
# 5) Event-time difference (treated − control) with symmetric 95% CIs
#    d_k = [ln Tot_T(k) − ln Tot_T(-1)] − [ln Tot_C(k) − ln Tot_C(-1)]
# -------------------------
# Helper to compute the whole d_k vector given a country sample
compute_diff_series <- function(df_country_q) {
  gq <- df_country_q |>
    group_by(ever, qidx) |>
    summarise(total = sum(value, na.rm = TRUE), .groups = "drop") |>
    mutate(k = qidx - shock_qidx) |>
    filter(k >= -W, k <= W)
  
  base2 <- gq |> filter(k == -1) |> select(ever, base_total = total)
  if (nrow(base2) < 2 || any(base2$base_total <= 0)) return(NULL)
  
  rel <- gq |>
    left_join(base2, by = "ever") |>
    mutate(ln_rel = log(total) - log(base_total)) |>
    select(ever, k, ln_rel) |>
    tidyr::pivot_wider(names_from = ever, values_from = ln_rel, names_prefix = "ever_")
  
  if (!all(c("ever_0","ever_1") %in% names(rel))) return(NULL)
  rel <- rel[order(rel$k), ]
  rel$diff <- rel$ever_1 - rel$ever_0
  rel[, c("k","diff")]
}

# Original (point estimate)
orig <- compute_diff_series(country_q)
stopifnot(!is.null(orig))
orig_diff <- orig

# Bootstrap by resampling countries within each group
set.seed(123)
B <- 500
countries_T <- unique(country_q$country[country_q$ever == 1])
countries_C <- unique(country_q$country[country_q$ever == 0])
nT <- length(countries_T); nC <- length(countries_C)

boot_mat <- matrix(NA_real_, nrow = nrow(orig_diff), ncol = B)
rownames(boot_mat) <- orig_diff$k

for (b in 1:B) {
  samp_T <- sample(countries_T, nT, replace = TRUE)
  samp_C <- sample(countries_C, nC, replace = TRUE)
  df_b <- dplyr::bind_rows(
    country_q %>% filter(country %in% samp_T),
    country_q %>% filter(country %in% samp_C)
  )
  d_b <- compute_diff_series(df_b)
  if (!is.null(d_b)) {
    # align by k
    idx <- match(orig_diff$k, d_b$k)
    boot_mat[, b] <- d_b$diff[idx]
  }
}

# Drop columns with any NA (failed resamples)
boot_mat <- boot_mat[, colSums(is.na(boot_mat)) == 0, drop = FALSE]
if (ncol(boot_mat) == 0) stop("All bootstrap resamples failed (NA columns). Check data/baseline coverage.")

# --- SYMMETRIC CIs via bootstrap SEs + delta method (percent scale) ---
z <- qnorm(0.975)

# log-point effect and its bootstrap SD
b_log <- orig_diff$diff
se_log <- apply(boot_mat, 1, sd, na.rm = TRUE)

# percent effect and symmetric SE (delta method)
estPct <- 100 * (exp(b_log) - 1)
sePct  <- 100 * exp(b_log) * se_log

ci_tab <- data.frame(
  k      = orig_diff$k,
  b_log  = b_log,
  se_log = se_log,
  estPct = estPct,
  sePct  = sePct,
  loPct  = estPct - z * sePct,
  hiPct  = estPct + z * sePct
)
ci_tab$sigpos <- (ci_tab$k >= 0) & (ci_tab$loPct > 0)


# =========================
# Placebos (non-EU sample) + post-period p-values
# =========================

# Reusable helper: event-study diff relative to an arbitrary placebo shock
compute_diff_series_with_shock <- function(df_country_q, shock_idx) {
  gq <- df_country_q |>
    dplyr::group_by(ever, qidx) |>
    dplyr::summarise(total = sum(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(k = qidx - shock_idx) |>
    dplyr::filter(k >= -W, k <= W)

  base2 <- gq |> dplyr::filter(k == -1) |> dplyr::select(ever, base_total = total)
  if (nrow(base2) < 2 || any(base2$base_total <= 0)) return(NULL)

  rel <- gq |>
    dplyr::left_join(base2, by = "ever") |>
    dplyr::mutate(ln_rel = log(total) - log(base_total)) |>
    dplyr::select(ever, k, ln_rel) |>
    tidyr::pivot_wider(names_from = ever, values_from = ln_rel, names_prefix = "ever_")
  if (!all(c("ever_0","ever_1") %in% names(rel))) return(NULL)

  rel <- rel[order(rel$k), ]
  data.frame(k = rel$k, estPct = 100*(exp(rel$ever_1 - rel_0) - 1))
}

# (typo fix) correct compute_diff_series_with_shock:
compute_diff_series_with_shock <- function(df_country_q, shock_idx) {
  gq <- df_country_q |>
    dplyr::group_by(ever, qidx) |>
    dplyr::summarise(total = sum(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(k = qidx - shock_idx) |>
    dplyr::filter(k >= -W, k <= W)

  base2 <- gq |> dplyr::filter(k == -1) |> dplyr::select(ever, base_total = total)
  if (nrow(base2) < 2 || any(base2$base_total <= 0)) return(NULL)

  rel <- gq |>
    dplyr::left_join(base2, by = "ever") |>
    dplyr::mutate(ln_rel = log(total) - log(base_total)) |>
    dplyr::select(ever, k, ln_rel) |>
    tidyr::pivot_wider(names_from = ever, values_from = ln_rel, names_prefix = "ever_")
  if (!all(c("ever_0","ever_1") %in% names(rel))) return(NULL)

  rel <- rel[order(rel$k), ]
  data.frame(k = rel$k, estPct = 100*(exp(rel$ever_1 - rel$ever_0) - 1))
}

# Build all valid placebos using ONLY pre-event quarters in the non-EU sample
pre_candidates <- sort(unique(country_q_neu$qidx[country_q_neu$qidx < shock_qidx]))

placebos_df <- dplyr::bind_rows(lapply(pre_candidates, function(pq) {
  ds <- compute_diff_series_with_shock(country_q_neu, pq)
  if (is.null(ds)) return(NULL)
  ds$shock_qidx <- pq
  ds
}))

if (nrow(placebos_df) == 0) {
  warning("No valid placebo shocks found (non-EU sample). Skipping placebo p-values.")
} else {
  cat("Placebo shocks built (non-EU):", length(unique(placebos_df$shock_qidx)), "\n")

  # Main post-period (from ci_tab already built above)
  post_main <- ci_tab |> dplyr::filter(k >= 0)

  # Post-period of all placebos
  placebo_post <- placebos_df |> dplyr::filter(k >= 0)

  # Per-k empirical p-values
  pvals_k <- placebo_post |>
    dplyr::group_by(k) |>
    dplyr::summarise(
      pval = mean(abs(estPct) >= abs(post_main$estPct[match(k, post_main$k)]), na.rm = TRUE),
      .groups = "drop"
    )
  cat("\nEmpirical placebo p-values by k (post only, non-EU):\n")
  print(pvals_k)

  # Family-wise (max-stat) p-value across all post-k
  max_null <- placebo_post |>
    dplyr::group_by(shock_qidx) |>
    dplyr::summarise(max_abs = max(abs(estPct), na.rm = TRUE), .groups = "drop")
  main_max_abs <- max(abs(post_main$estPct), na.rm = TRUE)
  pval_family <- mean(max_null$max_abs >= main_max_abs, na.rm = TRUE)
  cat(sprintf("\nFamily-wise placebo p-value (max |post effect|): %.4f\n", pval_family))
}

# Plot treated − control difference with symmetric 95% CIs
ggplot(ci_tab, aes(k, estPct)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = loPct, ymax = hiPct), width = 0) +
  labs(
    title = "Event Study: Treated − Control vs 2010Q4, 95% CI",
    x = "Quarters from 2011/Q1", y = "% Difference vs t=-1"
  ) +
  theme_minimal(base_size = 12)

cat("Post-treatment quarters with significantly positive Treated−Control difference (symmetric 95%):\n")
print(ci_tab$k[ci_tab$sigpos])






# ============================================================
# 6) Decompose Treated−Control into treated-country contributions
#     → stacked bar chart (same axes as event-study plot)
# ============================================================

# Helper: Treated and Control aggregate log changes by quarter (aligned)
g_T <- gq_rel |> dplyr::filter(ever == 1) |> dplyr::select(k, ln_rel) |> dplyr::rename(ln_T = ln_rel)
g_C <- gq_rel |> dplyr::filter(ever == 0) |> dplyr::select(k, ln_rel) |> dplyr::rename(ln_C = ln_rel)
agg_path <- dplyr::left_join(g_T, g_C, by = "k") |>
  dplyr::mutate(diff_log = ln_T - ln_C,
                diff_pct = 100*(exp(diff_log) - 1))

# Treated countries with a valid baseline (k = -1)
treated_base <- country_q |>
  dplyr::filter(ever == 1, k == -1, value > 0) |>
  dplyr::select(country, base = value)

stopifnot(nrow(treated_base) > 0)

# Baseline treated-group totals and shares
T_base_tot <- sum(treated_base$base, na.rm = TRUE)
treated_base <- treated_base |>
  dplyr::mutate(w_base = base / T_base_tot)

# Build country × quarter values for treated countries only
treated_panel <- country_q |>
  dplyr::filter(ever == 1) |>
  dplyr::select(country, k, value) |>
  dplyr::left_join(treated_base, by = "country")  # brings base & w_base

# Function to compute treated-country contributions for one quarter k
contrib_one_k <- function(kk) {
  # countries with baseline AND an observation this quarter
  dfk <- treated_panel |>
    dplyr::filter(k == kk, !is.na(base), base > 0, !is.na(value), value > 0)
  
  if (nrow(dfk) == 0) return(NULL)
  
  # Individual log change vs own baseline
  dfk <- dfk |>
    dplyr::mutate(ln_rel_i = log(value) - log(base))
  
  # Shares in treated total at baseline and at k
  T_k_tot <- sum(dfk$value, na.rm = TRUE)
  if (!is.finite(T_k_tot) || T_k_tot <= 0) return(NULL)
  
  dfk <- dfk |>
    dplyr::mutate(s_k   = value / T_k_tot,
                  w_avg = 0.5*(w_base + s_k))   # Törnqvist weights
  
  # Raw Törnqvist contribution (treated side)
  dfk <- dfk |>
    dplyr::mutate(c_T_raw = w_avg * ln_rel_i)
  
  # Scale to match aggregate treated log change exactly
  ln_T_k <- agg_path$ln_T[agg_path$k == kk]
  sum_raw <- sum(dfk$c_T_raw, na.rm = TRUE)
  scale_fac <- ifelse(is.finite(sum_raw) && abs(sum_raw) > 0, ln_T_k / sum_raw, 0)
  
  dfk <- dfk |>
    dplyr::mutate(c_T = scale_fac * c_T_raw)
  
  # Allocate the control group's log change across treated countries by treated baseline shares
  ln_C_k <- agg_path$ln_C[agg_path$k == kk]
  dfk <- dfk |>
    dplyr::mutate(c_C_alloc = w_base * ln_C_k)
  
  # Country contribution to Treated−Control (in log points)
  dfk <- dfk |>
    dplyr::mutate(c_diff_log = c_T - c_C_alloc) |>
    dplyr::select(country, k, c_diff_log)
  
  dfk
}

# Compute contributions for all k in your window
ks <- sort(unique(agg_path$k))
contrib_log <- dplyr::bind_rows(lapply(ks, contrib_one_k))

# If any k returned NULL (no data), drop them from plotting to keep things clean
contrib_log <- contrib_log |> dplyr::filter(!is.na(c_diff_log))

# Sanity check: summed contributions = aggregate diff (log)
check_sum <- contrib_log |>
  dplyr::group_by(k) |>
  dplyr::summarise(sum_c = sum(c_diff_log, na.rm = TRUE), .groups = "drop") |>
  dplyr::left_join(agg_path, by = "k") |>
  dplyr::mutate(err = sum_c - diff_log)

cat("Max absolute mismatch between sum of country contributions and aggregate diff (log pts): ",
    max(abs(check_sum$err), na.rm = TRUE), "\n")

# Map contributions to percent so stacked bars match the event-study percent line
# (distribute total percent effect proportionally to log contributions)
contrib_pct <- contrib_log |>
  dplyr::left_join(agg_path |> dplyr::select(k, diff_log, diff_pct), by = "k") |>
  dplyr::mutate(share = dplyr::if_else(abs(diff_log) > 0, c_diff_log / diff_log, 0),
                contr_pct = share * diff_pct)

# Optional: order legend by total absolute contribution across quarters
country_order <- contrib_pct |>
  dplyr::group_by(country) |>
  dplyr::summarise(tot = sum(abs(contr_pct), na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(tot)) |>
  dplyr::pull(country)

contrib_pct$country <- factor(contrib_pct$country, levels = country_order)

# ---- Plot: stacked contributions (percent), same axes as your event study ----
p_contrib <- ggplot(contrib_pct, aes(x = k, y = contr_pct, fill = country)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_col(width = 0.9, alpha = 0.9) +
  labs(
    title = "Treated-country contributions to Treated − Control",
    subtitle = "Törnqvist contributions (treated), control change allocated by treated baseline shares",
    x = "Quarters from 2011Q1",
    y = "Difference vs t=-1",
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.8, "cm")
  )

print(p_contrib)







# ========= Leave-one-out (treated countries) =========
treated_ctrs <- sort(unique(country_q$country[country_q$ever == 1]))
loo_mat <- matrix(NA_real_, nrow = nrow(orig_diff), ncol = length(treated_ctrs),
                  dimnames = list(orig_diff$k, treated_ctrs))

for (j in seq_along(treated_ctrs)) {
  drop_c <- treated_ctrs[j]
  d_b <- compute_diff_series(country_q |> dplyr::filter(country != drop_c))
  if (!is.null(d_b)) {
    # use percent mapping with your main b_log to keep axis comparable
    loo_mat[, j] <- 100*(exp(d_b$diff[match(orig_diff$k, d_b$k)]) - 1)
  }
}

loo_df <- data.frame(
  k    = orig_diff$k,
  est  = 100*(exp(orig_diff$diff) - 1),
  lo   = apply(loo_mat, 1, function(x) min(x, na.rm = TRUE)),
  hi   = apply(loo_mat, 1, function(x) max(x, na.rm = TRUE))
)

ggplot(loo_df, aes(k, est)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) +
  geom_line(size = 1) +
  labs(title = "Robustness: Leave-one-out Jackknife",
       subtitle = "Grey band = min/max path when excluding one treated country at a time",
       x = "Quarters from 2011Q1", y = "% Difference vs t=-1") +
  theme_minimal(base_size = 12)














# ========= Placebo shocks =========
placebo_qidx <- c(2009L*4 + 1L - 1L, 2010L*4 + 1L - 1L)

compute_diff_series_with_shock <- function(df_country_q, shock_idx) {
  gq <- df_country_q |>
    dplyr::group_by(ever, qidx) |>
    dplyr::summarise(total = sum(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(k = qidx - shock_idx) |>
    dplyr::filter(k >= -W, k <= W)
  
  base2 <- gq |> dplyr::filter(k == -1) |> dplyr::select(ever, base_total = total)
  if (nrow(base2) < 2 || any(base2$base_total <= 0)) return(NULL)
  
  rel <- gq |>
    dplyr::left_join(base2, by = "ever") |>
    dplyr::mutate(ln_rel = log(total) - log(base_total)) |>
    dplyr::select(ever, k, ln_rel) |>
    tidyr::pivot_wider(names_from = ever, values_from = ln_rel, names_prefix = "ever_")
  if (!all(c("ever_0","ever_1") %in% names(rel))) return(NULL)
  rel <- rel[order(rel$k), ]
  data.frame(k = rel$k, estPct = 100*(exp(rel$ever_1 - rel$ever_0) - 1))
}

placebos <- lapply(placebo_qidx, function(pq) compute_diff_series_with_shock(country_q, pq))
names(placebos) <- c("Placebo 2009Q1","Placebo 2010Q1")

# Overlay placebo lines on the main symmetric-CI plot values (ci_tab from above)
base_df <- data.frame(k = ci_tab$k, main = ci_tab$estPct)
p_placebo <- ggplot(base_df, aes(k, main)) +
  geom_hline(yintercept = 0, linetype = 2) + geom_vline(xintercept = 0, linetype = 3) +
  geom_line(size = 1.1, color = "black") +
  lapply(seq_along(placebos), function(i) {
    if (!is.null(placebos[[i]])) {
      geom_line(data = placebos[[i]], aes(k, estPct), inherit.aes = FALSE,
                size = 0.8, alpha = 0.7, color = "grey50")
    }
  }) +
  labs(title = "Robustness: Placebo shocks (pre-period)",
       subtitle = "Gray lines = placebo events; black = main estimate",
       x = "Quarters (k)", y = "Difference vs baseline (percent)") +
  theme_minimal(base_size = 12)
print(p_placebo)










# =========================
# MORE PLACEBOS (auto)
# =========================

# Reusable: recompute diff series for an arbitrary placebo shock index
compute_diff_series_with_shock <- function(df_country_q, shock_idx) {
  gq <- df_country_q |>
    dplyr::group_by(ever, qidx) |>
    dplyr::summarise(total = sum(value, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(k = qidx - shock_idx) |>
    dplyr::filter(k >= -W, k <= W)
  
  base2 <- gq |> dplyr::filter(k == -1) |> dplyr::select(ever, base_total = total)
  if (nrow(base2) < 2 || any(base2$base_total <= 0)) return(NULL)
  
  rel <- gq |>
    dplyr::left_join(base2, by = "ever") |>
    dplyr::mutate(ln_rel = log(total) - log(base_total)) |>
    dplyr::select(ever, k, ln_rel) |>
    tidyr::pivot_wider(names_from = ever, values_from = ln_rel, names_prefix = "ever_")
  if (!all(c("ever_0","ever_1") %in% names(rel))) return(NULL)
  rel <- rel[order(rel$k), ]
  data.frame(k = rel$k, estPct = 100*(exp(rel$ever_1 - rel$ever_0) - 1))
}

# Candidates = ALL pre-event quarters we have in the data
pre_candidates <- sort(unique(country_q$qidx[country_q$qidx < shock_qidx]))

# Build all valid placebo paths (keeps only those with a valid baseline @ k=-1)
placebos_df <- dplyr::bind_rows(lapply(pre_candidates, function(pq) {
  ds <- compute_diff_series_with_shock(country_q, pq)
  if (is.null(ds)) return(NULL)
  ds$shock_qidx <- pq
  ds
}))

# If nothing passed, bail early
if (nrow(placebos_df) == 0) stop("No valid placebo shocks were found. Try shrinking W or check data coverage.")

cat("Placebo shocks built:", length(unique(placebos_df$shock_qidx)), "\n")

# 95% placebo envelope across all placebos (by k)
env <- placebos_df |>
  dplyr::group_by(k) |>
  dplyr::summarise(
    lo = quantile(estPct, 0.025, na.rm = TRUE),
    hi = quantile(estPct, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Main series we already computed above (from ci_tab)
main_df <- data.frame(k = ci_tab$k, estPct = ci_tab$estPct)

# ----- Plot: spaghetti + envelope + main -----
p_placebos_many <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_ribbon(data = env, aes(k, ymin = lo, ymax = hi), fill = "grey85", alpha = 0.6) +
  geom_line(data = placebos_df, aes(k, estPct, group = shock_qidx),
            color = "grey60", alpha = 0.5, linewidth = 0.5) +
  geom_line(data = main_df, aes(k, estPct), color = "black", linewidth = 1.1) +
  labs(
    title = "Placebo shocks: Distribution vs. Main Event Study",
    subtitle = "Grey lines are valid placebos (pre-event); shaded area is 95% placebo envelope",
    x = "Quarters from placebo", y = "% Difference vs 2010/Q4"
  ) +
  theme_minimal(base_size = 12)
print(p_placebos_many)

# ----- Per-k empirical p-values (post-period only) -----
post_main <- main_df |> dplyr::filter(k >= 0)
keep_ids <- head(sort(unique(placebos_df$shock_qidx)), 40)   # show only first 40
placebos_plot <- placebos_df |> dplyr::filter(shock_qidx %in% keep_ids)
# then use placebos_plot in the geom_line for the spaghetti

pvals_k <- placebo_post |>
  dplyr::group_by(k) |>
  dplyr::summarise(
    pval = mean(abs(estPct) >= abs(post_main$estPct[match(k, post_main$k)]), na.rm = TRUE),
    .groups = "drop"
  )

cat("\nEmpirical placebo p-values by k (post only):\n")
print(pvals_k)

# ----- Family-wise (max-stat) p-value across all post-k -----
# For each placebo shock, take the max |effect| over post-k; compare to main max |effect|
max_null <- placebo_post |>
  dplyr::group_by(shock_qidx) |>
  dplyr::summarise(max_abs = max(abs(estPct), na.rm = TRUE), .groups = "drop")

main_max_abs <- max(abs(post_main$estPct), na.rm = TRUE)
pval_family <- mean(max_null$max_abs >= main_max_abs, na.rm = TRUE)

cat(sprintf("\nFamily-wise (max |effect| over post-k) placebo p-value: %.4f\n", pval_family))

















# ============================================================
# B) Decompose Treated−Control into reporting-country contributions
#     (mirror of the treated-country Törnqvist decomposition)
# ============================================================

# --- 0) Identify the reporting-country column (include 're_country') ---
report_col_candidates <- c(
  "re_country","reporting_country","reporter","reporter_country",
  "reportername","report","rp","rp_country","destination",
  "bank_center","reportingcountry"
)
report_col <- intersect(report_col_candidates, names(dat))[1]
if (is.na(report_col)) {
  stop(paste0(
    "Couldn't find a reporting-country column. Rename the appropriate column to one of: ",
    paste(report_col_candidates, collapse = ", ")
  ))
}
report_sym <- rlang::sym(report_col)

# --- 1) Build k (event time) on the raw rows and keep ever-treated flag ---
dat_k <- dat |>
  dplyr::left_join(ever_tab, by = "country") |>
  dplyr::mutate(k = qidx - shock_qidx)

# --- 2) Baseline (k = -1) reporting-country totals within treated origins
rp_base <- dat_k |>
  dplyr::filter(ever == 1L, k == -1, value > 0) |>
  dplyr::group_by(reporting = !!report_sym) |>
  dplyr::summarise(base = sum(value, na.rm = TRUE), .groups = "drop") |>
  dplyr::filter(is.finite(base) & base > 0)

stopifnot(nrow(rp_base) > 0)

T_base_tot_rp <- sum(rp_base$base, na.rm = TRUE)
rp_base <- rp_base |> dplyr::mutate(w_base = base / T_base_tot_rp)

# --- 3) Reporting-country panel over k for treated-origin flows ---
rp_panel <- dat_k |>
  dplyr::filter(ever == 1L) |>
  dplyr::group_by(reporting = !!report_sym, k) |>
  dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

# --- 4) Quarter-k contributions function (Törnqvist) ---
contrib_rp_one_k <- function(kk) {
  dfk <- rp_panel |>
    dplyr::filter(k == kk) |>
    dplyr::left_join(rp_base, by = "reporting") |>
    dplyr::filter(!is.na(base), base > 0, !is.na(value), value > 0)
  if (nrow(dfk) == 0) return(NULL)
  
  # Individual log change vs own baseline
  dfk <- dfk |>
    dplyr::mutate(ln_rel_i = log(value) - log(base))
  
  # Shares in treated total at baseline and at k
  T_k_tot <- sum(dfk$value, na.rm = TRUE)
  if (!is.finite(T_k_tot) || T_k_tot <= 0) return(NULL)
  
  dfk <- dfk |>
    dplyr::mutate(
      s_k   = value / T_k_tot,
      w_avg = 0.5*(w_base + s_k),      # Törnqvist weights
      c_T_raw = w_avg * ln_rel_i
    )
  
  # Scale treated-side contributions to match aggregate treated log change exactly
  ln_T_k <- agg_path$ln_T[agg_path$k == kk]
  sum_raw <- sum(dfk$c_T_raw, na.rm = TRUE)
  scale_fac <- ifelse(is.finite(sum_raw) && abs(sum_raw) > 0, ln_T_k / sum_raw, 0)
  
  dfk <- dfk |>
    dplyr::mutate(c_T = scale_fac * c_T_raw)
  
  # Allocate the control group's log change across reporting countries by baseline shares
  ln_C_k <- agg_path$ln_C[agg_path$k == kk]
  dfk <- dfk |>
    dplyr::mutate(c_C_alloc = w_base * ln_C_k)
  
  # Reporting-country contribution to Treated − Control (in log points)
  dfk |>
    dplyr::mutate(c_diff_log = c_T - c_C_alloc) |>
    dplyr::select(reporting, k, c_diff_log)
}

# --- 5) Compute contributions for all k in your window ---
ks <- sort(unique(agg_path$k))
rp_contrib_log <- dplyr::bind_rows(lapply(ks, contrib_rp_one_k)) |>
  dplyr::filter(!is.na(c_diff_log))

# Sanity check: summed contributions = aggregate diff (log)
rp_check <- rp_contrib_log |>
  dplyr::group_by(k) |>
  dplyr::summarise(sum_c = sum(c_diff_log, na.rm = TRUE), .groups = "drop") |>
  dplyr::left_join(agg_path, by = "k") |>
  dplyr::mutate(err = sum_c - diff_log)

cat("Reporting-country: max |mismatch| between sum of contributions and aggregate diff (log pts): ",
    max(abs(rp_check$err), na.rm = TRUE), "\n")

# --- 6) Map contributions to percent so stacked bars match the event-study percent line ---
rp_contrib_pct <- rp_contrib_log |>
  dplyr::left_join(agg_path |> dplyr::select(k, diff_log, diff_pct), by = "k") |>
  dplyr::mutate(
    share = dplyr::if_else(abs(diff_log) > 0, c_diff_log / diff_log, 0),
    contr_pct = share * diff_pct
  )

# Order legend by total absolute contribution across quarters
rp_order <- rp_contrib_pct |>
  dplyr::group_by(reporting) |>
  dplyr::summarise(tot = sum(abs(contr_pct), na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(tot)) |>
  dplyr::pull(reporting)

rp_contrib_pct$reporting <- factor(rp_contrib_pct$reporting, levels = rp_order)

# --- 7) Plot: stacked contributions (percent), same axes as your event study ---
p_contrib_reporting <- ggplot(rp_contrib_pct, aes(x = k, y = contr_pct, fill = reporting)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_col(width = 0.9, alpha = 0.9) +
  labs(
    title = "Reporting-country contributions to Treated − Control",
    subtitle = "Törnqvist contributions (treated origins), control change allocated by baseline shares",
    x = "Quarters from 2011Q1",
    y = "Difference vs t=-1",
    fill = "Reporting country"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.key.height = unit(0.4, "cm"),
    legend.key.width  = unit(0.8, "cm")
  )

print(p_contrib_reporting)
