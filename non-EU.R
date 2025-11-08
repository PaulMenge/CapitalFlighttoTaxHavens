suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(lubridate); library(zoo)
  library(ggplot2); library(tidyr)
})

# -------------------------
# 0) Parameters
# -------------------------
W <- 10
shock_qidx <- 2011L*4 + 1L - 1L   # 2011Q1
set.seed(123)

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
# 1a) Mark EU membership around the crisis year (2011)
# -------------------------
eu_2011 <- c(
  "Bulgaria","Czechia","Denmark","Estonia","Finland","France","Germany","Greece",
  "Hungary","Italy","Latvia","Lithuania","Poland","Portugal","Romania",
  "Slovakia","Slovenia","Spain","Sweden","United Kingdom"
)

dat <- dat |>
  mutate(Eur = as.integer(country %in% eu_2011))  # 1 = EU member in 2011

# quick sanity
eu_present     <- sort(unique(dat$country[dat$Eur == 1]))
non_eu_present <- sort(unique(dat$country[dat$Eur == 0]))
cat("EU countries found in data (2011 definition):", length(eu_present), "\n")
cat("Non-EU countries found in data:", length(non_eu_present), "\n")

# -------------------------
# 2) Collapse to one obs per (country, quarter), ignoring TH inside country
# -------------------------
country_q <- dat |>
  group_by(country, qidx) |>
  summarise(
    value = sum(value, na.rm = TRUE),      # sum within country×quarter across TH cells
    treat = as.integer(any(treat == 1L)),  # country treated that quarter?
    .groups = "drop"
  ) |>
  left_join(dat |> distinct(country, Eur), by = "country")

# Ever-treated flag per country
ever_tab <- country_q |>
  group_by(country) |>
  summarise(ever = as.integer(any(treat == 1L)), .groups = "drop")

country_q <- country_q |>
  left_join(ever_tab, by = "country") |>
  mutate(k = qidx - shock_qidx)

# -------------------------
# Restrict to non-EU sample
# -------------------------
country_q_neu <- country_q |> filter(Eur == 0)
cat("Countries used (Eur==0):", length(unique(country_q_neu$country)), "\n")

# -------------------------
# 3) AGGREGATE FIRST: sum across countries within each group×quarter (non-EU only)
# -------------------------
gq_sum <- country_q_neu |>
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
# 4) Two-line plot: % change from each group's baseline (aggregate-first, non-EU)
# -------------------------
ggplot(gq_rel, aes(k, pct, color = group)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_line(size = 1) + geom_point(size = 1.6) +
  scale_color_manual(values = c("Treated"="#1f9ac9","Control"="#fb8161")) +
  labs(
    title = "Non-EU sample • Aggregate totals: % change vs group’s 2010Q4",
    subtitle = sprintf("Aggregate-first normalization; window ±%dq; t(0)=2011Q1", W),
    x = "Quarters relative to shock",
    y = "% change from 2010Q4",
    color = ""
  ) +
  theme_minimal(base_size = 12)

# -------------------------
# 5) Event-time difference (treated − control) – symmetric 95% CIs
#    d_k = [ln Tot_T(k) − ln Tot_T(-1)] − [ln Tot_C(k) − ln Tot_C(-1)]
# -------------------------
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
    pivot_wider(names_from = ever, values_from = ln_rel, names_prefix = "ever_")
  
  if (!all(c("ever_0","ever_1") %in% names(rel))) return(NULL)
  rel <- rel[order(rel$k), ]
  transform(rel, diff = ever_1 - ever_0)[, c("k","diff")]
}

# Point estimate (non-EU)
orig <- compute_diff_series(country_q_neu)
stopifnot(!is.null(orig))
orig_diff <- orig

# Bootstrap by resampling countries within each group (non-EU)
B <- 500
countries_T <- unique(country_q_neu$country[country_q_neu$ever == 1])
countries_C <- unique(country_q_neu$country[country_q_neu$ever == 0])
nT <- length(countries_T); nC <- length(countries_C)
cat(sprintf("Non-EU countries: Treated=%d, Control=%d\n", nT, nC))

boot_mat <- matrix(NA_real_, nrow = nrow(orig_diff), ncol = B)
rownames(boot_mat) <- orig_diff$k

for (b in 1:B) {
  samp_T <- sample(countries_T, nT, replace = TRUE)
  samp_C <- sample(countries_C, nC, replace = TRUE)
  df_b <- dplyr::bind_rows(
    country_q_neu %>% filter(country %in% samp_T),
    country_q_neu %>% filter(country %in% samp_C)
  )
  d_b <- compute_diff_series(df_b)
  if (!is.null(d_b)) {
    idx <- match(orig_diff$k, d_b$k)   # align by k
    boot_mat[, b] <- d_b$diff[idx]
  }
}
boot_mat <- boot_mat[, colSums(is.na(boot_mat)) == 0, drop = FALSE]
if (ncol(boot_mat) == 0) stop("All bootstrap resamples failed (NA columns). Check data coverage around k = -1.")

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

# Plot treated − control difference with symmetric 95% CIs (non-EU)
ggplot(ci_tab, aes(k, estPct)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = loPct, ymax = hiPct), width = 0) +
  labs(
    title = "Non-EU sample Event Study: Treated − Control vs 2010Q4, 95% CI",
    x = "Quarters from 2011/Q1", y = "% Difference vs t=-1"
  ) +
  theme_minimal(base_size = 12)

cat("Post-treatment quarters with significantly positive Treated−Control difference (symmetric 95%) — Non-EU only:\n")
print(ci_tab$k[ci_tab$sigpos])










# ============================================================
# A) Treated-country Törnqvist contributions (non-EU origins)
#     → stacked bar chart (same axes as event-study plot)
# ============================================================

# Aggregate treated vs control paths (non-EU)
g_T <- gq_rel %>% dplyr::filter(ever == 1) %>%
  dplyr::select(k, ln_T = ln_rel)
g_C <- gq_rel %>% dplyr::filter(ever == 0) %>%
  dplyr::select(k, ln_C = ln_rel)

agg_path <- dplyr::left_join(g_T, g_C, by = "k") %>%
  dplyr::mutate(
    diff_log = ln_T - ln_C,
    diff_pct = 100 * (exp(diff_log) - 1)
  )

# Treated countries with a valid baseline (k = -1)
treated_base <- country_q_neu %>%
  dplyr::filter(ever == 1, k == -1, value > 0) %>%
  dplyr::select(country, base = value)

stopifnot(nrow(treated_base) > 0)
T_base_tot <- sum(treated_base$base, na.rm = TRUE)
treated_base <- treated_base %>%
  dplyr::mutate(w_base = base / T_base_tot)

# Treated-country panel (values over k), non-EU only
treated_panel <- country_q_neu %>%
  dplyr::filter(ever == 1) %>%
  dplyr::select(country, k, value) %>%
  dplyr::left_join(treated_base, by = "country")

# Quarter-k contributions function (Törnqvist on treated side, allocate control by baseline shares)
contrib_one_k <- function(kk) {
  dfk <- treated_panel %>%
    dplyr::filter(k == kk, !is.na(base), base > 0, !is.na(value), value > 0)
  if (nrow(dfk) == 0) return(NULL)
  
  # Individual log change vs own baseline
  dfk <- dfk %>% dplyr::mutate(ln_rel_i = log(value) - log(base))
  
  # Shares in treated total at baseline and at k
  T_k_tot <- sum(dfk$value, na.rm = TRUE)
  if (!is.finite(T_k_tot) || T_k_tot <= 0) return(NULL)
  
  dfk <- dfk %>%
    dplyr::mutate(
      s_k   = value / T_k_tot,
      w_avg = 0.5 * (w_base + s_k),
      c_T_raw = w_avg * ln_rel_i
    )
  
  # Scale treated-side contributions to match aggregate treated log change exactly
  ln_T_k <- agg_path$ln_T[agg_path$k == kk]
  sum_raw <- sum(dfk$c_T_raw, na.rm = TRUE)
  scale_fac <- ifelse(is.finite(sum_raw) && abs(sum_raw) > 0, ln_T_k / sum_raw, 0)
  dfk <- dfk %>% dplyr::mutate(c_T = scale_fac * c_T_raw)
  
  # Allocate the control-group log change by treated baseline shares
  ln_C_k <- agg_path$ln_C[agg_path$k == kk]
  dfk <- dfk %>% dplyr::mutate(c_C_alloc = w_base * ln_C_k)
  
  # Contribution to Treated − Control (log points)
  dfk %>%
    dplyr::mutate(c_diff_log = c_T - c_C_alloc) %>%
    dplyr::select(country, k, c_diff_log)
}

# Compute contributions across quarters
ks <- sort(unique(agg_path$k))
contrib_log <- dplyr::bind_rows(lapply(ks, contrib_one_k)) %>%
  dplyr::filter(!is.na(c_diff_log))

# Sanity: summed contributions match aggregate diff (log)
check_sum <- contrib_log %>%
  dplyr::group_by(k) %>%
  dplyr::summarise(sum_c = sum(c_diff_log, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(agg_path, by = "k") %>%
  dplyr::mutate(err = sum_c - diff_log)
cat("Non-EU treated-country: max |mismatch| (log pts) = ",
    max(abs(check_sum$err), na.rm = TRUE), "\n")

# Map to percent so stacked bars match event-study percent line
contrib_pct <- contrib_log %>%
  dplyr::left_join(agg_path %>% dplyr::select(k, diff_log, diff_pct), by = "k") %>%
  dplyr::mutate(
    share     = dplyr::if_else(abs(diff_log) > 0, c_diff_log / diff_log, 0),
    contr_pct = share * diff_pct
  )

# Order legend by total absolute contribution
country_order <- contrib_pct %>%
  dplyr::group_by(country) %>%
  dplyr::summarise(tot = sum(abs(contr_pct), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(tot)) %>%
  dplyr::pull(country)
contrib_pct$country <- factor(contrib_pct$country, levels = country_order)

# Plot: stacked treated-country contributions (percent)
p_contrib_countries <- ggplot(contrib_pct, aes(x = k, y = contr_pct, fill = country)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_col(width = 0.9, alpha = 0.9) +
  labs(
    title = "Treated-country Contributions to Treated − Control",
    subtitle = "Törnqvist contributions, control change allocated by treated baseline shares, excl. EU controls",
    x = "Quarters from 2011Q1",
    y = "% Difference vs t = −1"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "bottom",
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width  = grid::unit(0.8, "cm")
  )
print(p_contrib_countries)


# ============================================================
# B) Reporting-country contributions (non-EU origins)
#     mirror of treated-country decomposition using re_country
# ============================================================

# Identify reporting-country column
report_col <- if ("re_country" %in% names(dat)) "re_country" else {
  cand <- intersect(
    c("re_country","reporting_country","reporter","reporter_country",
      "reportername","report","rp","rp_country","destination",
      "bank_center","reportingcountry"),
    names(dat)
  )
  if (length(cand) == 0) stop("No reporting-country column found."); cand[1]
}
report_sym <- rlang::sym(report_col)

# Keep only non-EU origins, add ever flag, and event time k
ever_tab_neu <- country_q_neu %>% dplyr::distinct(country, ever)
dat_k_neu <- dat %>%
  dplyr::filter(Eur == 0) %>%
  dplyr::left_join(ever_tab_neu, by = "country") %>%
  dplyr::mutate(k = qidx - shock_qidx)

# Baseline (k = −1) reporting-country totals within treated origins
rp_base <- dat_k_neu %>%
  dplyr::filter(ever == 1L, k == -1, value > 0) %>%
  dplyr::group_by(reporting = !!report_sym) %>%
  dplyr::summarise(base = sum(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(is.finite(base) & base > 0)
stopifnot(nrow(rp_base) > 0)

T_base_tot_rp <- sum(rp_base$base, na.rm = TRUE)
rp_base <- rp_base %>% dplyr::mutate(w_base = base / T_base_tot_rp)

# Reporting-country panel over k for treated origins (non-EU)
rp_panel <- dat_k_neu %>%
  dplyr::filter(ever == 1L) %>%
  dplyr::group_by(reporting = !!report_sym, k) %>%
  dplyr::summarise(value = sum(value, na.rm = TRUE), .groups = "drop")

# Quarter-k contributions (reporting side)
contrib_rp_one_k <- function(kk) {
  dfk <- rp_panel %>%
    dplyr::filter(k == kk) %>%
    dplyr::left_join(rp_base, by = "reporting") %>%
    dplyr::filter(!is.na(base), base > 0, !is.na(value), value > 0)
  if (nrow(dfk) == 0) return(NULL)
  
  dfk <- dfk %>% dplyr::mutate(ln_rel_i = log(value) - log(base))
  
  T_k_tot <- sum(dfk$value, na.rm = TRUE)
  if (!is.finite(T_k_tot) || T_k_tot <= 0) return(NULL)
  
  dfk <- dfk %>%
    dplyr::mutate(
      s_k    = value / T_k_tot,
      w_avg  = 0.5 * (w_base + s_k),
      c_T_raw = w_avg * ln_rel_i
    )
  
  # Scale to match aggregate treated log change exactly (non-EU agg_path)
  ln_T_k <- agg_path$ln_T[agg_path$k == kk]
  sum_raw <- sum(dfk$c_T_raw, na.rm = TRUE)
  scale_fac <- ifelse(is.finite(sum_raw) && abs(sum_raw) > 0, ln_T_k / sum_raw, 0)
  dfk <- dfk %>% dplyr::mutate(c_T = scale_fac * c_T_raw)
  
  # Allocate control log change by baseline shares
  ln_C_k <- agg_path$ln_C[agg_path$k == kk]
  dfk <- dfk %>% dplyr::mutate(c_C_alloc = w_base * ln_C_k)
  
  dfk %>%
    dplyr::mutate(c_diff_log = c_T - c_C_alloc) %>%
    dplyr::select(reporting, k, c_diff_log)
}

# Compute reporting-country contributions over k
rp_contrib_log <- dplyr::bind_rows(lapply(ks, contrib_rp_one_k)) %>%
  dplyr::filter(!is.na(c_diff_log))

# Sanity: sum of contributions equals aggregate diff (log)
rp_check <- rp_contrib_log %>%
  dplyr::group_by(k) %>%
  dplyr::summarise(sum_c = sum(c_diff_log, na.rm = TRUE), .groups = "drop") %>%
  dplyr::left_join(agg_path, by = "k") %>%
  dplyr::mutate(err = sum_c - diff_log)
cat("Non-EU reporting-country: max |mismatch| (log pts) = ",
    max(abs(rp_check$err), na.rm = TRUE), "\n")

# Map to percent for stacked bars
rp_contrib_pct <- rp_contrib_log %>%
  dplyr::left_join(agg_path %>% dplyr::select(k, diff_log, diff_pct), by = "k") %>%
  dplyr::mutate(
    share     = dplyr::if_else(abs(diff_log) > 0, c_diff_log / diff_log, 0),
    contr_pct = share * diff_pct
  )

# Legend order by absolute contribution
rp_order <- rp_contrib_pct %>%
  dplyr::group_by(reporting) %>%
  dplyr::summarise(tot = sum(abs(contr_pct), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(tot)) %>%
  dplyr::pull(reporting)
rp_contrib_pct$reporting <- factor(rp_contrib_pct$reporting, levels = rp_order)

# Plot: stacked reporting-country contributions (percent)
p_contrib_reporting <- ggplot(rp_contrib_pct, aes(x = k, y = contr_pct, fill = reporting)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_col(width = 0.9, alpha = 0.9) +
  labs(
    title = "Reporting-country Contributions to Treated − Control",
    subtitle = "Törnqvist contributions, control change allocated by baseline shares, excl. EU controls",
    x = "Quarters from 2011Q1",
    y = "% Difference vs t = −1",
    fill = "Reporting country"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position   = "bottom",
    legend.key.height = grid::unit(0.4, "cm"),
    legend.key.width  = grid::unit(0.8, "cm")
  )
print(p_contrib_reporting)




# Helper: event-study diff for an arbitrary placebo shock index (returns % scale)
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

# Build all valid pre-event placebos in the non-EU sample
pre_candidates <- sort(unique(country_q_neu$qidx[country_q_neu$qidx < shock_qidx]))

placebos_df <- dplyr::bind_rows(lapply(pre_candidates, function(pq) {
  ds <- compute_diff_series_with_shock(country_q_neu, pq)
  if (is.null(ds)) return(NULL)
  ds$shock_qidx <- pq
  ds
}))

if (nrow(placebos_df) == 0) {
  warning("No valid placebo shocks found in non-EU sample. Skipping placebo plots and p-values.")
} else {
  cat("Placebo shocks built (non-EU):", length(unique(placebos_df$shock_qidx)), "\n")
  
  # 95% envelope across placebos (by k)
  env <- placebos_df |>
    dplyr::group_by(k) |>
    dplyr::summarise(
      lo = quantile(estPct, 0.025, na.rm = TRUE),
      hi = quantile(estPct, 0.975, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Main series on % scale (use ci_tab if present; else rebuild from country_q_neu)
  if (exists("ci_tab")) {
    main_df <- data.frame(k = ci_tab$k, estPct = ci_tab$estPct)
  } else {
    orig_tmp <- compute_diff_series(country_q_neu)
    stopifnot(!is.null(orig_tmp))
    main_df <- data.frame(k = orig_tmp$k, estPct = 100*(exp(orig_tmp$diff) - 1))
  }
  
  # Spaghetti subset for readability (first 40 placebos)
  keep_ids <- head(sort(unique(placebos_df$shock_qidx)), 40)
  placebos_plot <- placebos_df |> dplyr::filter(shock_qidx %in% keep_ids)
  
  # ---- Plot: spaghetti + 95% envelope + main line (assign + print) ----
  p_placebo_spaghetti <- ggplot() +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_ribbon(data = env, aes(k, ymin = lo, ymax = hi), fill = "grey85", alpha = 0.6) +
    geom_line(data = placebos_plot, aes(k, estPct, group = shock_qidx),
              color = "grey60", alpha = 0.5, linewidth = 0.5) +
    geom_line(data = main_df, aes(k, estPct), color = "black", linewidth = 1.1) +
    labs(
      title = "Placebo shocks: Distribution vs. Main Event Study",
      subtitle = "Grey lines: 40 placebos (pre-event); shaded area = 95% placebo envelope, non-EU control",
      x = "Quarters from placebo", y = "% Difference vs 2010Q4"
    ) +
    theme_minimal(base_size = 12)
  
  print(p_placebo_spaghetti)
  
  # ---- Empirical p-values by k (post only) ----
  post_main    <- main_df      |> dplyr::filter(k >= 0)
  placebo_post <- placebos_df  |> dplyr::filter(k >= 0)
  
  pvals_k <- placebo_post |>
    dplyr::group_by(k) |>
    dplyr::summarise(
      pval = mean(abs(estPct) >= abs(post_main$estPct[match(k, post_main$k)]), na.rm = TRUE),
      .groups = "drop"
    )
  cat("\nEmpirical placebo p-values by k (post only, non-EU):\n")
  print(pvals_k)
  
  # ---- Family-wise (max-stat) p-value across all post-k ----
  max_null <- placebo_post |>
    dplyr::group_by(shock_qidx) |>
    dplyr::summarise(max_abs = max(abs(estPct), na.rm = TRUE), .groups = "drop")
  main_max_abs <- max(abs(post_main$estPct), na.rm = TRUE)
  pval_family <- mean(max_null$max_abs >= main_max_abs, na.rm = TRUE)
  cat(sprintf("\nFamily-wise placebo p-value (max |post effect|): %.4f\n", pval_family))
}




# =========================
# Jackknife: leave-one-treated-country-out (non-EU)
# =========================

# Treated countries in the non-EU sample
treated_ctrs_neu <- sort(unique(country_q_neu$country[country_q_neu$ever == 1]))
stopifnot(length(treated_ctrs_neu) > 0)

# Main percent series (to keep axes identical)
main_k   <- orig_diff$k
main_pct <- 100 * (exp(orig_diff$diff) - 1)

# Build % paths when excluding one treated country at a time
loo_mat <- matrix(NA_real_, nrow = length(main_k), ncol = length(treated_ctrs_neu),
                  dimnames = list(main_k, treated_ctrs_neu))

for (j in seq_along(treated_ctrs_neu)) {
  drop_c <- treated_ctrs_neu[j]
  d_b <- compute_diff_series(country_q_neu |> dplyr::filter(country != drop_c))
  if (!is.null(d_b)) {
    # map to percent to match main axis
    pct_b <- 100 * (exp(d_b$diff) - 1)
    loo_mat[, j] <- pct_b[match(main_k, d_b$k)]
  }
}

# Min–max envelope across leave-one-out paths
loo_lo <- apply(loo_mat, 1, function(x) min(x, na.rm = TRUE))
loo_hi <- apply(loo_mat, 1, function(x) max(x, na.rm = TRUE))
loo_df <- data.frame(k = main_k, est = main_pct, lo = loo_lo, hi = loo_hi)

# ---- Plot: Jackknife envelope + main line ----
p_jackknife <- ggplot(loo_df, aes(k, est)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80", alpha = 0.6) +
  geom_line(color = "black", linewidth = 1) +
  labs(
    title = "Jackknife (non-EU): Leave-one-treated-country-out",
    subtitle = "Grey band = min/max path when excluding one treated country at a time",
    x = "Quarters from 2011Q1", y = "% Difference vs t = −1"
  ) +
  theme_minimal(base_size = 12)

print(p_jackknife)

# ---- Simple influence table (post period): who widens the band most? ----
post_idx <- which(main_k >= 0)
influence <- vapply(seq_along(treated_ctrs_neu), function(j) {
  x <- loo_mat[, j]
  if (all(is.na(x))) return(NA_real_)
  max(abs(x[post_idx] - main_pct[post_idx]), na.rm = TRUE)
}, numeric(1))
infl_tbl <- data.frame(country = treated_ctrs_neu, max_abs_dev_post = influence) |>
  dplyr::arrange(dplyr::desc(max_abs_dev_post))

cat("\nJackknife influence (top 10 by max |deviation| in post period, % points):\n")
print(utils::head(infl_tbl, 19))
