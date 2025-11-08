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
# 1a) EU membership in 2011 -> Eur = 1 (else 0)
# -------------------------
eu_2011 <- c(
  "Bulgaria","Czechia","Denmark","Estonia","Finland","France","Germany","Greece",
  "Hungary","Italy","Latvia","Lithuania","Poland","Portugal","Romania",
  "Slovakia","Slovenia","Spain","Sweden","United Kingdom"
)
dat <- dat |>
  mutate(
    th    = as.integer(th),
    treat = as.integer(treat),
    value = as.numeric(value),
    Eur   = as.integer(country %in% eu_2011)
  )

# -------------------------
# 2) Build country × TH × quarter totals (keep TH separate)
# -------------------------
cth_q <- dat |>
  group_by(country, th, qidx) |>
  summarise(
    value = sum(value, na.rm = TRUE),
    treat = as.integer(any(treat == 1L)),
    Eur   = first(Eur),
    .groups = "drop"
  )

# Ever-treated (by country, consistent with earlier scripts)
ever_tab <- cth_q |>
  group_by(country) |>
  summarise(ever = as.integer(any(treat == 1L)), .groups = "drop")

cth_q <- cth_q |>
  left_join(ever_tab, by = "country") |>
  mutate(k = qidx - shock_qidx)

# -------------------------
# Restrict to non-EU
# -------------------------
cth_q_neu <- cth_q |> filter(Eur == 0)
cat("Non-EU countries in data:", length(unique(cth_q_neu$country)), "\n")

# -------------------------
# Helper: compute treated−control diff series within a TH group
# -------------------------
compute_diff_series_TH <- function(df_cth, th_value) {
  gq <- df_cth |>
    filter(th == th_value) |>
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

# -------------------------
# Plot+bootstrap for one TH
# -------------------------
run_event_study_TH <- function(df_cth, th_value, B = 500) {
  
  # Two-line aggregate-first plot (within this TH)
  gq_sum <- df_cth |>
    filter(th == th_value) |>
    group_by(ever, qidx) |>
    summarise(total = sum(value, na.rm = TRUE), .groups = "drop") |>
    mutate(k = qidx - shock_qidx) |>
    filter(k >= -W, k <= W)
  
  base <- gq_sum |> filter(k == -1) |> select(ever, base_total = total)
  if (nrow(base) < 2 || any(base$base_total <= 0)) {
    stop(sprintf("TH==%d: missing/zero baseline at k=-1 for a group.", th_value))
  }
  
  gq_rel <- gq_sum |>
    left_join(base, by = "ever") |>
    mutate(
      ln_rel = log(total) - log(base_total),
      pct    = 100*(exp(ln_rel) - 1),
      group  = ifelse(ever == 1, "Treated", "Control")
    )
  
  p_lines <- ggplot(gq_rel, aes(k, pct, color = group)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_line(size = 1) + geom_point(size = 1.6) +
    scale_color_manual(values = c("Treated"="#1f9ac9","Control"="#fb8161")) +
    labs(
      title = sprintf("TH == %d • Non-EU • %% change vs group’s 2010Q4", th_value),
      subtitle = sprintf("Aggregate-first; window ±%dq; t(0)=2011Q1", W),
      x = "Quarters relative to shock", y = "% change from 2010Q4", color = ""
    ) +
    theme_minimal(base_size = 12)
  
  # Treated−Control diff with bootstrap CIs (within this TH)
  orig <- compute_diff_series_TH(df_cth, th_value)
  if (is.null(orig)) stop(sprintf("TH==%d: could not compute difference series.", th_value))
  
  countries_T <- unique(df_cth$country[df_cth$th == th_value & df_cth$ever == 1])
  countries_C <- unique(df_cth$country[df_cth$th == th_value & df_cth$ever == 0])
  nT <- length(countries_T); nC <- length(countries_C)
  cat(sprintf("TH==%d • Non-EU countries: Treated=%d, Control=%d\n", th_value, nT, nC))
  
  if (nT == 0 || nC == 0) stop(sprintf("TH==%d: need both treated and control countries.", th_value))
  
  boot_mat <- matrix(NA_real_, nrow = nrow(orig), ncol = B, dimnames = list(orig$k, NULL))
  
  for (b in 1:B) {
    samp_T <- sample(countries_T, nT, replace = TRUE)
    samp_C <- sample(countries_C, nC, replace = TRUE)
    df_b <- dplyr::bind_rows(
      df_cth %>% filter(th == th_value, country %in% samp_T),
      df_cth %>% filter(th == th_value, country %in% samp_C)
    )
    d_b <- compute_diff_series_TH(df_b, th_value)
    if (!is.null(d_b)) {
      idx <- match(orig$k, d_b$k)
      boot_mat[, b] <- d_b$diff[idx]
    }
  }
  boot_mat <- boot_mat[, colSums(is.na(boot_mat)) == 0, drop = FALSE]
  # --- SYMMETRIC CIs via bootstrap SEs + delta method (percent scale) ---
  z <- qnorm(0.975)
  
  # point estimate in log space
  b_log <- orig$diff                      # log-points: ln T − ln C (vs own baselines)
  
  # bootstrap SD in log space
  se_log <- apply(boot_mat, 1, sd, na.rm = TRUE)
  
  # transform to percent via delta method (symmetric in percent space)
  estPct <- 100 * (exp(b_log) - 1)
  sePct  <- 100 * exp(b_log) * se_log
  
  ci_tab <- data.frame(
    k      = orig$k,
    b_log  = b_log,
    se_log = se_log,
    estPct = estPct,
    sePct  = sePct,
    loPct  = estPct - z * sePct,
    hiPct  = estPct + z * sePct
  )
  ci_tab$sigpos <- (ci_tab$k >= 0) & (ci_tab$loPct > 0)
  
  # --- plot with symmetric bands ---
  p_diff <- ggplot(ci_tab, aes(k, estPct)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = loPct, ymax = hiPct), width = 0) +
    labs(
      title = sprintf("TH == %d • Non-EU • Treated − Control vs 2010Q4", th_value),
      subtitle = sprintf("Symmetric 95%% CIs: bootstrap SE + delta method (B=%d)", ncol(boot_mat)),
      x = "Quarters from 2011Q1 (k)", y = "Difference vs baseline (percent)"
    ) +
    theme_minimal(base_size = 12)
  
  
  list(two_line_plot = p_lines, diff_plot = p_diff,
       sigpos_k = ci_tab$k[ci_tab$sigpos], ci_table = ci_tab)
}

# -------------------------
# Run for TH==1 and TH==0 on the non-EU sample
# -------------------------
res_TH1 <- run_event_study_TH(cth_q_neu, th_value = 1, B = 500)
res_TH0 <- run_event_study_TH(cth_q_neu, th_value = 0, B = 500)

# Show plots
print(res_TH1$two_line_plot); print(res_TH1$diff_plot)
print(res_TH0$two_line_plot); print(res_TH0$diff_plot)

cat("\nTH==1 (Non-EU) post-treatment k significantly > 0 (95%):\n")
print(res_TH1$sigpos_k)
cat("\nTH==0 (Non-EU) post-treatment k significantly > 0 (95%):\n")
print(res_TH0$sigpos_k)
