suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(lubridate); library(zoo)
  library(ggplot2); library(tidyr)
  set.seed(123)
})

# -------------------------
# 0) Params
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
    th    = as.integer(th),
    treat = as.integer(treat),
    qidx  = lubridate::year(date)*4 + lubridate::quarter(date) - 1L
  )

# -------------------------
# 2) Build country × TH × quarter totals (sum within country across any subcells)
# -------------------------
cth_q <- dat |>
  group_by(country, th, qidx) |>
  summarise(value = sum(value, na.rm = TRUE),
            treat = as.integer(any(treat == 1L)), .groups = "drop")

# Ever-treated by country (independent of TH)
ever_tab <- cth_q |>
  group_by(country) |>
  summarise(ever = as.integer(any(treat == 1L)), .groups = "drop")

cth_q <- cth_q |>
  left_join(ever_tab, by = "country") |>
  mutate(k = qidx - shock_qidx)

# -------------------------
# 3) Helper: aggregate-first event study within a TH group
# -------------------------
plot_one_TH <- function(df, TH_value, B = 500) {
  
  df_th <- df %>% filter(th == TH_value)
  
  # Aggregate first: sum across countries within (ever, qidx)
  gq_sum <- df_th %>%
    group_by(ever, qidx) %>%
    summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
    mutate(k = qidx - shock_qidx) %>%
    filter(k >= -W, k <= W)
  
  # Ensure both groups have baseline at k = -1
  base <- gq_sum %>% filter(k == -1) %>% select(ever, base_total = total)
  if (nrow(base) < 2 || any(base$base_total <= 0)) {
    stop(sprintf("TH==%d: missing/zero baseline at k=-1 for a group.", TH_value))
  }
  
  # Relative-to-baseline (log points and percent), by group
  gq_rel <- gq_sum %>%
    left_join(base, by = "ever") %>%
    mutate(
      ln_rel = log(total) - log(base_total),
      pct    = 100*(exp(ln_rel) - 1),
      group  = ifelse(ever == 1, "Treated", "Control")
    )
  
  # ---------- Plot A: two-line plot ----------
  pA <- ggplot(gq_rel, aes(k, pct, color = group)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_line(size = 1) + geom_point(size = 1.6) +
    scale_color_manual(values = c("Treated"="#1f9ac9","Control"="#fb8161")) +
    labs(
      title = sprintf("TH == %d • Aggregate totals: %% change vs group’s 2010Q4", TH_value),
      subtitle = sprintf("Aggregate-first normalization; window ±%dq; t(0)=2011Q1", W),
      x = "Quarters relative to shock", y = "% change from 2010Q4", color = ""
    ) +
    theme_minimal(base_size = 12)
  
  # ---------- Difference series ----------
  compute_diff_series <- function(df_th_local) {
    gq <- df_th_local %>%
      group_by(ever, qidx) %>%
      summarise(total = sum(value, na.rm = TRUE), .groups = "drop") %>%
      mutate(k = qidx - shock_qidx) %>%
      filter(k >= -W, k <= W)
    base2 <- gq %>% filter(k == -1) %>% select(ever, base_total = total)
    if (nrow(base2) < 2 || any(base2$base_total <= 0)) return(NULL)
    rel <- gq %>%
      left_join(base2, by = "ever") %>%
      mutate(ln_rel = log(total) - log(base_total)) %>%
      select(ever, k, ln_rel) %>%
      pivot_wider(names_from = ever, values_from = ln_rel, names_prefix = "ever_")
    if (!all(c("ever_0","ever_1") %in% names(rel))) return(NULL)
    rel <- rel[order(rel$k), ]
    transform(rel, diff = ever_1 - ever_0)[, c("k","diff")]
  }
  
  # Point estimate
  orig <- compute_diff_series(df_th)
  if (is.null(orig)) stop(sprintf("TH==%d: could not compute difference series.", TH_value))
  
  # Bootstrap: resample countries that have k=-1 in this TH group
  countries_baseline <- df_th %>%
    group_by(country) %>%
    summarise(has_base = any(k == -1), .groups="drop") %>%
    filter(has_base) %>% pull(country)
  
  df_th_base <- df_th %>% filter(country %in% countries_baseline)
  
  boot_mat <- matrix(NA_real_, nrow = nrow(orig), ncol = B,
                     dimnames = list(orig$k, NULL))
  
  for (b in 1:B) {
    samp <- sample(countries_baseline, length(countries_baseline), replace = TRUE)
    d_b  <- compute_diff_series(df_th_base %>% filter(country %in% samp))
    if (!is.null(d_b)) {
      idx <- match(orig$k, d_b$k)
      boot_mat[, b] <- d_b$diff[idx]
    }
  }
  boot_mat <- boot_mat[, colSums(is.na(boot_mat)) == 0, drop = FALSE]
  if (ncol(boot_mat) == 0) stop(sprintf("TH==%d: all bootstrap resamples failed (NA columns).", TH_value))
  
  # ---------- SYMMETRIC CIs via bootstrap SEs + delta method ----------
  z <- qnorm(0.975)
  b_log <- orig$diff                                 # log-point effect
  se_log <- apply(boot_mat, 1, sd, na.rm = TRUE)     # bootstrap SD in log space
  
  estPct <- 100 * (exp(b_log) - 1)                   # percent effect
  sePct  <- 100 * exp(b_log) * se_log                # delta-method SE in percent
  
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
  
  # ---------- Plot B: difference with symmetric 95% CIs ----------
  pB <- ggplot(ci_tab, aes(k, estPct)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = loPct, ymax = hiPct), width = 0) +
    labs(
      title = sprintf("TH == %d • Treated − Control vs 2010Q4 (aggregate-first)", TH_value),
      subtitle = sprintf("Symmetric 95%% CIs: bootstrap SE + delta method (B=%d)", ncol(boot_mat)),
      x = "Quarters from 2011Q1 (k)", y = "Difference vs baseline (percent)"
    ) +
    theme_minimal(base_size = 12)
  
  list(two_line_plot = pA, diff_plot = pB,
       sigpos_k = ci_tab$k[ci_tab$sigpos], ci_table = ci_tab)
}

# -------------------------
# 4) Run for TH==1 and TH==0
# -------------------------
res_TH1 <- plot_one_TH(cth_q, TH_value = 1, B = 500)
res_TH0 <- plot_one_TH(cth_q, TH_value = 0, B = 500)

# Show plots:
print(res_TH1$two_line_plot); print(res_TH1$diff_plot)
print(res_TH0$two_line_plot); print(res_TH0$diff_plot)

cat("\nTH==1: post-treatment quarters with significantly positive difference (symmetric 95%):\n")
print(res_TH1$sigpos_k)
cat("\nTH==0: post-treatment quarters with significantly positive difference (symmetric 95%):\n")
print(res_TH0$sigpos_k)













