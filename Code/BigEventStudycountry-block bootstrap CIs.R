suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(lubridate); library(zoo)
  library(ggplot2); library(tidyr)
})

# Params
W <- 25
shock_qidx <- 2011L*4 + 1L - 1L     # 2011Q1
B <- 1999                            # bootstrap reps
set.seed(123)

# Load & normalize
dat <- read_excel("C:/Users/paulm/Desktop/Public Projekt/Long.xlsx", sheet = "Sheet1") |>
  rename_with(tolower) |>
  mutate(
    date  = as.Date(date),
    qidx  = year(date)*4 + quarter(date) - 1L
  )

# country-quarter totals (ignore TH inside country)
country_q <- dat |>
  group_by(country, qidx) |>
  summarise(value = sum(value, na.rm = TRUE),
            treat = as.integer(any(treat == 1L)), .groups = "drop")

# ever-treated by country
ever_tab <- country_q |>
  group_by(country) |>
  summarise(ever = as.integer(any(treat == 1L)), .groups = "drop")

country_q <- country_q |>
  left_join(ever_tab, by = "country") |>
  mutate(k = qidx - shock_qidx)

# Aggregate first: sum across countries within (ever, quarter)
gq_sum <- country_q |>
  group_by(ever, qidx) |>
  summarise(total = sum(value, na.rm = TRUE), .groups = "drop") |>
  mutate(k = qidx - shock_qidx) |>
  filter(k >= -W, k <= W)

# Baselines at k = -1 for both groups
base <- gq_sum |> filter(k == -1) |> select(ever, base_total = total)
stopifnot(nrow(base) == 2, all(base$base_total > 0))

# Two-line plot (percent from baseline for each group)
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
    title = "Aggregate totals: % change vs baseline (t=-1)",
    x = "Quarters from 2011Q1", y = "% change from 2010Q4", color = ""
  ) +
  theme_minimal(base_size = 12)
print(p_lines)

# Difference series Treated−Control vs baseline
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

orig <- compute_diff_series(country_q)
stopifnot(!is.null(orig))

# Country-block bootstrap on the micro data
countries_T <- unique(country_q$country[country_q$ever == 1])
countries_C <- unique(country_q$country[country_q$ever == 0])
nT <- length(countries_T); nC <- length(countries_C)

boot_mat <- matrix(NA_real_, nrow = nrow(orig), ncol = B,
                   dimnames = list(orig$k, NULL))

for (b in 1:B) {
  samp_T <- sample(countries_T, nT, replace = TRUE)
  samp_C <- sample(countries_C, nC, replace = TRUE)
  df_b <- dplyr::bind_rows(
    country_q %>% filter(country %in% samp_T),
    country_q %>% filter(country %in% samp_C)
  )
  d_b <- compute_diff_series(df_b)
  if (!is.null(d_b)) {
    idx <- match(orig$k, d_b$k)
    boot_mat[, b] <- d_b$diff[idx]
  }
}
boot_mat <- boot_mat[, colSums(is.na(boot_mat)) == 0, drop = FALSE]
if (ncol(boot_mat) == 0) stop("All bootstrap resamples failed (NA columns). Check data/baseline coverage.")

# --- SYMMETRIC CIs via bootstrap SEs + delta method (percent scale) ---
z <- qnorm(0.975)

# log-point effect and its bootstrap SD
b_log <- orig$diff
se_log <- apply(boot_mat, 1, sd, na.rm = TRUE)

# percent effect and symmetric SE (delta method)
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

ggplot(ci_tab, aes(k, estPct)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_point(size = 1.8) +
  geom_errorbar(aes(ymin = loPct, ymax = hiPct), width = 0) +
  labs(
    title = "Aggregate-first Event Study",
    x = "Quarters from 2011/Q1", y = "% Difference vs t=-1"
  ) +
  theme_minimal(base_size = 12)

cat("Post-treatment k with significantly positive Treated−Control difference (symmetric 95%):\n")
print(ci_tab$k[ci_tab$sigpos])
