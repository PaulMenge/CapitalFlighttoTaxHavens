suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(lubridate); library(zoo)
  library(ggplot2); library(tidyr)
  set.seed(123)
})

## --------------------------------------------
## 0) Load & normalize
## --------------------------------------------
if (!exists("shock_qidx")) shock_qidx <- 2011L*4 + 1L - 1L  # 2011Q1 baseline, not used below

dat <- read_excel("C:/Users/paulm/Desktop/Public Projekt/Long.xlsx", sheet = "Sheet1") |>
  rename_with(tolower) |>
  mutate(
    date  = as.Date(date),
    th    = as.integer(th),
    treat = as.integer(treat),
    qidx  = lubridate::year(date)*4 + lubridate::quarter(date) - 1L
  )

## --------------------------------------------
## 1) Build country × TH × quarter totals
## --------------------------------------------
cth_q <- dat |>
  group_by(country, th, qidx) |>
  summarise(value = sum(value, na.rm = TRUE),
            treat = as.integer(any(treat == 1L)), .groups = "drop")

ever_tab <- cth_q |>
  group_by(country) |>
  summarise(ever = as.integer(any(treat == 1L)), .groups = "drop")

cth_q <- cth_q |>
  left_join(ever_tab, by = "country") |>
  mutate(k = qidx - shock_qidx)

## --------------------------------------------
## 2) Construct 2010 (pre) vs 2012 (post) panel
## --------------------------------------------
df_1012 <- cth_q %>%
  mutate(
    ln_y  = ifelse(value > 0, log(value), NA_real_),
    year  = floor(qidx / 4),
    qtr   = (qidx %% 4) + 1L,
    yq    = paste(year, qtr, sep = "Q"),
    post  = as.integer(year == 2012L),  # post = 2012Q1–Q4; pre = 2010Q1–Q4
    treat = as.integer(treat),
    th    = as.integer(th),
    group = interaction(treat, th, drop = TRUE)
  ) %>%
  filter(year %in% c(2010L, 2012L), is.finite(ln_y))

stopifnot(dplyr::n_distinct(df_1012$yq) >= 2L,
          dplyr::n_distinct(df_1012$group) >= 2L)

## =========================================================
## Non-FE DiDs & DDD with tighter (HC1) 95% CIs
## =========================================================
suppressPackageStartupMessages({ library(sandwich); library(ggplot2); library(dplyr) })

# Build the minimal regression frame (2010 vs 2012 already in df_1012)
dat_reg <- df_1012 %>%
  transmute(
    ln_y = ln_y,
    treated = as.integer(treat),
    time    = as.integer(year == 2012L),  # post
    th      = as.integer(th)
  )

# Non-FE DDD regression
ddd_mod <- lm(ln_y ~ treated * time * th, data = dat_reg)

# HC1 robust vcov (tighter than country-bootstrap on aggregates)
V_full <- sandwich::vcovHC(ddd_mod, type = "HC1")
b_full <- coef(ddd_mod)
crit   <- qt(0.975, df = df.residual(ddd_mod))

# Name-proof contrasts via model matrix on 8 cells
grid <- expand.grid(treated = c(0,1), time = c(0,1), th = c(0,1), KEEP.OUT.ATTRS = FALSE)
TT0   <- delete.response(terms(ddd_mod))
X_all <- model.matrix(TT0, data = grid)

cn <- names(b_full)
X  <- X_all[, cn, drop = FALSE]
b  <- b_full[cn]
V  <- V_full[cn, cn, drop = FALSE]

idx <- function(tr, ti, th1) which(grid$treated==tr & grid$time==ti & grid$th==th1)

# Cells
iC0_TH <- idx(0,0,1); iC1_TH <- idx(0,1,1); iT0_TH <- idx(1,0,1); iT1_TH <- idx(1,1,1)
iC0_NT <- idx(0,0,0); iC1_NT <- idx(0,1,0); iT0_NT <- idx(1,0,0); iT1_NT <- idx(1,1,0)

# Log-scale DiDs
L_DiD_TH <- matrix( (X[iT1_TH,]-X[iT0_TH,]) - (X[iC1_TH,]-X[iC0_TH,]), 1 )
L_DiD_NT <- matrix( (X[iT1_NT,]-X[iT0_NT,]) - (X[iC1_NT,]-X[iC0_NT,]), 1 )
L_DDD    <- matrix( L_DiD_TH - L_DiD_NT, 1 )

est_DiD_TH <- as.numeric(L_DiD_TH %*% b)
est_DiD_NT <- as.numeric(L_DiD_NT %*% b)
est_DDD    <- as.numeric(L_DDD    %*% b)

se_DiD_TH <- sqrt(as.numeric(L_DiD_TH %*% V %*% t(L_DiD_TH)))
se_DiD_NT <- sqrt(as.numeric(L_DiD_NT %*% V %*% t(L_DiD_NT)))
se_DDD    <- sqrt(as.numeric(L_DDD    %*% V %*% t(L_DDD)))

to_pct <- function(est, se, crit) {
  lo <- est - crit*se; hi <- est + crit*se
  c(pct = (exp(est)-1)*100, lo = (exp(lo)-1)*100, hi = (exp(hi)-1)*100)
}
th_res <- to_pct(est_DiD_TH, se_DiD_TH, crit)
nt_res <- to_pct(est_DiD_NT, se_DiD_NT, crit)
ddd_res<- to_pct(est_DDD,    se_DDD,    crit)

# Print results
cat(sprintf("DiD (TH):     %+.2f%%  [%+.2f, %+.2f]\n", th_res["pct"], th_res["lo"], th_res["hi"]))
cat(sprintf("DiD (NonTH):  %+.2f%%  [%+.2f, %+.2f]\n", nt_res["pct"], nt_res["lo"], nt_res["hi"]))
cat(sprintf("DDD:          %+.2f%%  [%+.2f, %+.2f]\n", ddd_res["pct"], ddd_res["lo"], ddd_res["hi"]))

# Plots (compact, 95% CIs)
did_by_block <- data.frame(
  block = c("TH","NonTH"),
  pct   = c(th_res["pct"], nt_res["pct"]),
  lo    = c(th_res["lo"],  nt_res["lo"]),
  hi    = c(th_res["hi"],  nt_res["hi"])
)

p_blocks <- ggplot(did_by_block, aes(x = block, y = pct)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6) +
  labs(title = "Within-block DiD (Post:2012Q1–Q4 vs Pre:2010Q1–Q4; 95% HC1 CI)",
       x = NULL, y = "Percent change (%)") +
  theme_classic() +
  theme(legend.position = "none")
print(p_blocks)

ddd_df <- data.frame(
  stat = "DDD", pct = ddd_res["pct"], lo = ddd_res["lo"], hi = ddd_res["hi"]
)
p_ddd <- ggplot(ddd_df, aes(x = stat, y = pct)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6) +
  coord_flip() +
  labs(title = "DDD: (T−C)_TH − (T−C)_NonTH (95% HC1 CI)",
       x = NULL, y = "Percent change (%)") +
  theme_classic()
print(p_ddd)

# Full regression table for reference
print(summary(ddd_mod))






## ============================================
## 4) FE DDD (your “real” spec) — unchanged
## ============================================
## Helpers for robust SEs (base R only)
pinv <- function(A, tol = sqrt(.Machine$double.eps)) {
  s <- svd(A); d <- s$d; r <- sum(d > tol * max(d))
  if (r == 0) return(matrix(0, nrow(A), ncol(A)))
  s$u[, 1:r, drop=FALSE] %*% (diag(1/d[1:r], r)) %*% t(s$v[, 1:r, drop=FALSE])
}
safe_XtX_inv <- function(X) {
  XtX <- crossprod(X)
  out <- try(chol2inv(chol(XtX)), silent = TRUE)
  if (inherits(out, "try-error")) out <- pinv(XtX)
  out
}
near_psd <- function(V) {
  V <- (V + t(V)) / 2
  eg <- eigen(V, symmetric = TRUE)
  vals <- pmax(eg$values, 0)
  Vpsd <- eg$vectors %*% (vals * t(eg$vectors))
  dimnames(Vpsd) <- dimnames(V)
  Vpsd
}
vcov_cluster_1way <- function(X, u, cluster_id) {
  M  <- safe_XtX_inv(X)
  cl <- as.factor(cluster_id); G <- nlevels(cl); ids <- levels(cl)
  meat <- matrix(0, ncol(X), ncol(X))
  for (g in ids) {
    idx <- which(cl == g); Xg <- X[idx, , drop=FALSE]; ug <- u[idx]
    Sg <- crossprod(Xg, ug); meat <- meat + Sg %*% t(Sg)
  }
  N <- nrow(X); K <- ncol(X)
  dfc <- (N - 1) / pmax(N - K, 1)
  c_g <- if (G > 1) G / (G - 1) else 1
  V <- M %*% (dfc * c_g * meat) %*% M
  dimnames(V) <- list(colnames(X), colnames(X))
  V
}
vcov_cluster_2way_or_time_fallback <- function(X, u, cl_group, cl_time, target_name) {
  Vg  <- vcov_cluster_1way(X, u, cl_group)
  Vt  <- vcov_cluster_1way(X, u, cl_time)
  Vgt <- vcov_cluster_1way(X, u, interaction(cl_group, cl_time, drop = TRUE))
  V2  <- near_psd(Vg + Vt - Vgt)
  se2 <- sqrt(pmax(diag(V2), 0))
  se_target <- se2[match(target_name, colnames(V2))]
  if (!is.finite(se_target) || se_target < 1e-10) list(V = near_psd(Vt), used = "time-only")
  else                                            list(V = V2,          used = "two-way (group & time)")
}

# FE model
f_core <- ln_y ~ post:treat + post:th + post:treat:th + factor(group) + factor(yq)
mod    <- lm(f_core, data = df_1012)

cat("\n================ OLS with FE: Full output ================\n")
print(summary(mod))
cat("==========================================================\n")

# Cluster-robust DDD (two-way with fallback)
mf <- model.frame(mod); TT <- delete.response(terms(mod))
X  <- model.matrix(TT, mf); b <- coef(mod); X <- X[, names(b), drop = FALSE]; u <- resid(mod)

beta_name <- "post:treat:th"
Vinfo <- vcov_cluster_2way_or_time_fallback(
  X, u, cl_group = df_1012$group, cl_time = df_1012$yq, target_name = beta_name
)
V  <- Vinfo$V
se <- sqrt(pmax(diag(V), 0))

b123  <- b[beta_name]
se123 <- se[match(beta_name, names(b))]
z     <- qnorm(0.975)
pct   <- 100 * (exp(b123) - 1)
loPct <- 100 * (exp(b123 - z*se123) - 1)
hiPct <- 100 * (exp(b123 + z*se123) - 1)

cat("\n================ DDD (cluster-robust) =====================\n")
cat(sprintf("VCOV used: %s\n", Vinfo$used))
cat(sprintf("DDD (post×treat×th): b = %.6f, se = %.6f  →  %.2f%%  [%.2f, %.2f] (95%%)\n",
            b123, se123, pct, loPct, hiPct))
cat("==========================================================\n")

# Compact robust table (post terms)
terms_of_interest <- c("post:treat", "post:th", "post:treat:th")
keep <- intersect(terms_of_interest, names(b))
rob_tab <- data.frame(term = keep, estimate = b[keep], se = se[match(keep, names(b))])
rob_tab$t     <- with(rob_tab, estimate / se)
rob_tab$p_N   <- 2 * (1 - pnorm(abs(rob_tab$t)))
rob_tab$pct   <- 100 * (exp(rob_tab$estimate) - 1)
rob_tab$lo95  <- 100 * (exp(rob_tab$estimate - qnorm(0.975)*rob_tab$se) - 1)
rob_tab$hi95  <- 100 * (exp(rob_tab$estimate + qnorm(0.975)*rob_tab$se) - 1)
print(rob_tab, row.names = FALSE, digits = 4)
