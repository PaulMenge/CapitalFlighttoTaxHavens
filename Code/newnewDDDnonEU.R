# ===== Non-FE DiDs & DDD with HC1 CIs | Long.xlsx (Sheet1), COLLAPSED TO 4 CELLS/QTR
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(sandwich)
  library(zoo)
  library(tidyr)
})

# ---------- helper: safely parse European date strings ----------
parse_date_safely <- function(x) {
  if (inherits(x, "Date")) return(x)
  d1 <- suppressWarnings(as.Date(x))                 # ISO first
  d2 <- suppressWarnings(as.Date(x, "%d.%m.%Y"))     # dd.mm.yyyy
  if (sum(is.na(d1)) <= sum(is.na(d2))) d1 else d2
}

# ---------- file path ----------
EXCEL_PATH  <- "C:/Users/paulm/Desktop/Public Projekt/Long.xlsx"
EXCEL_SHEET <- "Sheet1"

# ---- 1) Load & normalize columns ----
raw <- read_excel(EXCEL_PATH, sheet = EXCEL_SHEET) %>% rename_with(tolower)
stopifnot(all(c("re_country","country","date","value","th","treat") %in% names(raw)))

dat0 <- raw %>%
  mutate(
    date        = parse_date_safely(date),
    outcome     = as.numeric(value),
    log_outcome = log(outcome),
    treated     = as.integer(treat),
    th          = as.integer(th),
    yq          = factor(as.character(as.yearqtr(date)))
  ) %>%
  filter(is.finite(log_outcome))

# ---- 1a) Mark EU membership (2011 definition) & DROP EU from the CONTROL (treated == 0) ----
eu_2011 <- c(
  "Bulgaria","Czechia","Denmark","Estonia","Finland","France","Germany","Greece",
  "Hungary","Italy","Latvia","Lithuania","Poland","Portugal","Romania",
  "Slovakia","Slovenia","Spain","Sweden","United Kingdom"
)

# add EU flag
dat0 <- dat0 %>% mutate(eu_2011 = as.integer(country %in% eu_2011))

# quick sanity (before filtering controls)
eu_present     <- sort(unique(dat0$country[dat0$eu_2011 == 1]))
non_eu_present <- sort(unique(dat0$country[dat0$eu_2011 == 0]))
cat("EU countries in data (2011 def.):", length(eu_present), "\n")
cat("Non-EU countries in data:        ", length(non_eu_present), "\n")

# *** CORE CHANGE: exclude EU countries from the *control* group ***
dat0 <- dat0 %>% filter(!(treated == 0 & eu_2011 == 1))

# ---- 2) Event window & POST flag (as before) ----
t0 <- as.Date("2011-03-31")  # 2011Q1 reference
dat0 <- dat0 %>%
  mutate(
    event_q = as.integer((as.yearqtr(date) - as.yearqtr(t0)) * 4L),
    time    = as.integer(event_q >= 4L)  # POST = 2012Q1+
  )

# ---- 2a) DROP 2011 and KEEP analysis window; THEN COLLAPSE (log of sums) ----
dat <- dat0 %>%
  filter(event_q >= -4L, event_q <= 8L, event_q <= -1L | event_q >= 4L) %>%
  group_by(yq, event_q, time, treated, th) %>%
  summarize(outcome_sum = sum(outcome, na.rm = TRUE), .groups = "drop") %>%
  mutate(log_outcome = log(outcome_sum))

cat("\n[Diag] Pre/Post means by block (collapsed, log of sums):\n")
prepost_tbl <- dat %>%
  mutate(period = ifelse(event_q >= 4L, "Post", "Pre")) %>%
  group_by(th, period) %>%
  summarize(mu = mean(log_outcome), .groups = "drop_last") %>%
  tidyr::pivot_wider(names_from = period, values_from = mu) %>%
  mutate(change_pct = 100*(exp(Post - Pre) - 1))
print(prepost_tbl)

# quick sanity: should be 4 rows per quarter
cat("\n[Check] Rows per quarter after collapse (should be 4):\n")
print(dat %>% count(yq) %>% arrange(yq))

# For FE DDD: strict 2010 vs 2012 quarters (collapsed)
dat_win <- dat %>%
  mutate(post2012 = as.integer(event_q %in% 4:7)) %>%
  filter(event_q %in% c(-4,-3,-2,-1, 4,5,6,7))

stopifnot(any(dat_win$post2012==0), any(dat_win$post2012==1))

# ---- 3) Non-FE DDD regression, HC1 CIs (reference) ----
ddd_mod <- lm(log_outcome ~ treated * time * th, data = dat)
V_full  <- sandwich::vcovHC(ddd_mod, type = "HC1")
b_full  <- coef(ddd_mod)
crit    <- qt(0.975, df = df.residual(ddd_mod))

grid <- expand.grid(treated = c(0,1), time = c(0,1), th = c(0,1), KEEP.OUT.ATTRS = FALSE)
TT0   <- delete.response(terms(ddd_mod))
X_all <- model.matrix(TT0, data = grid)
cn <- names(b_full)
X  <- X_all[, cn, drop = FALSE]
b  <- b_full[cn]
V  <- V_full[cn, cn, drop = FALSE]

cell <- function(tr, ti, h) which(grid$treated==tr & grid$time==ti & grid$th==h)
iC0_TH <- cell(0,0,1); iC1_TH <- cell(0,1,1); iT0_TH <- cell(1,0,1); iT1_TH <- cell(1,1,1)
iC0_NT <- cell(0,0,0); iC1_NT <- cell(0,1,0); iT0_NT <- cell(1,0,0); iT1_NT <- cell(1,1,0)

L_DiD_TH <- matrix((X[iT1_TH,]-X[iT0_TH,]) - (X[iC1_TH,]-X[iC0_TH,]), 1)
L_DiD_NT <- matrix((X[iT1_NT,]-X[iT0_NT,]) - (X[iC1_NT,]-X[iC0_NT,]), 1)
L_DDD    <- matrix(L_DiD_TH - L_DiD_NT, 1)

est <- function(L) as.numeric(L %*% b)
se  <- function(L) sqrt(as.numeric(L %*% V %*% t(L)))
to_pct <- function(e, s, crit) {
  lo <- e - crit*s; hi <- e + crit*s
  c(pct = 100*(exp(e)-1), lo = 100*(exp(lo)-1), hi = 100*(exp(hi)-1))
}

th_res <- to_pct(est(L_DiD_TH), se(L_DiD_TH), crit)
nt_res <- to_pct(est(L_DiD_NT), se(L_DiD_NT), crit)
ddd_res<- to_pct(est(L_DDD),    se(L_DDD),    crit)

cat(sprintf("\nNon-FE DiD (TH):     %+.2f%%  [%+.2f, %+.2f]\n", th_res["pct"], th_res["lo"], th_res["hi"]))
cat(sprintf("Non-FE DiD (NonTH):  %+.2f%%  [%+.2f, %+.2f]\n", nt_res["pct"], nt_res["lo"], nt_res["hi"]))
cat(sprintf("Non-FE DDD:          %+.2f%%  [%+.2f, %+.2f]\n\n", ddd_res["pct"], ddd_res["lo"], ddd_res["hi"]))

# =====================================================================
# DESIGN A — Minimal FE DDD: quarter FE + group FE (group = treated×th)
# (on the same collapsed cells)
# =====================================================================
datA <- dat_win %>% mutate(group = interaction(treated, th, drop = TRUE))  # 4 groups: 0.0, 1.0, 0.1, 1.1

fA <- as.formula(
  "log_outcome ~ post2012:treated + post2012:th + post2012:treated:th + factor(yq) + factor(group)"
)

modA <- lm(fA, data = datA)
VA   <- sandwich::vcovHC(modA, type = "HC1")
bA   <- coef(modA); cnA <- names(bA); critA <- qt(0.975, df = df.residual(modA))

# Linear combos (name-proof)
ixA <- function(nm) match(nm, cnA)
L_nonTH <- rep(0, length(bA)); L_nonTH[ixA("post2012:treated")] <- 1
L_TH    <- rep(0, length(bA)); L_TH[ixA("post2012:treated")] <- 1; L_TH[ixA("post2012:treated:th")] <- 1
L_DDD   <- rep(0, length(bA)); L_DDD[ixA("post2012:treated:th")] <- 1

estA <- function(L) as.numeric(crossprod(L, bA))
seA  <- function(L) sqrt(as.numeric(L %*% VA %*% L))
to_pctA <- function(e, s) {
  lo <- e - critA*s; hi <- e + critA*s
  c(pct=100*(exp(e)-1), lo=100*(exp(lo)-1), hi=100*(exp(hi)-1))
}

res_nonTH_A <- to_pctA(estA(L_nonTH), seA(L_nonTH))
res_TH_A    <- to_pctA(estA(L_TH),    seA(L_TH))
res_DDD_A   <- to_pctA(estA(L_DDD),   seA(L_DDD))

cat("DESIGN A (Minimal FE: quarter FE + group FE = treated×th)\n")
cat(sprintf("DiD (NonTH):  %+.2f%%  [%+.2f, %+.2f]\n", res_nonTH_A["pct"], res_nonTH_A["lo"], res_nonTH_A["hi"]))
cat(sprintf("DiD (TH):     %+.2f%%  [%+.2f, %+.2f]\n", res_TH_A["pct"],    res_TH_A["lo"],    res_TH_A["hi"]))
cat(sprintf("DDD:          %+.2f%%  [%+.2f, %+.2f]\n\n", res_DDD_A["pct"],  res_DDD_A["lo"],  res_DDD_A["hi"]))

# ---------- helpers: robust (HC1) coefficient table printed like summary ----------
robust_coef_table <- function(model, Vhc1, digits = 4) {
  b   <- coef(model)
  cn  <- names(b)
  if (is.null(colnames(Vhc1)) || is.null(rownames(Vhc1))) {
    colnames(Vhc1) <- rownames(Vhc1) <- cn[seq_len(nrow(Vhc1))]
  }
  common <- intersect(cn, colnames(Vhc1))
  b      <- b[common]
  Vhc1   <- Vhc1[common, common, drop = FALSE]
  se     <- sqrt(pmax(diag(Vhc1), 0))
  tval   <- b / se
  pval   <- 2 * (1 - pt(abs(tval), df = df.residual(model)))
  out <- data.frame(
    Estimate = b,
    `Std. Error (HC1)` = se,
    `t value` = tval,
    `Pr(>|t|)` = pval,
    check.names = FALSE
  )
  print(round(out, digits))
}

# =========================
# DESIGN A — FULL TABLES
# =========================
cat("\n================ DESIGN A — OLS summary (with FEs) ================\n")
print(summary(modA))
cat("===================================================================\n")

cat("\n================ DESIGN A — HC1-robust coefficients ================\n")
robust_coef_table(modA, VA)
cat("===================================================================\n")

cat("===================================================================\n")
## ================================================================
## FE interpretation (print once, then suppress in paper tables)
## ================================================================
coefA <- coef(modA)
fe_idx <- grepl("^factor\\(group\\)", names(coefA))
fe_coefs <- coefA[fe_idx]
if (length(fe_coefs)) {
  fe_pct <- 100 * (exp(fe_coefs) - 1)
  cat("\n[Info] Group FE magnitudes (percent vs base group level):\n")
  print(round(fe_pct, 1))
  cat("Note: These are time-invariant level gaps; they do NOT affect within-group post–pre contrasts.\n\n")
}

## ================================================================
## Plot: pre-FE (non-FE) within-block DiDs for TH and NonTH (HC1 95% CI)
## ================================================================
did_by_block_pre <- data.frame(
  block = c("TH", "NonTH"),
  pct   = c(th_res["pct"], nt_res["pct"]),
  lo    = c(th_res["lo"],  nt_res["lo"]),
  hi    = c(th_res["hi"],  nt_res["hi"])
)

p_blocks_pre <- ggplot(did_by_block_pre, aes(x = block, y = pct)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6) +
  labs(
    title = "Within-block DiD (Pre-FE, 95% HC1 CI), excl. EU from control",
    x = NULL, y = "% Difference to t=-1"
  ) +
  theme_classic() +
  theme(legend.position = "none")

print(p_blocks_pre)

## ================================================================
## Window sweep: DDD vs different pre/post choices (HC1; 95% CI) — on collapsed panel
## ================================================================
# Build an all-quarters collapsed panel with seasonality tags (using the EU-filtered dat0)
dat_all <- dat0 %>%
  mutate(
    year = as.integer(format(date, "%Y")),
    qtr  = factor(quarters(date), levels = c("Q1","Q2","Q3","Q4"))
  ) %>%
  group_by(yq, event_q, treated, th, year, qtr) %>%
  summarize(outcome_sum = sum(outcome, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    log_outcome = log(outcome_sum),
    time = as.integer(event_q >= 4L)
  )

est_ddd <- function(pre_len = 8, post_len = 8, add_seasonality = TRUE) {
  k_pre  <- seq.int(-pre_len, -1L)
  k_post <- seq.int(4L, 4L + post_len - 1L)
  keep <- dat_all %>%
    filter(event_q %in% c(k_pre, k_post)) %>%
    filter(event_q <= -1L | event_q >= 4L)
  
  if (n_distinct(keep$time) < 2 ||
      n_distinct(keep$treated) < 2 ||
      n_distinct(keep$th) < 2) {
    return(data.frame(pre_len, post_len,
                      DiD_TH=NA_real_, DiD_TH_lo=NA_real_, DiD_TH_hi=NA_real_,
                      DiD_NT=NA_real_, DiD_NT_lo=NA_real_, DiD_NT_hi=NA_real_,
                      DDD=NA_real_,    DDD_lo=NA_real_,    DDD_hi=NA_real_))
  }
  
  qtr_f  <- droplevels(keep$qtr);  include_qtr  <- add_seasonality && nlevels(qtr_f)  >= 2
  year_f <- droplevels(factor(keep$year)); include_year <- add_seasonality && nlevels(year_f) >= 2
  
  f_str <- "log_outcome ~ treated * time * th"
  if (include_qtr)  f_str <- paste0(f_str, " + qtr_f")
  if (include_year) f_str <- paste0(f_str, " + year_f")
  
  keep$qtr_f  <- qtr_f
  keep$year_f <- year_f
  
  mod  <- lm(as.formula(f_str), data = keep)
  Vhc1 <- sandwich::vcovHC(mod, type = "HC1")
  b    <- coef(mod)
  crit <- qt(0.975, df = df.residual(mod))
  
  TT0  <- delete.response(terms(mod))
  mf   <- model.frame(mod)
  lev_q <- if ("qtr_f" %in% names(mf)) levels(mf$qtr_f) else NULL
  lev_y <- if ("year_f"%in% names(mf)) levels(mf$year_f) else NULL
  
  grid <- expand.grid(treated=c(0,1), time=c(0,1), th=c(0,1), KEEP.OUT.ATTRS=FALSE)
  if (!is.null(lev_q)) grid$qtr_f  <- factor(rep(lev_q[1], nrow(grid)), levels=lev_q)
  if (!is.null(lev_y)) grid$year_f <- factor(rep(lev_y[1], nrow(grid)), levels=lev_y)
  X_all <- model.matrix(TT0, data = grid)
  
  k <- Reduce(intersect, list(names(b), colnames(Vhc1), rownames(Vhc1), colnames(X_all)))
  if (!length(k)) {
    return(data.frame(pre_len, post_len,
                      DiD_TH=NA_real_, DiD_TH_lo=NA_real_, DiD_TH_hi=NA_real_,
                      DiD_NT=NA_real_, DiD_NT_lo=NA_real_, DiD_NT_hi=NA_real_,
                      DDD=NA_real_,    DDD_lo=NA_real_,    DDD_hi=NA_real_))
  }
  b <- b[k]; V <- Vhc1[k,k,drop=FALSE]; X <- X_all[,k,drop=FALSE]
  
  cx <- function(tr,ti,hh) which(grid$treated==tr & grid$time==ti & grid$th==hh)
  iC0_TH <- cx(0,0,1); iC1_TH <- cx(0,1,1); iT0_TH <- cx(1,0,1); iT1_TH <- cx(1,1,1)
  iC0_NT <- cx(0,0,0); iC1_NT <- cx(0,1,0); iT0_NT <- cx(1,0,0); iT1_NT <- cx(1,1,0)
  
  L_DiD_TH <- matrix((X[iT1_TH,]-X[iT0_TH,])-(X[iC1_TH,]-X[iC0_TH,]),1)
  L_DiD_NT <- matrix((X[iT1_NT,]-X[iT0_NT,])-(X[iC1_NT,]-X[iC0_NT,]),1)
  L_DDD    <- matrix(L_DiD_TH - L_DiD_NT, 1)
  
  est <- function(L) as.numeric(L %*% b)
  se  <- function(L) sqrt(as.numeric(L %*% V %*% t(L)))
  to_pct <- function(e,s){ lo <- e - crit*s; hi <- e + crit*s; c(pct=100*(exp(e)-1), lo=100*(exp(lo)-1), hi=100*(exp(hi)-1)) }
  
  th <- to_pct(est(L_DiD_TH), se(L_DiD_TH))
  nt <- to_pct(est(L_DiD_NT), se(L_DiD_NT))
  dd <- to_pct(est(L_DDD),    se(L_DiD_NT))  # note: was se(L_DDD) originally; keep consistent
  dd <- to_pct(est(L_DDD),    se(L_DDD))     # correct line
  data.frame(pre_len, post_len,
             DiD_TH=th["pct"], DiD_TH_lo=th["lo"], DiD_TH_hi=th["hi"],
             DiD_NT=nt["pct"], DiD_NT_lo=nt["lo"], DiD_NT_hi=nt["hi"],
             DDD=dd["pct"],    DDD_lo=dd["lo"],    DDD_hi=dd["hi"])
}

# sweep grid of windows (with 6 x-axis ticks)
grid_windows <- expand.grid(
  pre_len  = c(2,4,6,8),
  post_len = c(2,4,6,8,10,12),
  KEEP.OUT.ATTRS = FALSE
)

robust_tbl <- dplyr::bind_rows(lapply(seq_len(nrow(grid_windows)), function(i) {
  pl <- grid_windows$pre_len[i]
  po <- grid_windows$post_len[i]
  est_ddd(pre_len = pl, post_len = po, add_seasonality = TRUE)
}))

# (optional) keep rows with finite estimates
robust_tbl <- robust_tbl[is.finite(robust_tbl$DDD), , drop = FALSE]

# plot: DDD vs post_len (colored by pre_len) — SAME STYLE, just x-axis breaks
p_window <- ggplot(robust_tbl, aes(x = post_len, y = DDD, group = pre_len, color = factor(pre_len))) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_point() + 
  geom_line() +
  geom_errorbar(aes(ymin = DDD_lo, ymax = DDD_hi), width = 0.15) +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
  labs(title = "DDD window choice (95% HC1 CI), excl. EU controls",
       x = "Post length (quarters)", color = "Pre length",
       y = "DDD (% change)") +
  theme_classic()

print(p_window)





## ====== TWFE DDD • Parallel Trends (pre-2011Q1) diagnostics ===================
{
  # ---------- helpers: rank-safe clustered vcov (no extra pkgs) ----------------
  vcov_cluster_1way <- function(mod, cluster) {
    X_full <- model.matrix(mod)
    b_all  <- coef(mod)                  # aliased coefs are NA
    keep   <- !is.na(b_all)              # non-aliased columns only
    X      <- X_full[, keep, drop = FALSE]
    e      <- residuals(mod)
    
    N <- nrow(X); K <- ncol(X)
    G <- length(unique(cluster))
    if (G < 2L) stop("Need at least 2 clusters for cluster-robust SEs.")
    
    U    <- X * as.numeric(e)                          # N x K
    Sg   <- rowsum(U, group = as.factor(cluster))      # G x K
    XtX  <- crossprod(X)
    Xinv <- solve(XtX)                                  # non-singular after dropping aliased cols
    meat <- t(Sg) %*% Sg
    cfac <- (G/(G-1)) * ((N-1)/(N-K))                  # finite-sample correction
    
    V <- cfac * (Xinv %*% meat %*% Xinv)
    colnames(V) <- rownames(V) <- names(b_all)[keep]
    V
  }
  vcov_cluster_2way <- function(mod, cl1, cl2) {
    V1  <- vcov_cluster_1way(mod, cl1)
    V2  <- vcov_cluster_1way(mod, cl2)
    V12 <- vcov_cluster_1way(mod, interaction(as.factor(cl1), as.factor(cl2), drop = TRUE))
    # by construction, dimnames align to non-aliased parameter names
    V1 + V2 - V12
  }
  # Moore–Penrose pseudoinverse via eigen (for possibly singular covariance)
  pinv <- function(S, tol = 1e-10) {
    S <- (S + t(S)) / 2
    ee <- eigen(S, symmetric = TRUE)
    val <- ee$values; vec <- ee$vectors
    keep <- val > max(val, 0) * tol
    if (!any(keep)) return(matrix(0, nrow(S), ncol(S), dimnames = dimnames(S)))
    vec[, keep, drop = FALSE] %*% diag(1/val[keep], sum(keep)) %*% t(vec[, keep, drop = FALSE])
  }
  
  # ---------- Build TWFE DDD event study on your collapsed panel ----------------
  # FE for quarter and for the 4 groups (treated x th)
  dat$yq_f        <- if (!"yq_f" %in% names(dat)) factor(dat$yq) else dat$yq_f
  dat$grp         <- if (!"grp"  %in% names(dat)) interaction(dat$treated, dat$th, drop = TRUE) else dat$grp
  dat$event_q_f   <- factor(dat$event_q)
  dat$event_q_f   <- stats::relevel(dat$event_q_f, ref = "-1")  # omit k = -1
  
  mod_es <- lm(log_outcome ~ yq_f + grp + treated:th:event_q_f, data = dat)
  
  # ---------- Two-way clustered VCOV (time × group) with safe fallbacks --------
  mf_es <- model.frame(mod_es)
  cl_t  <- mf_es$yq_f
  cl_g  <- mf_es$grp
  if (length(unique(cl_t)) >= 2 && length(unique(cl_g)) >= 2) {
    V_es <- vcov_cluster_2way(mod_es, cl_t, cl_g)
  } else if (length(unique(cl_t)) >= 2) {
    V_es <- vcov_cluster_1way(mod_es, cl_t)
    cat("[Info] ES VCOV: time clustering only.\n")
  } else if (length(unique(cl_g)) >= 2) {
    V_es <- vcov_cluster_1way(mod_es, cl_g)
    cat("[Info] ES VCOV: group clustering only.\n")
  } else {
    V_es <- sandwich::vcovHC(mod_es, type = "HC1")
    cat("[Info] ES VCOV: fallback to HC1.\n")
  }
  
  # ---------- Extract event-time coefficients (treated:th:event_q_f*) ----------
  b_all <- coef(mod_es)
  ok_nm <- names(b_all)[!is.na(b_all)]                    # non-aliased names
  V_ok  <- V_es[ok_nm, ok_nm, drop = FALSE]
  b_ok  <- b_all[ok_nm]
  
  evt_idx <- grepl("^treated:th:event_q_f", names(b_ok))
  beta_evt <- b_ok[evt_idx]
  Sigma_evt <- V_ok[evt_idx, evt_idx, drop = FALSE]
  
  # Parse k from names like "treated:th:event_q_f-4"
  get_k <- function(nm) as.integer(sub(".*event_q_f", "", nm))
  k_vec <- vapply(names(beta_evt), get_k, 0L)
  ord   <- order(k_vec)
  k_vec <- k_vec[ord]; beta_evt <- beta_evt[ord]; Sigma_evt <- Sigma_evt[ord, ord, drop = FALSE]
  
  # Identify PRE event times (< 0). Baseline -1 is omitted, so pre set is typically {-4,-3,-2}.
  pre_idx <- which(k_vec < 0)
  stopifnot(length(pre_idx) >= 1)
  
  beta_pre  <- beta_evt[pre_idx]
  Sigma_pre <- Sigma_evt[pre_idx, pre_idx, drop = FALSE]
  k_pre     <- k_vec[pre_idx]
  
  # ---------- (1) Strict Parallel Trends test: H0: all pre betas = 0 -----------
  # Wald statistic with clustered covariance
  W_strict <- as.numeric(t(beta_pre) %*% pinv(Sigma_pre) %*% beta_pre)
  q_strict <- length(beta_pre)
  p_strict <- 1 - pchisq(W_strict, df = q_strict)
  
  cat("\n=== Parallel Trends (strict) — joint pre-period test ===\n")
  cat(sprintf("Wald χ^2(%d) = %.3f,   p = %.4f\n", q_strict, W_strict, p_strict))
  
  # ---------- (2) “Linear-allowed” Parallel Trends (curvature only) ------------
  # Allow a common linear drift in pre: project out span{1, k} from beta_pre, test remaining curvature = 0
  Z      <- cbind(Intercept = 1, k = as.numeric(k_pre))
  M      <- diag(length(k_pre)) - Z %*% solve(crossprod(Z)) %*% t(Z)  # residual-maker
  b_curv <- as.numeric(M %*% beta_pre)
  S_curv <- M %*% Sigma_pre %*% t(M)
  
  # Rank of M (number of tested curvature components)
  r_curv <- qr(M)$rank
  if (r_curv >= 1) {
    W_curv <- as.numeric(t(b_curv) %*% pinv(S_curv) %*% b_curv)
    p_curv <- 1 - pchisq(W_curv, df = r_curv)
    cat("\n=== Parallel Trends (allow linear drift) — curvature test ===\n")
    cat(sprintf("Wald χ^2(%d) = %.3f,   p = %.4f\n", r_curv, W_curv, p_curv))
  } else {
    cat("\n[Curvature test] Not enough distinct pre periods to identify curvature (need ≥ 3 pre k).\n")
  }
  
  # ---------- Plot: pre-event coefficients with clustered 95% CIs --------------
  se_pre <- sqrt(pmax(diag(Sigma_pre), 0))
  to_pct <- function(x) 100 * (exp(x) - 1)
  pre_df <- data.frame(
    k     = k_pre,
    est   = beta_pre,
    lo    = beta_pre - 1.96 * se_pre,
    hi    = beta_pre + 1.96 * se_pre
  ) |>
    mutate(
      estPct = to_pct(est),
      loPct  = to_pct(lo),
      hiPct  = to_pct(hi)
    )
  
  # Fit & overlay the OLS linear trend over pre betas (visual check for linearity)
  lin_fit <- lm(estPct ~ k, data = pre_df)
  k_grid  <- seq(min(pre_df$k), max(pre_df$k), by = 1)
  lin_df  <- data.frame(k = k_grid, fit = predict(lin_fit, newdata = data.frame(k = k_grid)))
  
  p_pre <- ggplot(pre_df, aes(x = k, y = estPct)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_errorbar(aes(ymin = loPct, ymax = hiPct), width = 0) +
    geom_point(size = 2) +
    geom_line(data = lin_df, aes(k, fit), linewidth = 0.8, linetype = 3) +
    scale_x_continuous(breaks = sort(unique(pre_df$k))) +
    labs(
      title = "TWFE DDD — Pre-2011Q1 event-time coefficients",
      subtitle = "95% two-way clustered CIs; dashed line = fitted pre linear trend",
      x = "Event time k (pre only)", y = "Percent vs baseline (k = −1)"
    ) +
    theme_classic(base_size = 12)
  print(p_pre)
  
  # ---------- Compact console summary ------------------------------------------
  cat("\n[Summary] Pre-period coefficients (percent, clustered 95% CI):\n")
  print(pre_df[order(pre_df$k), c("k","estPct","loPct","hiPct")], row.names = FALSE)
}



## ===== TWFE DDD — Pre-trends tests + Honest PT sensitivity (append below) =====

# --- Helpers: rank-safe clustered vcov (1- & 2-way) and PSD fix ----------------
.make_psd <- function(S, eps = 1e-10){
  Ssym <- 0.5*(S + t(S))
  ev   <- eigen(Ssym, symmetric = TRUE)
  vals <- pmax(ev$values, eps)
  Spsd <- ev$vectors %*% diag(vals, nrow = length(vals)) %*% t(ev$vectors)
  dimnames(Spsd) <- dimnames(S)
  Spsd
}
.vcov_cluster_1way <- function(mod, cluster) {
  Xfull <- model.matrix(mod)
  b_all <- coef(mod)                # aliased coefs are NA
  keep  <- !is.na(b_all)            # non-aliased only
  X     <- Xfull[, keep, drop = FALSE]
  e     <- residuals(mod)
  
  N <- nrow(X); K <- ncol(X)
  G <- length(unique(cluster))
  if (G < 2L) stop("Need ≥2 clusters for cluster-robust SEs.")
  
  U    <- X * as.numeric(e)                       # N x K
  Sg   <- rowsum(U, group = as.factor(cluster))   # G x K
  XtX  <- crossprod(X)
  Xinv <- solve(XtX)
  meat <- t(Sg) %*% Sg
  cfac <- (G/(G-1)) * ((N-1)/(N-K))              # finite-sample correction
  
  V <- cfac * (Xinv %*% meat %*% Xinv)
  colnames(V) <- rownames(V) <- names(b_all)[keep]
  .make_psd(V)
}
.vcov_cluster_2way <- function(mod, cl1, cl2) {
  V1  <- .vcov_cluster_1way(mod, cl1)
  V2  <- .vcov_cluster_1way(mod, cl2)
  V12 <- .vcov_cluster_1way(mod, interaction(as.factor(cl1), as.factor(cl2), drop = TRUE))
  .make_psd(V1 + V2 - V12)
}

# --- Build TWFE DDD event-study on your collapsed panel ------------------------
if (!"yq_f" %in% names(dat))        dat$yq_f        <- factor(dat$yq)
if (!"grp"  %in% names(dat))        dat$grp         <- interaction(dat$treated, dat$th, drop = TRUE)
if (!"event_q_f" %in% names(dat))   dat$event_q_f   <- factor(dat$event_q)
dat$event_q_f <- stats::relevel(dat$event_q_f, ref = "-1")  # baseline k = -1

# TWFE DDD ES
mod_es <- lm(log_outcome ~ yq_f + grp + treated:th:event_q_f, data = dat)

# Two-way clustered VCOV (time × group), safe fallbacks
mf_es <- model.frame(mod_es)
cl_t  <- mf_es$yq_f
cl_g  <- mf_es$grp
V_es  <- try(.vcov_cluster_2way(mod_es, cl_t, cl_g), silent = TRUE)
if (inherits(V_es, "try-error")) {
  V_es <- try(.vcov_cluster_1way(mod_es, cl_t), silent = TRUE)
  if (inherits(V_es, "try-error")) V_es <- sandwich::vcovHC(mod_es, type = "HC1")
}
V_es <- .make_psd(V_es)

# --- Extract event-time coefficients and their covariance ----------------------
b_all <- coef(mod_es)
ok_nm <- names(b_all)[!is.na(b_all)]
b_ok  <- b_all[ok_nm]
V_ok  <- V_es[ok_nm, ok_nm, drop = FALSE]

is_evt   <- grepl("^treated:th:event_q_f", names(b_ok))
beta_evt <- b_ok[is_evt]
Sigma    <- V_ok[is_evt, is_evt, drop = FALSE]

# names like "treated:th:event_q_f-4" -> parse k as integer
get_k <- function(nm) as.integer(sub(".*event_q_f", "", nm))
k_all <- vapply(names(beta_evt), get_k, 0L)

# Order by k and split pre/post (baseline -1 omitted by releveling)
ord        <- order(k_all)
k_all      <- k_all[ord]
beta_evt   <- beta_evt[ord]
Sigma      <- Sigma[ord, ord, drop = FALSE]

pre_idx    <- which(k_all < 0)
post_idx   <- which(k_all >= 0)
Sigma_pre  <- Sigma[pre_idx, pre_idx, drop = FALSE]
beta_pre   <- beta_evt[pre_idx]
k_pre      <- k_all[pre_idx]

# --- (1) Joint Wald test: all pre coefficients = 0 (parallel levels) ----------
if (length(pre_idx) >= 1) {
  L      <- diag(length(pre_idx))
  W_stat <- as.numeric(t(beta_pre) %*% solve(.make_psd(Sigma_pre)) %*% beta_pre)
  df     <- length(pre_idx)
  p_wald <- 1 - pchisq(W_stat, df = df)
  cat(sprintf("\n[Pre-period joint test] Wald χ²(%d) = %.2f,  p = %.3f\n", df, W_stat, p_wald))
}

# --- (2) GLS slope test: linear pre-trend = 0 (parallel trends) ---------------
if (length(pre_idx) >= 2) {
  W  <- solve(.make_psd(Sigma_pre))
  k0 <- k_pre - mean(k_pre)                    # center k for interpretability
  denom   <- as.numeric(t(k0) %*% W %*% k0)
  c_vec   <- drop(W %*% k0) / denom           # linear functional selecting slope
  slope   <- as.numeric(crossprod(c_vec, beta_pre))
  var_sl  <- as.numeric(t(c_vec) %*% Sigma_pre %*% c_vec)
  se_sl   <- sqrt(max(var_sl, 0))
  t_sl    <- slope / se_sl
  p_sl    <- 2*(1 - pnorm(abs(t_sl)))
  cat(sprintf("[Pre-period trend] GLS slope (log-pts/quarter) = %.4f, SE = %.4f,  p = %.3f\n",
              slope, se_sl, p_sl))
}

# --- Plot: pre-period coefficients with 95% two-way clustered CIs -------------
pct_ <- function(x) 100*(exp(x) - 1)
pre_df <- data.frame(
  k    = k_pre,
  est  = pct_(beta_pre),
  lo   = pct_(beta_pre - 1.96*sqrt(pmax(diag(Sigma_pre), 0))),
  hi   = pct_(beta_pre + 1.96*sqrt(pmax(diag(Sigma_pre), 0)))
)

p_pre <- ggplot(pre_df, aes(k, est)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 1.9) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", formula = y ~ x) +
  labs(
    title = "TWFE DDD — Pre-2011Q1 event-time coefficients",
    subtitle = "95% two-way clustered CIs; dashed line = fitted pre linear trend",
    x = "Event time k (pre only)", y = "Percent vs baseline (k = −1)"
  ) +
  theme_classic(base_size = 12)
print(p_pre)
## ===== Within-block DiD (no FEs): pre-trend tests + Honest PT sensitivity =====

# small helper: make a covariance PSD
.make_psd <- function(S, eps = 1e-10){
  Ssym <- 0.5*(S + t(S))
  ev   <- eigen(Ssym, symmetric = TRUE)
  vals <- pmax(ev$values, eps)
  Spsd <- ev$vectors %*% diag(vals, nrow = length(vals)) %*% t(ev$vectors)
  dimnames(Spsd) <- dimnames(S)
  Spsd
}
`%||%` <- function(a,b) if (!is.null(a)) a else b
.pct <- function(x) 100*(exp(x) - 1)

run_block <- function(block_flag = 1L, block_name = "TH") {
  # --- 1) Build block panel and event-study (no FEs) --------------------------
  df <- dat %>% dplyr::filter(th == block_flag)
  stopifnot(dplyr::n_distinct(df$treated) >= 2)
  
  df$eq_f <- factor(df$event_q)
  df$eq_f <- stats::relevel(df$eq_f, ref = "-1")
  
  # classic within-block DiD event-study (controls time with eq_f main effects)
  mod <- lm(log_outcome ~ treated + eq_f + treated:eq_f, data = df)
  Vhc1 <- sandwich::vcovHC(mod, type = "HC1")
  
  # --- 2) Extract dynamic DiD betas: treated:eq_f{k} --------------------------
  b_all <- coef(mod)
  keep  <- grepl("^treated:eq_f", names(b_all))
  beta  <- b_all[keep]
  Sig   <- Vhc1[keep, keep, drop = FALSE]
  
  # parse k from names like "treated:eq_f-3"
  k      <- as.integer(sub("treated:eq_f", "", names(beta)))
  ord    <- order(k)
  k      <- k[ord]; beta <- beta[ord]; Sig <- Sig[ord, ord, drop = FALSE]
  
  pre    <- which(k < 0)              # baseline -1 is omitted; pre = {-4,-3,-2,...}
  post   <- which(k >= 0)             # includes k = 0, 1, ...
  
  if (length(pre) == 0L) {
    message("[", block_name, "] No pre periods found. Skipping tests.")
    return(invisible(NULL))
  }
  
  Sig_pre <- .make_psd(Sig[pre, pre, drop = FALSE])
  b_pre   <- beta[pre]
  k_pre   <- k[pre]
  
  # --- 3) Pre-period diagnostics ----------------------------------------------
  # (a) Joint Wald test: all pre coefficients = 0
  W  <- as.numeric(t(b_pre) %*% solve(Sig_pre) %*% b_pre)
  df <- length(pre)
  pW <- 1 - pchisq(W, df = df)
  
  # (b) GLS slope test (linear pre-trend)
  res <- NA; se <- NA; pT <- NA
  if (length(pre) >= 2) {
    Winv  <- solve(Sig_pre)
    k0    <- k_pre - mean(k_pre)
    denom <- as.numeric(t(k0) %*% Winv %*% k0)
    cvec  <- drop(Winv %*% k0) / denom
    res   <- as.numeric(crossprod(cvec, b_pre))
    se    <- sqrt(as.numeric(t(cvec) %*% Sig_pre %*% cvec))
    pT    <- 2*(1 - pnorm(abs(res/se)))
  }
  
  cat(sprintf("\n[%s] Pre-period joint Wald χ²(%d)=%.2f, p=%.3f", block_name, df, W, pW))
  if (!is.na(res)) cat(sprintf(" | GLS slope=%.4f (log-pts/qtr), SE=%.4f, p=%.3f\n", res, se, pT)) else cat("\n")
  
  # --- 4) Plot pre coefficients (percent) with 95%% HC1 CIs + dashed linear fit
  pre_df <- data.frame(
    k  = k_pre,
    est = .pct(b_pre),
    lo  = .pct(b_pre - 1.96*sqrt(pmax(diag(Sig_pre), 0))),
    hi  = .pct(b_pre + 1.96*sqrt(pmax(diag(Sig_pre), 0)))
  )
  p_pre <- ggplot(pre_df, aes(k, est)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_point(size = 1.9) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0) +
    geom_smooth(method = "lm", se = FALSE, linetype = "dashed", formula = y ~ x, color = "steelblue") +
    labs(
      title = sprintf("Within-block DiD (no FEs) — %s pre-event coefficients", block_name),
      subtitle = "95% HC1 CIs; dashed = fitted pre linear trend",
      x = "Event time k (pre only)", y = "Percent vs baseline (k = −1)"
    ) +
    theme_classic(base_size = 12)
  print(p_pre)
  
  # --- 5) Honest Parallel Trends FLCIs for a post horizon (if pkg works) ------
  has_hd <- requireNamespace("HonestDiD", quietly = TRUE)
  if (has_hd && length(post) >= 1) {
    # reorder to [pre, post] for HonestDiD
    beta_stk <- c(b_pre, beta[post])
    Sig_stk  <- .make_psd(rbind(
      cbind(Sig_pre,        Sig[pre, post, drop = FALSE]),
      cbind(Sig[post, pre, drop = FALSE], Sig[post, post, drop = FALSE])
    ))
    numPre  <- length(pre)
    numPost <- length(post)
    
    # choose horizon: last available post k*
    k_post <- k[post]
    j_star <- which(k_post == max(k_post))
    
    # l-vector that picks beta_{k*} among post entries
    l_vec <- { v <- rep(0, numPre + numPost); v[numPre + j_star] <- 1; v }
    
    # OLS CI (blue) for reference
    b_hat <- beta[post][j_star]
    se_b  <- sqrt(max(Sig[post, post, drop = FALSE][j_star, j_star], 0))
    ols_ci <- data.frame(M = 0, loPct = .pct(b_hat - 1.96*se_b), estPct = .pct(b_hat), hiPct = .pct(b_hat + 1.96*se_b))
    
    hd_exports <- getNamespaceExports("HonestDiD")
    M_grid <- c(0, 0.25, 0.5, 1, 2)
    
    get_hd <- function(){
      if ("honest_did" %in% hd_exports) {
        HonestDiD::honest_did(
          beta_hat = beta_stk,
          Sigma    = Sig_stk,
          numPre   = numPre,
          numPost  = numPost,
          l_vec    = l_vec,
          Mvec     = M_grid,
          alpha    = 0.05,
          method   = "FLCI"
        )
      } else if ("createSensitivityResults" %in% hd_exports) {
        HonestDiD::createSensitivityResults(
          betahat = beta_stk,
          sigma   = Sig_stk,
          numPre  = numPre,
          numPost = numPost,
          l_vec   = l_vec,
          Mvec    = M_grid,
          alpha   = 0.05,
          method  = "FLCI"
        )
      } else if ("createSensitivityBounds" %in% hd_exports) {
        HonestDiD::createSensitivityBounds(
          betahat = beta_stk,
          sigma   = Sig_stk,
          numPre  = numPre,
          numPost = numPost,
          l_vec   = l_vec,
          Mvec    = M_grid,
          alpha   = 0.05,
          method  = "FLCI"
        )
      } else stop("No suitable HonestDiD interface.")
    }
    
    hd <- try(get_hd(), silent = TRUE)
    
    if (!inherits(hd, "try-error")) {
      M_out <- hd$Mvec %||% hd$M %||% hd$Mgrid
      L_out <- hd$L    %||% hd$Lower %||% hd$lower %||% hd$ci_lower
      U_out <- hd$U    %||% hd$Upper %||% hd$upper %||% hd$ci_upper
      
      flci_df <- data.frame(
        M      = as.numeric(M_out),
        loPct  = .pct(as.numeric(L_out)),
        estPct = .pct(rep(b_hat, length(M_out))),
        hiPct  = .pct(as.numeric(U_out))
      )
      
      p_flci <- ggplot() +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(data = flci_df, aes(M, estPct), color = "red", size = 1.7) +
        geom_errorbar(data = flci_df, aes(M, ymin = loPct, ymax = hiPct), color = "red", width = 0.03) +
        geom_point(data = ols_ci, aes(M, estPct), color = "steelblue", size = 2) +
        geom_errorbar(data = ols_ci, aes(M, ymin = loPct, ymax = hiPct), color = "steelblue", width = 0.03) +
        scale_x_continuous(breaks = M_grid) +
        labs(
          title = sprintf("Honest PT (Within-block %s, no FEs) at k* = %d", block_name, max(k_post)),
          subtitle = "Blue = OLS 95% CI; Red = FLCIs under curvature bound M",
          x = "M (non-linearity of PT violations)", y = "Percent"
        ) +
        theme_classic(base_size = 12)
      print(p_flci)
      
      cat(sprintf("[%s] Honest PT target k*=%d | OLS CI: [%+.2f, %+.2f] %%  point=%+.2f %%\n",
                  block_name, max(k_post), ols_ci$loPct, ols_ci$hiPct, ols_ci$estPct))
      print(flci_df[,c("M","loPct","hiPct")])
    } else {
      message("[HonestDiD] failed for block ", block_name, " — showing OLS CI only.")
      p_ols <- ggplot(ols_ci, aes(M, estPct)) +
        geom_hline(yintercept = 0, linetype = 2) +
        geom_point(color = "steelblue", size = 2) +
        geom_errorbar(aes(ymin = loPct, ymax = hiPct), color = "steelblue", width = 0.03) +
        scale_x_continuous(breaks = 0, labels = sprintf("k=%d", max(k_post))) +
        labs(title = sprintf("OLS CI (Within-block %s) at k*=%d", block_name, max(k_post)),
             x = "", y = "Percent") +
        theme_classic(base_size = 12)
      print(p_ols)
    }
  } else {
    message("[", block_name, "] HonestDiD not available or no post periods — skipped sensitivity.")
  }
}

# ---- Run for TH and NonTH blocks ---------------------------------------------
run_block(1L, "TH")
run_block(0L, "NonTH")
