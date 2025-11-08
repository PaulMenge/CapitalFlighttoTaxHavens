suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(sandwich)
  library(zoo)
})

## ---- 1) Load & normalize columns ----
dat <- read_excel("C:/Users/paulm/Desktop/Public Projekt/VERLGEICHALL.xlsx", sheet = "Long")
names(dat) <- tolower(names(dat))

# outcome: value or number
if ("value" %in% names(dat)) {
  dat$outcome <- dat$value
} else if ("number" %in% names(dat)) {
  dat$outcome <- dat$number
} else {
  stop("Need an outcome column named 'value' or 'number'.")
}

# treated: from group / treated / treat
if ("group" %in% names(dat)) {
  dat$treated <- as.integer(dat$group %in% c("Treat","T"))
} else if ("treated" %in% names(dat)) {
  dat$treated <- as.integer(dat$treated)
} else if ("treat" %in% names(dat)) {
  dat$treated <- as.integer(dat$treat)
} else {
  stop("Need 'group' (with 'Treat') or a treated indicator ('treated'/'treat').")
}

# TH vs NonTH: from block / th
if ("block" %in% names(dat)) {
  dat$th <- as.integer(dat$block == "TH")
} else if ("th" %in% names(dat)) {
  dat$th <- as.integer(dat$th)
} else {
  stop("Need 'block' (TH/NonTH) or a 'th' indicator (0/1).")
}
# date to Date
dat$date <- as.Date(dat$date)

## ---- 2) Event window ±20q and time dummy ----
t0 <- as.Date("2011-03-31")
dat <- dat %>%
  mutate(
    log_outcome = log(outcome),
    event_q = as.integer((as.yearqtr(date) - as.yearqtr(t0)) * 4)
  ) %>%
  filter(event_q >= -2, event_q <= 20) %>%
  mutate(time = as.integer(event_q >= 0))

# Safety: drop non-positive outcomes if present
dat <- dat %>% filter(is.finite(log_outcome))

stopifnot(any(dat$time==0), any(dat$time==1))  # need both pre and post

## ---- 3) DDD regression (log scale) ----
ddd_mod <- lm(log_outcome ~ treated * time * th, data = dat)

# Robust vcov (HC1). For clustering by country, use vcovCL(ddd_mod, cluster = dat$country)
V_full <- sandwich::vcovHC(ddd_mod, type = "HC1")
b_full <- coef(ddd_mod)
crit   <- qt(0.975, df = df.residual(ddd_mod))

## ---- 4) Model-matrix contrasts (name-proof) ----
# Build the 8 cells and align design to the model’s coefficient order
grid <- expand.grid(treated = c(0,1), time = c(0,1), th = c(0,1))
TT0   <- delete.response(terms(ddd_mod))
X_all <- model.matrix(TT0, data = grid)

cn <- names(b_full)
X  <- X_all[, cn, drop = FALSE]
b  <- b_full[cn]
V  <- V_full[cn, cn, drop = FALSE]

idx <- function(tr, ti, th1) which(grid$treated==tr & grid$time==ti & grid$th==th1)

# TH cells
iC0_TH <- idx(0,0,1); iC1_TH <- idx(0,1,1); iT0_TH <- idx(1,0,1); iT1_TH <- idx(1,1,1)
# NonTH cells
iC0_NT <- idx(0,0,0); iC1_NT <- idx(0,1,0); iT0_NT <- idx(1,0,0); iT1_NT <- idx(1,1,0)

# Log-scale contrasts: DiD within block = (T post-pre) - (C post-pre)
L_DiD_TH <- matrix( (X[iT1_TH,]-X[iT0_TH,]) - (X[iC1_TH,]-X[iC0_TH,]), 1 )
L_DiD_NT <- matrix( (X[iT1_NT,]-X[iT0_NT,]) - (X[iC1_NT,]-X[iC0_NT,]), 1 )
L_DDD    <- matrix( L_DiD_TH - L_DiD_NT, 1 )

# Estimates & robust SEs
est_DiD_TH <- as.numeric(L_DiD_TH %*% b)
est_DiD_NT <- as.numeric(L_DiD_NT %*% b)
est_DDD    <- as.numeric(L_DDD    %*% b)

se_DiD_TH <- sqrt(as.numeric(L_DiD_TH %*% V %*% t(L_DiD_TH)))
se_DiD_NT <- sqrt(as.numeric(L_DiD_NT %*% V %*% t(L_DiD_NT)))
se_DDD    <- sqrt(as.numeric(L_DDD    %*% V %*% t(L_DDD)))

to_pct <- function(est, se) {
  lo <- est - crit*se; hi <- est + crit*se
  c(pct = (exp(est)-1)*100, lo = (exp(lo)-1)*100, hi = (exp(hi)-1)*100)
}
th_res <- to_pct(est_DiD_TH, se_DiD_TH)
nt_res <- to_pct(est_DiD_NT, se_DiD_NT)
ddd_res<- to_pct(est_DDD,   se_DDD)

## ---- 5) Plots: in-group DiDs & DDD (percent, 95% CI) ----
did_by_block <- data.frame(
  block = c("TH","NonTH"),
  pct   = c(th_res["pct"], nt_res["pct"]),
  lo    = c(th_res["lo"],  nt_res["lo"]),
  hi    = c(th_res["hi"],  nt_res["hi"])
)

ddd_df <- data.frame(
  stat = "DDD = DiD(TH) − DiD(NonTH)",
  pct  = ddd_res["pct"], lo = ddd_res["lo"], hi = ddd_res["hi"]
)

#DiD within blocks
p_blocks <- ggplot(did_by_block, aes(x = block, y = pct)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6) +
  labs(title = "Within-group DiD Treated - Control",
       x = NULL, y = "Percent change (%)") +
  theme_classic() +
  theme(legend.position = "none")
print(p_blocks)

#Triple difference
p_ddd <- ggplot(ddd_df, aes(x = stat, y = pct)) +
  geom_hline(yintercept = 0, color = "grey70") +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, linewidth = 0.6) +
  coord_flip() +
  labs(title = "DDD: (T−C)_TH − (T−C)_NonTH",
       x = NULL, y = "Percent change (%)") +
  theme_classic()
print(p_ddd)

cat(sprintf("DiD (TH):     %+.2f%%  [%+.2f, %+.2f]\n", th_res["pct"], th_res["lo"], th_res["hi"]))
cat(sprintf("DiD (NonTH):  %+.2f%%  [%+.2f, %+.2f]\n", nt_res["pct"], nt_res["lo"], nt_res["hi"]))
cat(sprintf("DDD:          %+.2f%%  [%+.2f, %+.2f]\n", ddd_res["pct"], ddd_res["lo"], ddd_res["hi"]))
summary(ddd_mod)
