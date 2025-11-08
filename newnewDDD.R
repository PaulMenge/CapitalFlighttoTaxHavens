  # ===== Non-FE DiDs & DDD with HC1 CIs | Long.xlsx (Sheet1), COLLAPSED TO 4 CELLS/QTR =====
  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(ggplot2)
    library(sandwich)
    library(zoo)
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
  
  # For FE DDD: strict 2010 vs 2012 quarters (8 quarters total) — still collapsed
  dat_win <- dat %>%
    mutate(post2012 = as.integer(event_q %in% 4:7)) %>%
    filter(event_q %in% c(-4,-3,-2,-1, 4,5,6,7))
  
  #stopifnot(nrow(dat_win) == 32)  # 4 cells × 8 quarters
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
  # (on the same 32 collapsed cells as VERGLEICHALL)
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
      title = "Within-block DiD (Non-FE, 2010 vs 2012; 95% HC1 CI)",
      subtitle = "Estimated from collapsed lm(log_outcome ~ treated*time*th)",
      x = NULL, y = "Percent change (%)"
    ) +
    theme_classic() +
    theme(legend.position = "none")
  
  print(p_blocks_pre)
  
  ## ================================================================
  ## Window sweep: DDD vs different pre/post choices (HC1; 95% CI) — on collapsed panel
  ## ================================================================
  # Build an all-quarters collapsed panel with seasonality tags
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
    dd <- to_pct(est(L_DDD),    se(L_DDD))
    data.frame(pre_len, post_len,
               DiD_TH=th["pct"], DiD_TH_lo=th["lo"], DiD_TH_hi=th["hi"],
               DiD_NT=nt["pct"], DiD_NT_lo=nt["lo"], DiD_NT_hi=nt["hi"],
               DDD=dd["pct"],    DDD_lo=dd["lo"],    DDD_hi=dd["hi"])
  }
  
  # sweep grid of windows (same choices as you used)
  grid_windows <- expand.grid(
    pre_len  = c(2,4,6,8),
    post_len = c(2,4,6,8,10,12),
    KEEP.OUT.ATTRS = FALSE
  )
  # (optional) keep rows with finite estimates
  robust_tbl <- robust_tbl[is.finite(robust_tbl$DDD), , drop = FALSE]
  
  # plot: DDD vs post_len (colored by pre_len) — SAME STYLE, just x-axis breaks
  p_window <- ggplot(robust_tbl, aes(x = post_len, y = DDD, group = pre_len, color = factor(pre_len))) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_point() + 
    geom_line() +
    geom_errorbar(aes(ymin = DDD_lo, ymax = DDD_hi), width = 0.15) +
    scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12)) +
    labs(title = "DDD window choice",
         x = "Post length (quarters)", color = "Pre length",
         y = "DDD (% change, 95% HC1 CI)") +
    theme_classic()
  
  print(p_window)
