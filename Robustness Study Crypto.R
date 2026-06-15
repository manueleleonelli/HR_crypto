# =============================================================================
# crypto_extreme_robustness_v2.R
# Companion pipeline to crypto_extreme_ggm_pipeline_v18.R
#
# v2 changes (FIGURES ONLY; all statistics identical to the original):
# - Shared publication theme theme_pub() with larger fonts, matching v18.
# - Thicker lines, larger points, bigger and clearer legends, dpi 200.
# - SHOW_TITLES toggle: in-figure titles/subtitles OFF by default, since the
#   LaTeX captions already carry them. Set TRUE for standalone viewing.
# - chi(q) figure: solid/dashed (empirical vs Gaussian) is now a proper
#   linetype legend rather than text in a subtitle, so it survives the toggle.
#
# Produces the robustness and economic-significance analyses agreed for the
# IRFA submission. Reuses the SAME data download and GJR-GARCH filtering as
# v18, so it runs standalone, and caches the standardised residuals so it can
# be re-run cheaply.
#
# MODULES (toggle with the RUN_* flags below):
#   A. Threshold robustness   p in {0.10, 0.15, 0.20}          -> Appendix B
#   B. Window robustness      window in {500, 750, 1000} days  -> Appendix B
#   C. Economic significance  joint-exceedance, chi(q), portfolio tail
#                                                              -> Results 4.x
#   D. Asymmetry test         within-window block bootstrap CI on chi^- - chi^+
#                             + periphery edge-stability       -> Results / App B
#   E. Gaussian benchmark     frequentist graphical lasso vs BGGM density
#                                                              -> Appendix B
#
# All numeric outputs go to output/crypto_robustness_results.xlsx;
# figures go to figures/.
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
required_pkgs <- c(
  "quantmod", "rugarch", "graphicalExtremes", "igraph",
  "openxlsx", "ggplot2", "dplyr", "tidyr", "purrr",
  "xts", "MASS", "huge"
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
set.seed(42)

# ---- Master switches --------------------------------------------------------
RUN_THRESHOLD   <- TRUE
RUN_WINDOW      <- TRUE
RUN_ECON        <- TRUE
RUN_ASYMMETRY   <- TRUE
RUN_GAUSS_BENCH <- TRUE
USE_CACHE       <- TRUE     # reuse cached residuals if present

# ---- 1. Configuration (mirrors v18) -----------------------------------------
COINS <- c("BTC-USD","ETH-USD","XRP-USD","BNB-USD","SOL-USD",
           "ADA-USD","DOGE-USD","TRX-USD","LINK-USD","DOT-USD",
           "AVAX-USD","LTC-USD","XLM-USD")
SHORT_NAMES <- c("BTC","ETH","XRP","BNB","SOL","ADA","DOGE",
                 "TRX","LINK","DOT","AVAX","LTC","XLM")
SECTORS <- c(BTC="Currency", LTC="Currency", DOGE="Currency",
             ETH="SmartContract", SOL="SmartContract", ADA="SmartContract",
             AVAX="SmartContract", XRP="Payment", XLM="Payment", TRX="Payment",
             BNB="Infrastructure", LINK="Infrastructure", DOT="Infrastructure")

START_DATE <- as.Date("2020-10-01")
END_DATE   <- Sys.Date()

STEP_DAYS  <- 15L
RHO_GRID   <- seq(0.01, 0.15, length.out = 15)

THRESHOLD_GRID <- c(0.10, 0.15, 0.20)   # p; 0.05 deliberately excluded (too few joint exceedances per window)
WINDOW_GRID    <- c(500L, 750L, 1000L)
BASE_WINDOW    <- 750L
BASE_P         <- 0.20

B_BOOT     <- 300L          # bootstrap replicates for the asymmetry test
BLOCK_LEN  <- 20L           # moving-block length (trading days)
N_SIM      <- 200000L       # simulation size for the Gaussian-copula comparisons

OUT_DIR <- "output"; FIG_DIR <- "figures"
dir.create(OUT_DIR, showWarnings = FALSE); dir.create(FIG_DIR, showWarnings = FALSE)
CACHE <- file.path(OUT_DIR, "Z_residuals.rds")

# ---- 1b. Figure styling (matches crypto_extreme_ggm_pipeline_v18.R) ---------
# Keep in-figure titles off for the paper (captions cover them); flip to TRUE
# to show titles when viewing the PNGs on their own.
SHOW_TITLES <- FALSE
FIG_DPI     <- 200

t_if <- function(x) if (SHOW_TITLES) x else NULL

theme_pub <- function(base_size = 15) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title      = element_text(face = "bold", size = rel(1.15),
                                     margin = margin(b = 6)),
      plot.subtitle   = element_text(size = rel(0.95), colour = "grey35",
                                     margin = margin(b = 8)),
      axis.title      = element_text(size = rel(1.0)),
      axis.title.x    = element_text(margin = margin(t = 6)),
      axis.title.y    = element_text(margin = margin(r = 6)),
      axis.text       = element_text(size = rel(0.85), colour = "grey25"),
      strip.text      = element_text(face = "bold", size = rel(0.95),
                                     margin = margin(4, 4, 4, 4)),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold", size = rel(0.9)),
      legend.text     = element_text(size = rel(0.85)),
      legend.key.width  = unit(1.6, "lines"),
      legend.key.height = unit(1.0, "lines"),
      panel.grid.minor = element_blank(),
      panel.spacing    = unit(1.0, "lines")
    )
}

save_fig <- function(name, plot, width, height) {
  ggsave(file.path(FIG_DIR, name), plot, width = width, height = height,
         dpi = FIG_DPI, bg = "white")
}

# ---- 2. Data + GJR-GARCH filtering (reused from v18) ------------------------
get_residuals <- function() {
  if (USE_CACHE && file.exists(CACHE)) {
    message("Loading cached residuals from ", CACHE)
    return(readRDS(CACHE))
  }
  message("Downloading prices ...")
  price_list <- list()
  for (i in seq_along(COINS)) {
    nm <- SHORT_NAMES[i]
    dat <- tryCatch(getSymbols(COINS[i], src="yahoo", from=START_DATE,
                               to=END_DATE, auto.assign=FALSE),
                    error=function(e) NULL)
    if (!is.null(dat)) price_list[[nm]] <- Ad(dat)
  }
  prices <- do.call(merge, price_list); colnames(prices) <- names(price_list)
  prices <- na.omit(na.locf(prices, na.rm = FALSE))
  rets   <- na.omit(diff(log(prices))[-1, ])
  
  message("GJR-GARCH(1,1)-t filtering ...")
  spec_gjr <- ugarchspec(variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
                         mean.model=list(armaOrder=c(1,0), include.mean=TRUE),
                         distribution.model="std")
  spec_fb  <- ugarchspec(variance.model=list(model="sGARCH", garchOrder=c(1,1)),
                         mean.model=list(armaOrder=c(1,0), include.mean=TRUE),
                         distribution.model="std")
  n <- nrow(rets); d <- ncol(rets); nms <- colnames(rets)
  Zmat <- matrix(NA, n, d, dimnames=list(NULL, nms))
  for (j in seq_len(d)) {
    fit <- tryCatch(ugarchfit(spec_gjr, as.numeric(rets[,j]), solver="hybrid"),
                    error=function(e) ugarchfit(spec_fb, as.numeric(rets[,j]), solver="hybrid"))
    Zmat[,j] <- as.numeric(residuals(fit, standardize=TRUE))
  }
  out <- list(Z = Zmat, dates = index(rets))
  saveRDS(out, CACHE)
  out
}

filtered <- get_residuals()
Z <- filtered$Z[, SHORT_NAMES]      # enforce column order
dates <- filtered$dates
d <- ncol(Z); n_total <- nrow(Z)
message(sprintf("Residuals: %d obs x %d assets, %s to %s",
                n_total, d, format(min(dates)), format(max(dates))))

# ---- 3. Shared helpers ------------------------------------------------------

# Husler-Reiss pairwise chi from a variogram matrix
gamma_to_chi <- function(G) { chi <- 2 - 2*pnorm(sqrt(G)/2); diag(chi) <- 1; chi }

# 10-fold CV selection of the global penalty at a given threshold p
select_rho <- function(Zmat, p, grid = RHO_GRID, K = 10) {
  pick <- function(Y) {
    folds <- sample(rep(1:K, length.out = nrow(Y)))
    ll <- matrix(NA, length(grid), K)
    for (k in 1:K) {
      tr <- Y[folds != k, , drop=FALSE]; te <- Y[folds == k, , drop=FALSE]
      fit <- tryCatch(eglearn(tr, rholist=grid, reg_method="glasso", complete_Gamma=TRUE),
                      error=function(e) NULL)
      if (is.null(fit)) next
      for (j in seq_along(grid)) if (!is.null(fit$Gamma[[j]]))
        ll[j,k] <- tryCatch(loglik_HR(data=te, Gamma=fit$Gamma[[j]], graph=fit$graph[[j]])[1],
                            error=function(e) NA)
    }
    avg <- suppressWarnings(rowMeans(ll, na.rm=TRUE))
    ok  <- which(is.finite(avg))
    if (!length(ok)) return(grid[round(length(grid)/2)])
    grid[ok[which.max(avg[ok])]]
  }
  list(rho_pos = pick(data2mpareto(Zmat, p)),
       rho_neg = pick(data2mpareto(-Zmat, p)))
}

# Per-window summary metrics at a given (window, p, rho); optionally add a
# frequentist graphical-lasso Gaussian density for the benchmark module.
graph_density   <- function(g) igraph::edge_density(g)
graph_modularity<- function(g) { if (igraph::ecount(g) == 0) return(0)
  igraph::modularity(igraph::cluster_louvain(g)) }

glasso_gauss_density <- function(Zw) {
  h  <- tryCatch(huge(as.matrix(Zw), method="glasso", nlambda=30, verbose=FALSE),
                 error=function(e) NULL)
  if (is.null(h)) return(NA_real_)
  hs <- tryCatch(huge.select(h, criterion="ebic", verbose=FALSE),
                 error=function(e) NULL)
  if (is.null(hs)) return(NA_real_)
  A <- as.matrix(hs$refit); diag(A) <- 0
  (sum(A != 0)/2) / choose(ncol(Zw), 2)
}

run_dynamic_summary <- function(Zmat, dts, window_days, p, rho_pos, rho_neg,
                                with_glasso = FALSE) {
  starts <- seq(1, nrow(Zmat) - window_days + 1, by = STEP_DAYS)
  rows <- vector("list", length(starts))
  for (w in seq_along(starts)) {
    idx <- starts[w]:(starts[w] + window_days - 1)
    Zw  <- Zmat[idx, ]; dw <- dts[idx]; mid <- dw[round(length(dw)/2)]
    fp <- tryCatch(eglearn(data2mpareto(Zw, p),  rholist=rho_pos, reg_method="glasso", complete_Gamma=TRUE), error=function(e) NULL)
    fn <- tryCatch(eglearn(data2mpareto(-Zw, p), rholist=rho_neg, reg_method="glasso", complete_Gamma=TRUE), error=function(e) NULL)
    if (is.null(fp) || is.null(fn)) next
    gp <- fp$graph[[1]]; gn <- fn$graph[[1]]
    cp <- gamma_to_chi(fp$Gamma[[1]]); cn <- gamma_to_chi(fn$Gamma[[1]])
    dimnames(cp) <- dimnames(cn) <- list(SHORT_NAMES, SHORT_NAMES)
    rows[[w]] <- data.frame(
      window = window_days, p = p, Window_n = w, Mid_date = mid,
      dens_neg = graph_density(gn), dens_pos = graph_density(gp),
      mod_neg  = graph_modularity(gn), mod_pos = graph_modularity(gp),
      chi_be_neg = cn["BTC","ETH"], chi_be_pos = cp["BTC","ETH"],
      glasso_gauss_dens = if (with_glasso) glasso_gauss_density(Zw) else NA_real_,
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

results_wb <- createWorkbook()

# =============================================================================
# MODULE A. THRESHOLD ROBUSTNESS
# =============================================================================
if (RUN_THRESHOLD) {
  message("\n=== Module A: threshold robustness ===")
  thr_list <- list()
  for (p in THRESHOLD_GRID) {
    message(sprintf("  p = %.2f : selecting rho ...", p))
    r <- select_rho(Z, p)
    message(sprintf("    rho_pos=%.3f rho_neg=%.3f ; running windows ...", r$rho_pos, r$rho_neg))
    thr_list[[as.character(p)]] <- run_dynamic_summary(Z, dates, BASE_WINDOW, p, r$rho_pos, r$rho_neg)
  }
  thr_df <- do.call(rbind, thr_list)
  addWorksheet(results_wb, "ThresholdRobustness"); writeData(results_wb, "ThresholdRobustness", thr_df)
  
  thr_summary <- thr_df %>% group_by(p) %>%
    summarise(across(c(dens_neg,dens_pos,mod_neg,mod_pos,chi_be_neg,chi_be_pos),
                     ~round(mean(.x, na.rm=TRUE),3)), .groups="drop")
  addWorksheet(results_wb, "ThresholdSummary"); writeData(results_wb, "ThresholdSummary", thr_summary)
  print(thr_summary)
  
  thr_long <- thr_df %>%
    pivot_longer(c(dens_neg,dens_pos,mod_neg,mod_pos,chi_be_neg,chi_be_pos),
                 names_to="key", values_to="value") %>%
    mutate(Tail = ifelse(grepl("neg$", key), "Lower (crash)", "Upper (rally)"),
           Metric = case_when(grepl("^dens", key) ~ "Edge density",
                              grepl("^mod",  key) ~ "Modularity",
                              TRUE ~ "BTC-ETH chi"),
           p = factor(p))
  p_thr <- ggplot(thr_long, aes(Mid_date, value, colour=p)) +
    geom_line(linewidth=1.0) +
    facet_grid(Metric ~ Tail, scales="free_y") +
    labs(title=t_if("Robustness to the exceedance threshold"),
         x=NULL, y=NULL, colour="Threshold p") +
    theme_pub(14) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2)))
  save_fig("robustness_threshold.png", p_thr, 12, 9)
}

# =============================================================================
# MODULE B. WINDOW-LENGTH ROBUSTNESS  (p fixed at BASE_P)
# =============================================================================
if (RUN_WINDOW) {
  message("\n=== Module B: window-length robustness ===")
  r0 <- select_rho(Z, BASE_P)             # penalty is window-length-independent
  message(sprintf("  rho_pos=%.3f rho_neg=%.3f", r0$rho_pos, r0$rho_neg))
  win_list <- list()
  for (wd in WINDOW_GRID) {
    message(sprintf("  window = %d days ; running ...", wd))
    win_list[[as.character(wd)]] <- run_dynamic_summary(
      Z, dates, wd, BASE_P, r0$rho_pos, r0$rho_neg,
      with_glasso = (wd == BASE_WINDOW))   # compute glasso once, at the base window
  }
  win_df <- do.call(rbind, win_list)
  addWorksheet(results_wb, "WindowRobustness"); writeData(results_wb, "WindowRobustness", win_df)
  
  win_summary <- win_df %>% group_by(window) %>%
    summarise(across(c(dens_neg,dens_pos,mod_neg,mod_pos,chi_be_neg,chi_be_pos),
                     ~round(mean(.x, na.rm=TRUE),3)), .groups="drop")
  addWorksheet(results_wb, "WindowSummary"); writeData(results_wb, "WindowSummary", win_summary)
  print(win_summary)
  
  win_long <- win_df %>%
    pivot_longer(c(dens_neg,dens_pos,mod_neg,mod_pos,chi_be_neg,chi_be_pos),
                 names_to="key", values_to="value") %>%
    mutate(Tail = ifelse(grepl("neg$", key), "Lower (crash)", "Upper (rally)"),
           Metric = case_when(grepl("^dens", key) ~ "Edge density",
                              grepl("^mod",  key) ~ "Modularity",
                              TRUE ~ "BTC-ETH chi"),
           window = factor(window))
  p_win <- ggplot(win_long, aes(Mid_date, value, colour=window)) +
    geom_line(linewidth=1.0) +
    facet_grid(Metric ~ Tail, scales="free_y") +
    labs(title=t_if("Robustness to the window length"),
         x=NULL, y=NULL, colour="Window (days)") +
    theme_pub(14) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2)))
  save_fig("robustness_window.png", p_win, 12, 9)
}

# =============================================================================
# MODULE C. ECONOMIC SIGNIFICANCE
#   All comparisons hold the marginal distributions fixed (empirical) and vary
#   only the dependence, so they isolate the effect of tail dependence.
# =============================================================================
if (RUN_ECON) {
  message("\n=== Module C: economic significance ===")
  R   <- cor(Z)
  MVN <- MASS::mvrnorm(N_SIM, mu=rep(0,d), Sigma=R)   # Gaussian copula sample
  U_g <- pnorm(MVN)                                   # margin-free uniforms
  
  ## C1. chi as a function of the quantile level (crash / lower tail) --------
  ## empirical chi(q) stays high; the Gaussian copula chi(q) decays toward 0.
  qs <- c(0.90, 0.95, 0.975, 0.99)
  pairs <- list(c("BTC","ETH"), c("BTC","TRX"))
  chi_q_rows <- list()
  for (pr in pairs) {
    ie <- match(pr, SHORT_NAMES)
    emp_full <- sapply(qs, function(q) graphicalExtremes::emp_chi(-Z[, ie], p=q)[1,2])
    gau_full <- sapply(qs, function(q) graphicalExtremes::emp_chi(-MVN[, ie], p=q)[1,2])
    chi_q_rows[[paste(pr, collapse="-")]] <- data.frame(
      Pair = paste(pr, collapse="-"), q = qs,
      emp_chi = round(emp_full,3), gaussian_chi = round(gau_full,3))
  }
  chi_q_df <- do.call(rbind, chi_q_rows)
  addWorksheet(results_wb, "ChiByQuantile"); writeData(results_wb, "ChiByQuantile", chi_q_df)
  
  # Long form purely for plotting: turns solid/dashed into a proper legend.
  chi_q_long <- chi_q_df %>%
    pivot_longer(c(emp_chi, gaussian_chi), names_to="model", values_to="chi") %>%
    mutate(model = ifelse(model == "emp_chi", "Empirical", "Gaussian copula"))
  p_chiq <- ggplot(chi_q_long, aes(q, chi, colour=Pair, linetype=model)) +
    geom_line(linewidth=1.1) +
    geom_point(size=2.2) +
    scale_linetype_manual(values=c("Empirical"="solid", "Gaussian copula"="dashed")) +
    labs(title=t_if("Extremal correlation vs quantile level (lower tail)"),
         x="Quantile level q", y=expression(hat(chi)(q)),
         colour="Pair", linetype=NULL) +
    theme_pub(15) +
    guides(colour   = guide_legend(override.aes = list(linewidth = 2)),
           linetype = guide_legend(override.aes = list(linewidth = 0.8)))
  save_fig("chi_by_quantile.png", p_chiq, 9.5, 6)
  
  ## C2. joint-exceedance counts (margin-free; depends only on the copula) ----
  joint_tab <- list()
  for (q in c(0.05, 0.01)) {
    emp_thr <- apply(Z, 2, quantile, probs=q)               # empirical lower q-quantile
    emp_cnt <- rowSums(sweep(Z, 2, emp_thr, "<"))           # # assets crashing per day
    gau_cnt <- rowSums(U_g < q)                             # Gaussian copula
    ks <- 1:d
    emp_p <- sapply(ks, function(k) mean(emp_cnt >= k))
    gau_p <- sapply(ks, function(k) mean(gau_cnt >= k))
    joint_tab[[as.character(q)]] <- data.frame(
      q = q, k = ks,
      empirical = round(emp_p, 5), gaussian = round(gau_p, 5),
      ratio = round(emp_p / pmax(gau_p, 1e-9), 2))
  }
  joint_df <- do.call(rbind, joint_tab)
  addWorksheet(results_wb, "JointExceedance"); writeData(results_wb, "JointExceedance", joint_df)
  
  # headline number: P(>= ceil(d/2)+1 assets crash together) at q = 0.05
  kk <- ceiling(d/2) + 1
  hl <- subset(joint_df, q==0.05 & k==kk)
  message(sprintf("  HEADLINE: P(>=%d of %d crash together | q=0.05): empirical %.3f vs Gaussian %.3f (%.1fx)",
                  kk, d, hl$empirical, hl$gaussian, hl$ratio))
  
  p_joint <- ggplot(subset(joint_df, q==0.05),
                    aes(k)) +
    geom_line(aes(y=empirical, colour="Empirical"), linewidth=1.0) +
    geom_line(aes(y=gaussian,  colour="Gaussian copula"), linewidth=1.0, linetype="dashed") +
    geom_point(aes(y=empirical, colour="Empirical"), size=2.2) +
    geom_point(aes(y=gaussian,  colour="Gaussian copula"), size=2.2) +
    scale_y_log10() +
    scale_colour_manual(values=c("Empirical"="#e76f51","Gaussian copula"="#457b9d")) +
    labs(title=t_if("Probability that k or more of the 13 assets crash on the same day"),
         x="k (assets crashing together)", y="P(count >= k)", colour=NULL) +
    theme_pub(15) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2)))
  save_fig("joint_exceedance.png", p_joint, 9, 6)
  
  ## C3. equal-weight portfolio tail (empirical vs Gaussian copula vs independence)
  ## All three share the empirical margins; only the dependence differs.
  to_margins <- function(Umat) sapply(seq_len(d), function(j)
    quantile(Z[,j], probs=Umat[,j], type=8))
  port_emp   <- rowMeans(Z)
  port_gauss <- rowMeans(to_margins(U_g))
  port_indep <- rowMeans(sapply(seq_len(d), function(j) sample(Z[,j], N_SIM, replace=TRUE)))
  
  tail_stats <- function(x) {
    loss <- -x
    v95 <- quantile(loss, 0.95); v99 <- quantile(loss, 0.99)
    c(VaR95=v95, ES95=mean(loss[loss>=v95]),
      VaR99=v99, ES99=mean(loss[loss>=v99]))
  }
  port_df <- rbind(
    data.frame(model="Empirical dependence",    t(round(tail_stats(port_emp),3))),
    data.frame(model="Gaussian copula",         t(round(tail_stats(port_gauss),3))),
    data.frame(model="Independence",            t(round(tail_stats(port_indep),3))))
  addWorksheet(results_wb, "PortfolioTail"); writeData(results_wb, "PortfolioTail", port_df)
  print(port_df)
  message(sprintf("  Portfolio ES99 inflation, empirical / Gaussian = %.2fx",
                  tail_stats(port_emp)["ES99"] / tail_stats(port_gauss)["ES99"]))
}

# =============================================================================
# MODULE D. FORMAL ASYMMETRY TEST  (within-window block bootstrap)
#   Tests chi^- > chi^+ within representative windows, avoiding the
#   pseudo-replication of testing across 98%-overlapping windows.
#   Also records the stability of the distinguishing (TRON-periphery) edges.
# =============================================================================
if (RUN_ASYMMETRY) {
  message("\n=== Module D: asymmetry bootstrap ===")
  r0 <- if (exists("r0")) r0 else select_rho(Z, BASE_P)
  starts <- seq(1, n_total - BASE_WINDOW + 1, by = STEP_DAYS)
  mids   <- dates[round((starts + starts + BASE_WINDOW - 1)/2)]
  targets <- as.Date(c("2022-06-01","2023-09-01","2024-06-01"))
  sel <- sapply(targets, function(t) which.min(abs(as.numeric(mids - t))))
  
  moving_block <- function(M, block=BLOCK_LEN) {
    n <- nrow(M); nb <- ceiling(n/block)
    s <- sample(1:(n-block+1), nb, replace=TRUE)
    idx <- as.vector(vapply(s, function(x) x:(x+block-1), integer(block)))[1:n]
    M[idx, , drop=FALSE]
  }
  offdiag_mean <- function(M) mean(M[upper.tri(M)])
  trx_pairs_mean <- function(M) { i <- which(SHORT_NAMES=="TRX"); mean(M[i, -i]) }
  
  asym_rows <- list()
  for (s in sel) {
    idx <- starts[s]:(starts[s] + BASE_WINDOW - 1)
    Zw  <- Z[idx, ]; mid <- mids[s]
    message(sprintf("  window mid %s : %d replicates ...", format(mid), B_BOOT))
    gap <- be_neg <- be_pos <- trx_neg <- core_neg <- numeric(B_BOOT)
    for (b in 1:B_BOOT) {
      Zb <- moving_block(Zw)
      fp <- tryCatch(eglearn(data2mpareto(Zb, BASE_P),  rholist=r0$rho_pos, reg_method="glasso", complete_Gamma=TRUE), error=function(e) NULL)
      fn <- tryCatch(eglearn(data2mpareto(-Zb, BASE_P), rholist=r0$rho_neg, reg_method="glasso", complete_Gamma=TRUE), error=function(e) NULL)
      if (is.null(fp) || is.null(fn)) { gap[b] <- NA; next }
      cp <- gamma_to_chi(fp$Gamma[[1]]); cn <- gamma_to_chi(fn$Gamma[[1]])
      dimnames(cp) <- dimnames(cn) <- list(SHORT_NAMES, SHORT_NAMES)
      gap[b]    <- offdiag_mean(cn) - offdiag_mean(cp)
      be_neg[b] <- cn["BTC","ETH"]; be_pos[b] <- cp["BTC","ETH"]
      trx_neg[b]  <- trx_pairs_mean(cn)
      core_neg[b] <- mean(cn[upper.tri(cn)][ !grepl("TRX", outer(SHORT_NAMES,SHORT_NAMES,paste,sep="-")[upper.tri(cn)]) ])
    }
    ci <- function(x) quantile(x, c(.025,.5,.975), na.rm=TRUE)
    asym_rows[[length(asym_rows)+1]] <- data.frame(
      Mid_date = mid,
      metric = c("chi_neg_minus_pos (panel mean)","BTC-ETH chi (lower)","BTC-ETH chi (upper)",
                 "TRON-pair chi (lower)","core-pair chi (lower)"),
      lo  = round(c(ci(gap)[1], ci(be_neg)[1], ci(be_pos)[1], ci(trx_neg)[1], ci(core_neg)[1]),3),
      med = round(c(ci(gap)[2], ci(be_neg)[2], ci(be_pos)[2], ci(trx_neg)[2], ci(core_neg)[2]),3),
      hi  = round(c(ci(gap)[3], ci(be_neg)[3], ci(be_pos)[3], ci(trx_neg)[3], ci(core_neg)[3]),3))
  }
  asym_df <- do.call(rbind, asym_rows)
  addWorksheet(results_wb, "AsymmetryBootstrap"); writeData(results_wb, "AsymmetryBootstrap", asym_df)
  print(asym_df)
}

# =============================================================================
# MODULE E. GAUSSIAN BENCHMARK  (frequentist graphical lasso vs BGGM density)
#   Shows the density gap is not an artefact of Bayesian-vs-penalised
#   estimation: a like-for-like L1 Gaussian graph is just as sparse as BGGM.
# =============================================================================
if (RUN_GAUSS_BENCH) {
  message("\n=== Module E: Gaussian benchmark (graphical lasso) ===")
  # glasso density per window comes from the base-window run if available,
  # else compute here on the fly.
  if (exists("win_df") && any(win_df$window == BASE_WINDOW & is.finite(win_df$glasso_gauss_dens))) {
    glasso_mean <- mean(win_df$glasso_gauss_dens[win_df$window==BASE_WINDOW], na.rm=TRUE)
    ext_neg <- mean(win_df$dens_neg[win_df$window==BASE_WINDOW], na.rm=TRUE)
    ext_pos <- mean(win_df$dens_pos[win_df$window==BASE_WINDOW], na.rm=TRUE)
  } else {
    starts <- seq(1, n_total - BASE_WINDOW + 1, by = STEP_DAYS)
    gl <- sapply(starts, function(s) glasso_gauss_density(Z[s:(s+BASE_WINDOW-1), ]))
    glasso_mean <- mean(gl, na.rm=TRUE); ext_neg <- NA; ext_pos <- NA
  }
  # BGGM density read from the main results workbook if present
  bggm_mean <- NA_real_
  main_wb <- file.path(OUT_DIR, "crypto_extreme_results.xlsx")
  if (file.exists(main_wb)) {
    gs <- tryCatch(read.xlsx(main_wb, "GlobalStats"), error=function(e) NULL)
    if (!is.null(gs)) {
      bggm_mean <- mean(gs$Edge_density[gs$Type == "BGGM Baseline"], na.rm=TRUE)
      if (is.na(ext_neg)) ext_neg <- mean(gs$Edge_density[gs$Type=="Negative Extreme"], na.rm=TRUE)
      if (is.na(ext_pos)) ext_pos <- mean(gs$Edge_density[gs$Type=="Positive Extreme"], na.rm=TRUE)
    }
  }
  bench_df <- data.frame(
    Network = c("Lower-tail extremal (HR)","Upper-tail extremal (HR)",
                "Gaussian, BGGM (Bayesian)","Gaussian, graphical lasso (frequentist L1)"),
    Mean_edge_density = round(c(ext_neg, ext_pos, bggm_mean, glasso_mean), 3))
  addWorksheet(results_wb, "GaussianBenchmark"); writeData(results_wb, "GaussianBenchmark", bench_df)
  print(bench_df)
}

# =============================================================================
# Save workbook
# =============================================================================
saveWorkbook(results_wb, file.path(OUT_DIR, "crypto_robustness_results.xlsx"), overwrite = TRUE)
message("\nDone. Workbook -> output/crypto_robustness_results.xlsx ; figures -> figures/")