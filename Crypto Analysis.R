# =============================================================================
# crypto_extreme_ggm_pipeline_v18.R
# Sliding-window Hüsler-Reiss extreme graphical models for crypto markets
#
# v18 changes (FIGURES ONLY; all statistics identical to v17):
# - Shared publication theme theme_pub() with larger fonts throughout.
# - Human-readable axis labels (no more "Mid_date", "Edge_density", etc.).
# - Thicker lines, larger points, bigger and clearer legends.
# - Higher output resolution (dpi 200) and tuned figure dimensions.
# - Larger text on the base-R ACF panels.
# - SHOW_TITLES toggle: in-figure titles/subtitles are OFF by default,
#   since the LaTeX captions already carry them. Set TRUE for standalone use.
# =============================================================================

# ---- 0. Packages ------------------------------------------------------------
required_pkgs <- c(
  "quantmod", "rugarch", "graphicalExtremes", "igraph",
  "ergm", "network", "intergraph", "openxlsx", "ggplot2",
  "dplyr", "tidyr", "purrr", "lubridate", "patchwork",
  "xts", "BGGM", "tseries", "FinTS", "moments", "av", "cluster"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

set.seed(42)

# ---- 0b. Figure styling -----------------------------------------------------
# Master switch: keep in-figure titles off for the paper (captions cover them),
# flip to TRUE if you want titles when viewing the PNGs on their own.
SHOW_TITLES <- FALSE
FIG_DPI     <- 200

# Helper so each plot can request a title that only appears when SHOW_TITLES.
t_if <- function(x) if (SHOW_TITLES) x else NULL

# Shared publication theme: large, legible, minimal.
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

# Consistent saver.
save_fig <- function(name, plot, width, height) {
  ggsave(file.path(FIG_DIR, name), plot, width = width, height = height,
         dpi = FIG_DPI, bg = "white")
}

# ---- 1. Configuration -------------------------------------------------------

COINS <- c("BTC-USD", "ETH-USD", "XRP-USD", "BNB-USD", "SOL-USD",
           "ADA-USD", "DOGE-USD", "TRX-USD", "LINK-USD", "DOT-USD",
           "AVAX-USD", "LTC-USD", "XLM-USD")

SHORT_NAMES <- c("BTC", "ETH", "XRP", "BNB", "SOL",
                 "ADA", "DOGE", "TRX", "LINK", "DOT",
                 "AVAX", "LTC", "XLM")

SECTORS <- c(
  BTC = "Currency", LTC = "Currency", DOGE = "Currency",
  ETH = "SmartContract", SOL = "SmartContract",
  ADA = "SmartContract", AVAX = "SmartContract",
  XRP = "Payment", XLM = "Payment", TRX = "Payment",
  BNB = "Infrastructure", LINK = "Infrastructure", DOT = "Infrastructure"
)

START_DATE <- as.Date("2020-10-01")
END_DATE <- Sys.Date()

# Window configuration
WINDOW_DAYS <- 750L     
STEP_DAYS <- 15L        
EXTREME_P <- 0.20       

# Grid for Full-Sample CV
RHO_GRID <- seq(0.01, 0.15, length.out = 15) 

OUT_DIR <- "output"
FIG_DIR <- "figures"
dir.create(OUT_DIR, showWarnings = FALSE)
dir.create(FIG_DIR, showWarnings = FALSE)

EVENTS <- data.frame(
  date = as.Date(c("2021-05-19", "2022-05-09", "2022-11-11",
                   "2024-01-11", "2024-11-13")),
  label = c("BTC Crash", "Terra/Luna", "FTX Collapse", "BTC ETF", "BTC ATH"),
  stringsAsFactors = FALSE
)

ERGM_SEED <- 42


# ---- 2. Download & align daily prices ---------------------------------------

download_crypto <- function(coins, from, to) {
  price_list <- list()
  for (i in seq_along(coins)) {
    ticker <- coins[i]; nm <- SHORT_NAMES[i]
    message(sprintf("  Downloading %s ...", ticker))
    tryCatch({
      dat <- getSymbols(ticker, src = "yahoo", from = from, to = to, auto.assign = FALSE)
      price_list[[nm]] <- Ad(dat)
    }, error = function(e) warning(sprintf("Failed: %s", ticker)))
  }
  prices <- do.call(merge, price_list)
  colnames(prices) <- names(price_list)
  prices <- na.locf(prices, na.rm = FALSE)
  na.omit(prices)
}

message("=== [1/8] Downloading price data ===")
prices <- download_crypto(COINS, START_DATE, END_DATE)
log_returns <- na.omit(diff(log(prices))[-1, ])


# ============================================================================
# ---- 3. GJR-GARCH filtering (FULL DATASET) ---------------------------------
# ============================================================================

message("\n=== [2/8] GJR-GARCH(1,1)-t filtering on full dataset ===")

garch_filter <- function(returns_xts) {
  d <- ncol(returns_xts)
  nms <- colnames(returns_xts)
  n <- nrow(returns_xts)
  resid_mat <- matrix(NA, n, d, dimnames = list(NULL, nms))
  
  spec_gjr <- ugarchspec(
    variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = "std"
  )
  spec_fallback <- ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(1, 0), include.mean = TRUE),
    distribution.model = "std"
  )
  
  for (j in seq_len(d)) {
    coin <- nms[j]
    fit <- tryCatch(
      ugarchfit(spec_gjr, data = as.numeric(returns_xts[, j]), solver = "hybrid"),
      error = function(e) ugarchfit(spec_fallback, data = as.numeric(returns_xts[, j]), solver = "hybrid")
    )
    resid_mat[, j] <- as.numeric(residuals(fit, standardize = TRUE))
  }
  list(residuals = resid_mat, dates = index(returns_xts))
}

filtered <- garch_filter(log_returns)
Z <- filtered$residuals
dates <- filtered$dates


# ============================================================================
# ---- 4. Residual diagnostics & Plots ---------------------------------------
# ============================================================================

message("\n=== [3/8] Residual diagnostics & Plots ===")

d <- ncol(Z); nms <- colnames(Z)

# 1. ACF of residuals  (base R; larger text + higher resolution)
png(file.path(FIG_DIR, "acf_residuals.png"), width = 1700, height = 1500, res = 150)
par(mfrow = c(4, 4), mar = c(3.4, 3.4, 2.6, 1),
    cex.main = 1.6, cex.axis = 1.15, cex.lab = 1.2, mgp = c(2.1, 0.7, 0))
for (j in seq_len(d)) acf(Z[, j], lag.max = 30, main = nms[j], col = "steelblue", lwd = 1.6)
dev.off()

# 2. ACF of squared residuals
png(file.path(FIG_DIR, "acf_squared_residuals.png"), width = 1700, height = 1500, res = 150)
par(mfrow = c(4, 4), mar = c(3.4, 3.4, 2.6, 1),
    cex.main = 1.6, cex.axis = 1.15, cex.lab = 1.2, mgp = c(2.1, 0.7, 0))
for (j in seq_len(d)) acf(Z[, j]^2, lag.max = 30, main = paste0(nms[j], "\u00b2"), col = "tomato", lwd = 1.6)
dev.off()

# 3. QQ plots
qq_data <- do.call(rbind, lapply(seq_len(d), function(j) {
  x <- sort(Z[, j]); n <- length(x)
  data.frame(Theoretical = qnorm(ppoints(n)), Sample = x, Coin = nms[j])
}))
p_qq <- ggplot(qq_data, aes(x = Theoretical, y = Sample)) +
  geom_point(size = 0.6, alpha = 0.5, colour = "steelblue") +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 0.9) +
  facet_wrap(~ Coin, ncol = 4) +
  labs(title = t_if("QQ vs Normal — standardised residuals"),
       x = "Theoretical quantile", y = "Sample quantile") +
  theme_pub(14)
save_fig("qq_plots.png", p_qq, 13, 11)

# 4. Density plots
dens_data <- do.call(rbind, lapply(seq_len(d), function(j) {
  de <- density(Z[, j], n = 512)
  data.frame(x = de$x, y = de$y, Coin = nms[j])
}))
x_seq <- seq(-6, 6, length.out = 300)
norm_df <- data.frame(x = x_seq, y = dnorm(x_seq))
p_dens <- ggplot(dens_data, aes(x = x, y = y)) +
  geom_area(fill = "steelblue", alpha = 0.35) +
  geom_line(colour = "steelblue", linewidth = 0.8) +
  geom_line(data = norm_df, aes(x = x, y = y), colour = "red", linetype = "dashed", linewidth = 0.9, inherit.aes = FALSE) +
  facet_wrap(~ Coin, ncol = 4) +
  coord_cartesian(xlim = c(-6, 6)) +
  labs(title = t_if("Kernel density vs N(0,1)"),
       x = "Standardised residual", y = "Density") +
  theme_pub(14)
save_fig("density_plots.png", p_dens, 13, 11)

# 5. Hill plots
hill_estimate <- function(x, k_max = NULL) {
  x_pos <- sort(abs(x), decreasing = TRUE)
  n <- length(x_pos)
  if (is.null(k_max)) k_max <- floor(n / 2)
  k_seq <- 2:k_max
  alpha <- sapply(k_seq, function(k) 1 / mean(log(x_pos[1:k]) - log(x_pos[k + 1])))
  data.frame(k = k_seq, alpha = alpha)
}
hill_data <- do.call(rbind, lapply(seq_len(d), function(j) {
  rbind(
    data.frame(hill_estimate(Z[, j][Z[, j] > 0]), Tail = "Right", Coin = nms[j]),
    data.frame(hill_estimate(-Z[, j][Z[, j] < 0]), Tail = "Left", Coin = nms[j])
  )
}))
p_hill <- ggplot(hill_data, aes(x = k, y = alpha, colour = Tail)) +
  geom_line(linewidth = 0.9, alpha = 0.9) +
  geom_hline(yintercept = 2, colour = "red", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 4, colour = "orange", linetype = "dashed", linewidth = 0.8) +
  facet_wrap(~ Coin, ncol = 4) +
  scale_colour_manual(values = c(Right = "steelblue", Left = "tomato")) +
  labs(title = t_if("Hill plots — tail index"),
       x = "Number of upper order statistics (k)", y = expression(hat(alpha)),
       colour = "Tail") +
  theme_pub(14) +
  guides(colour = guide_legend(override.aes = list(linewidth = 2)))
save_fig("hill_plots.png", p_hill, 13, 11)

# Statistical test table
diag_rows <- lapply(seq_len(d), function(j) {
  x <- Z[, j]; n <- length(x)
  lb_z <- Box.test(x, lag = 10, type = "Ljung-Box")
  lb_z2 <- Box.test(x^2, lag = 10, type = "Ljung-Box")
  arch <- tryCatch(FinTS::ArchTest(x, lags = 5)$p.value, error = function(e) NA)
  jb <- tseries::jarque.bera.test(x)$p.value
  k_hill <- floor(n * 0.10)
  x_sort <- sort(abs(x), decreasing = TRUE)
  alpha_h <- 1 / mean(log(x_sort[1:k_hill]) - log(x_sort[k_hill + 1]))
  
  data.frame(
    Coin = nms[j], n = n,
    Excess_kurtosis = round(moments::kurtosis(x) - 3, 2),
    Hill_alpha_10pct = round(alpha_h, 3),
    LB_Z_p10 = round(lb_z$p.value, 4),
    LB_Z2_p10 = round(lb_z2$p.value, 4),
    ARCHLM_p5 = round(arch, 4),
    JB_p = round(jb, 6),
    Pass_LB_Z2 = lb_z2$p.value > 0.05,
    Pass_ARCH = !is.na(arch) && arch > 0.05,
    stringsAsFactors = FALSE
  )
})
diag_df <- do.call(rbind, diag_rows)

wb_diag <- createWorkbook()
addWorksheet(wb_diag, "Diagnostics"); writeData(wb_diag, "Diagnostics", diag_df)
saveWorkbook(wb_diag, file.path(OUT_DIR, "residual_diagnostics.xlsx"), overwrite = TRUE)

pca_out <- prcomp(Z, scale. = FALSE)
eig_vals <- pca_out$sdev^2
var_exp <- round(eig_vals / sum(eig_vals) * 100, 1)


# ============================================================================
# ---- 5. Global Penalty Selection via 10-Fold CV (Full Sample) --------------
# ============================================================================

message("\n=== [4/8] Global Penalty Selection (10-Fold CV on Full Sample) ===")

Y_pos_full <- data2mpareto(Z, p = EXTREME_P)
Y_neg_full <- data2mpareto(-Z, p = EXTREME_P)

K <- 10
n_pos <- nrow(Y_pos_full)
n_neg <- nrow(Y_neg_full)

folds_pos  <- sample(rep(1:K, length.out = n_pos))
folds_neg  <- sample(rep(1:K, length.out = n_neg))

cv_loglik_pos <- matrix(NA, length(RHO_GRID), K)
cv_loglik_neg <- matrix(NA, length(RHO_GRID), K)

for (k in 1:K) {
  # Positive Extremes
  train_idx_pos <- which(folds_pos != k); test_idx_pos  <- which(folds_pos == k)
  Y_train_pos <- Y_pos_full[train_idx_pos, , drop = FALSE]; Y_test_pos  <- Y_pos_full[test_idx_pos, , drop = FALSE]
  
  fit_pos <- tryCatch(eglearn(data = Y_train_pos, rholist = RHO_GRID, reg_method = "glasso", complete_Gamma = TRUE), error = function(e) NULL)
  if (!is.null(fit_pos)) {
    for (j in seq_along(RHO_GRID)) {
      if (!is.null(fit_pos$Gamma[[j]])) cv_loglik_pos[j, k] <- tryCatch(loglik_HR(data = Y_test_pos, Gamma = fit_pos$Gamma[[j]], graph = fit_pos$graph[[j]])[1], error = function(e) NA)
    }
  }
  
  # Negative Extremes
  train_idx_neg <- which(folds_neg != k); test_idx_neg  <- which(folds_neg == k)
  Y_train_neg <- Y_neg_full[train_idx_neg, , drop = FALSE]; Y_test_neg  <- Y_neg_full[test_idx_neg, , drop = FALSE]
  
  fit_neg <- tryCatch(eglearn(data = Y_train_neg, rholist = RHO_GRID, reg_method = "glasso", complete_Gamma = TRUE), error = function(e) NULL)
  if (!is.null(fit_neg)) {
    for (j in seq_along(RHO_GRID)) {
      if (!is.null(fit_neg$Gamma[[j]])) cv_loglik_neg[j, k] <- tryCatch(loglik_HR(data = Y_test_neg, Gamma = fit_neg$Gamma[[j]], graph = fit_neg$graph[[j]])[1], error = function(e) NA)
    }
  }
}

get_best_rho <- function(mat, grid) {
  avg <- suppressWarnings(rowMeans(mat, na.rm = TRUE))
  valid_idx <- which(!is.na(avg) & !is.nan(avg))
  if (length(valid_idx) == 0) return(grid[round(length(grid)/2)]) 
  best_idx <- valid_idx[which.max(avg[valid_idx])]
  return(grid[best_idx])
}

GLOBAL_RHO_POS <- get_best_rho(cv_loglik_pos, RHO_GRID)
GLOBAL_RHO_NEG <- get_best_rho(cv_loglik_neg, RHO_GRID)

message(sprintf("  -> Optimal Pos rho: %.4f", GLOBAL_RHO_POS))
message(sprintf("  -> Optimal Neg rho: %.4f", GLOBAL_RHO_NEG))


# ============================================================================
# ---- 6. Helper functions for Networks & ERGM -------------------------------
# ============================================================================

compute_global_stats <- function(g, label) {
  g_uw <- g; E(g_uw)$weight <- rep(1, ecount(g_uw))
  data.frame(
    Window = label, Nodes = vcount(g), Edges = ecount(g),
    Edge_density = round(edge_density(g), 4),
    Avg_degree = round(mean(degree(g)), 2),
    Diameter = ifelse(is.connected(g_uw), diameter(g_uw), NA_integer_),
    Avg_path_length = round(mean_distance(g_uw, directed = FALSE), 3),
    Global_transitivity = round(transitivity(g, type = "global"), 3),
    Louvain_modularity = round(modularity(cluster_louvain(g_uw)), 3),
    Degree_assortativity = round(assortativity_degree(g), 3),
    stringsAsFactors = FALSE
  )
}

compute_centrality <- function(g, label) {
  g_uw <- g; E(g_uw)$weight <- rep(1, ecount(g_uw))
  data.frame(
    Window = label, Node = V(g)$name, Sector = SECTORS[V(g)$name],
    Degree = degree(g_uw), Betweenness = round(betweenness(g_uw), 2),
    Closeness = round(closeness(g_uw), 4),
    Eigenvector = round(eigen_centrality(g_uw)$vector, 4),
    PageRank = round(page.rank(g_uw)$vector, 4),
    stringsAsFactors = FALSE
  )
}

fit_ergm_mple <- function(g, lbl) {
  ne <- ecount(g); nv <- vcount(g)
  max_edges <- nv * (nv - 1) / 2
  if (ne < 1 || ne >= max_edges) return(NULL)
  
  net <- tryCatch(intergraph::asNetwork(g), error = function(e) NULL)
  if (is.null(net)) return(NULL)
  set.vertex.attribute(net, "sector", SECTORS[network.vertex.names(net)])
  
  ctrl <- control.ergm(seed = ERGM_SEED)
  tryCatch({
    m0 <- ergm(net ~ edges, estimate = "MPLE", control = ctrl)
    m1 <- ergm(net ~ edges + nodefactor("sector"), estimate = "MPLE", control = ctrl)
    m2 <- ergm(net ~ edges + nodefactor("sector") + nodematch("sector") + gwesp(decay = 0.5, fixed = TRUE),
               estimate = "MPLE", control = ctrl)
    list(
      label = lbl, m0 = m0, m1 = m1, m2 = m2,
      aic = c(M0 = AIC(m0), M1 = AIC(m1), M2 = AIC(m2)),
      bic = c(M0 = BIC(m0), M1 = BIC(m1), M2 = BIC(m2))
    )
  }, error = function(e) NULL)
}

# Hüsler-Reiss theoretical mapping from Gamma to Chi
gamma_to_chi <- function(Gamma_mat) {
  chi <- 2 - 2 * pnorm(sqrt(Gamma_mat) / 2)
  diag(chi) <- 1
  return(chi)
}

# ============================================================================
# ---- 7. Sliding-window loop ------------------------------------------------
# ============================================================================

n_total <- nrow(Z)
starts <- seq(1, n_total - WINDOW_DAYS + 1, by = STEP_DAYS)
n_windows <- length(starts)
message(sprintf("\n=== [5/8] Sliding window: %d windows, %d days, step %d ===",
                n_windows, WINDOW_DAYS, STEP_DAYS))

all_global <- list(); all_centrality <- list(); all_rho <- list()
all_chi_pos <- list(); all_chi_neg <- list(); all_ergm <- list()
gamma_trajectories <- list()

for (w in seq_len(n_windows)) {
  idx <- starts[w]:(starts[w] + WINDOW_DAYS - 1)
  Z_w <- Z[idx, ]
  d_w <- dates[idx]
  mid <- d_w[round(length(d_w) / 2)]
  label <- sprintf("W%02d [%s]", w, format(mid, "%Y-%m-%d"))
  message(sprintf("  Window %02d / %02d [%s]", w, n_windows, format(mid, "%Y-%m-%d")))
  
  Y_pos_w <- data2mpareto(Z_w, p = EXTREME_P)
  Y_neg_w <- data2mpareto(-Z_w, p = EXTREME_P)
  
  all_rho[[w]] <- data.frame(
    Window_n = w, Mid_date = mid, Label = label,
    rho_pos = GLOBAL_RHO_POS, rho_neg = GLOBAL_RHO_NEG,
    stringsAsFactors = FALSE
  )
  
  fp <- tryCatch(eglearn(data = Y_pos_w, rholist = GLOBAL_RHO_POS, reg_method = "glasso", complete_Gamma = TRUE), error=function(e) NULL)
  fn <- tryCatch(eglearn(data = Y_neg_w, rholist = GLOBAL_RHO_NEG, reg_method = "glasso", complete_Gamma = TRUE), error=function(e) NULL)
  
  if (is.null(fp) || is.null(fn)) next
  
  active_coins <- colnames(Z_w)
  g_pos <- fp$graph[[1]]; V(g_pos)$name <- active_coins
  g_neg <- fn$graph[[1]]; V(g_neg)$name <- active_coins
  
  # Model-implied chi from estimated Gamma
  chi_pos <- gamma_to_chi(fp$Gamma[[1]]); rownames(chi_pos) <- colnames(chi_pos) <- active_coins
  chi_neg <- gamma_to_chi(fn$Gamma[[1]]); rownames(chi_neg) <- colnames(chi_neg) <- active_coins
  all_chi_pos[[w]] <- chi_pos; all_chi_neg[[w]] <- chi_neg
  
  # Assign edge widths proportional to chi
  set_edge_width <- function(g, chi_mat) {
    el <- as_edgelist(g, names = FALSE)
    wts <- numeric(nrow(el))
    for(e in seq_len(nrow(el))) wts[e] <- chi_mat[el[e,1], el[e,2]]
    E(g)$width <- wts * 2
    g
  }
  g_pos <- set_edge_width(g_pos, chi_pos)
  g_neg <- set_edge_width(g_neg, chi_neg)
  
  # Extract Gamma entries
  extract_gamma <- function(Gam, tail_name) {
    if (is.null(Gam)) return(NULL)
    rownames(Gam) <- colnames(Gam) <- active_coins
    ut <- upper.tri(Gam)
    data.frame(
      Window_n = w, Mid_date = mid, Tail = tail_name,
      Pair = paste(rownames(Gam)[row(Gam)[ut]], colnames(Gam)[col(Gam)[ut]], sep = "-"),
      Gamma = Gam[ut], stringsAsFactors = FALSE
    )
  }
  gamma_trajectories[[length(gamma_trajectories) + 1]] <- rbind(
    extract_gamma(fp$Gamma[[1]], "Positive"),
    extract_gamma(fn$Gamma[[1]], "Negative")
  )
  
  # BGGM baseline
  fit_est <- tryCatch(BGGM::estimate(Z_w, type = "continuous"), error=function(e) NULL)
  fit_sel <- tryCatch(BGGM::select(fit_est), error=function(e) NULL)
  
  if (!is.null(fit_sel)) {
    adj <- fit_sel$adj; diag(adj) <- 0
    pcor <- fit_sel$pcor
    rownames(pcor) <- colnames(pcor) <- active_coins
    
    g_gl <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    V(g_gl)$name <- active_coins
    
    el <- as_edgelist(g_gl, names = FALSE)
    wts <- numeric(nrow(el))
    for(e in seq_len(nrow(el))) wts[e] <- pcor[el[e,1], el[e,2]]
    E(g_gl)$weight <- wts
  } else {
    g_gl <- make_empty_graph(n = length(active_coins), directed = FALSE); V(g_gl)$name <- active_coins
  }
  
  all_global[[w]] <- rbind(
    compute_global_stats(g_pos, paste(label, "Pos")),
    compute_global_stats(g_neg, paste(label, "Neg")),
    compute_global_stats(g_gl, paste(label, "BGGM"))
  )
  all_centrality[[w]] <- rbind(
    compute_centrality(g_pos, paste(label, "Pos")),
    compute_centrality(g_neg, paste(label, "Neg")),
    compute_centrality(g_gl, paste(label, "BGGM"))
  )
  
  all_ergm[[w]] <- list(
    pos = fit_ergm_mple(g_pos, paste(label, "Pos")),
    neg = fit_ergm_mple(g_neg, paste(label, "Neg")),
    bggm = fit_ergm_mple(g_gl, paste(label, "BGGM"))
  )
  
  # Save Networks  (movie frames; larger node labels and titles)
  for (info in list(
    list(g = g_pos, suf = "pos", col = "palegreen3"),
    list(g = g_neg, suf = "neg", col = "tomato"),
    list(g = g_gl, suf = "bggm", col = "white")
  )) {
    fname <- file.path(FIG_DIR, sprintf("net_W%02d_%s.png", w, info$suf))
    png(fname, width = 900, height = 900, res = 120)
    par(cex.main = 1.4)
    lay <- layout_in_circle(info$g)
    
    if (info$suf == "bggm") {
      plot(info$g, layout = lay, main = sprintf("%s – %s", label, toupper(info$suf)),
           vertex.size = 24, vertex.label.cex = 1.05, vertex.label.color = "black",
           vertex.label.family = "sans", vertex.label.font = 2,
           vertex.color = "white", vertex.frame.color = "gray40",
           edge.color = ifelse(E(info$g)$weight > 0, "darkgreen", "darkorange"),
           edge.width = abs(E(info$g)$weight) * 2)
    } else {
      plot(info$g, layout = lay, main = sprintf("%s – %s", label, toupper(info$suf)),
           vertex.size = 24, vertex.label.cex = 1.05, vertex.label.color = "black",
           vertex.label.family = "sans", vertex.label.font = 2,
           vertex.color = info$col, vertex.frame.color = "black",
           edge.color = "gray40", edge.width = E(info$g)$width)
    }
    dev.off()
  }
}


# ============================================================================
# ---- 8. Compile results & temporal plots -----------------------------------
# ============================================================================

message("\n=== [6/8] Compiling results & plots ===")

global_df <- do.call(rbind, all_global)
centrality_df <- do.call(rbind, all_centrality)
rho_df <- do.call(rbind, all_rho)
gamma_df <- do.call(rbind, gamma_trajectories)

global_df <- global_df %>%
  mutate(
    Window_n = as.integer(gsub("W(\\d+).*", "\\1", Window)),
    Type = case_when(
      grepl("Pos$", Window) ~ "Positive Extreme",
      grepl("Neg$", Window) ~ "Negative Extreme",
      grepl("BGGM$", Window) ~ "BGGM Baseline",
      TRUE ~ "Other"
    )
  ) %>%
  left_join(rho_df %>% dplyr::select(Window_n, Mid_date), by = "Window_n")

# Pretty y-axis label per metric (no raw column names on the figure).
plot_metric <- function(df, metric, panel_title, y_lab) {
  ggplot(df %>% filter(Type != "Other"),
         aes(x = Mid_date, y = .data[[metric]], colour = Type, linetype = Type)) +
    geom_line(linewidth = 1.0) + geom_point(size = 1.8) +
    geom_vline(data = EVENTS, aes(xintercept = date), linetype = "dashed", colour = "grey55", alpha = 0.6) +
    geom_text(data = EVENTS, aes(x = date, y = Inf, label = label), angle = 90,
              hjust = 1.05, vjust = -0.4, size = 3.2, colour = "grey35", inherit.aes = FALSE) +
    scale_colour_manual(values = c("Positive Extreme" = "#2a9d8f", "Negative Extreme" = "#e76f51", "BGGM Baseline" = "#457b9d")) +
    labs(title = panel_title, x = NULL, y = y_lab) +
    theme_pub(14) +
    theme(legend.title = element_blank()) +
    guides(colour = guide_legend(override.aes = list(linewidth = 2, size = 2.4)))
}

p_ev <- (plot_metric(global_df, "Edge_density", "Network density", "Edge density") +
           plot_metric(global_df, "Global_transitivity", "Clustering", "Transitivity")) /
  (plot_metric(global_df, "Louvain_modularity", "Community modularity", "Modularity") +
     plot_metric(global_df, "Avg_degree", "Average degree", "Average degree")) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
p_ev <- p_ev + plot_annotation(
  title = t_if("Dynamic Conditional Extremal Dependence in Cryptocurrency Markets"),
  theme = theme(plot.title = element_text(face = "bold", size = 17, hjust = 0.5)))
save_fig("temporal_evolution.png", p_ev, 15, 10)

# Extract BTC and ETH chi
extract_chi <- function(target_coin) {
  do.call(rbind, lapply(seq_along(all_chi_pos), function(w) {
    if (is.null(all_chi_pos[[w]])) return(NULL)
    active_coins <- colnames(all_chi_pos[[w]])
    if (!target_coin %in% active_coins) return(NULL)
    others <- setdiff(active_coins, target_coin)
    rbind(
      data.frame(Window_n = w, Mid_date = rho_df$Mid_date[w], Coin = others, Sector = SECTORS[others],
                 Chi = all_chi_pos[[w]][target_coin, others], Tail = "Positive", stringsAsFactors = FALSE),
      data.frame(Window_n = w, Mid_date = rho_df$Mid_date[w], Coin = others, Sector = SECTORS[others],
                 Chi = all_chi_neg[[w]][target_coin, others], Tail = "Negative", stringsAsFactors = FALSE)
    )
  }))
}

btc_chi_long <- extract_chi("BTC")
eth_chi_long <- extract_chi("ETH")

plot_chi <- function(df, coin_name) {
  for (tdir in c("Positive", "Negative")) {
    tail_lab <- if (tdir == "Negative") "lower (crash)" else "upper (rally)"
    p_chi <- ggplot(df %>% filter(Tail == tdir), aes(x = Mid_date, y = Chi, colour = Coin)) +
      geom_line(linewidth = 1.0) + geom_point(size = 1.4) + facet_wrap(~ Sector) +
      geom_vline(data = EVENTS, aes(xintercept = date), linetype = "dashed", alpha = 0.5, colour = "grey55") +
      labs(title = t_if(sprintf("%s pairwise extremal dependence — %s tail", coin_name, tail_lab)),
           x = NULL, y = expression(hat(chi)), colour = NULL) +
      theme_pub(15) +
      guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                                   override.aes = list(linewidth = 2, size = 2.4)))
    ggsave(file.path(FIG_DIR, sprintf("%s_chi_%s.png", tolower(coin_name), tolower(tdir))),
           p_chi, width = 11, height = 8, dpi = FIG_DPI, bg = "white")
  }
}
plot_chi(btc_chi_long, "BTC")
plot_chi(eth_chi_long, "ETH")


# Gamma Clustering
gamma_long <- do.call(rbind, gamma_trajectories)
gamma_wide <- gamma_long %>% pivot_wider(names_from = Window_n, values_from = Gamma, id_cols = c(Pair, Tail))

cluster_assignments <- data.frame()
for (t_dir in c("Positive", "Negative")) {
  gw_sub <- gamma_wide %>% filter(Tail == t_dir)
  if(nrow(gw_sub) == 0) next
  mat <- as.matrix(gw_sub %>% dplyr::select(-Tail, -Pair))
  
  n_unq <- nrow(unique(mat))
  if (n_unq < 2) {
    gw_sub$Cluster <- as.factor(1)
  } else {
    max_k <- min(6, n_unq - 1)
    if (max_k < 2) {
      gw_sub$Cluster <- as.factor(1)
    } else {
      best_k <- 2; best_sil <- -1
      for (k in 2:max_k) {
        km <- tryCatch(kmeans(mat, centers = k, nstart = 25), error=function(e) NULL)
        if (!is.null(km)) {
          sil <- mean(cluster::silhouette(km$cluster, dist(mat))[, 3])
          if (sil > best_sil) { best_sil <- sil; best_k <- k }
        }
      }
      km_final <- kmeans(mat, centers = best_k, nstart = 25)
      gw_sub$Cluster <- as.factor(km_final$cluster)
    }
  }
  cluster_assignments <- rbind(cluster_assignments, gw_sub %>% dplyr::select(Tail, Pair, Cluster))
}

gamma_long <- gamma_long %>% left_join(cluster_assignments, by = c("Tail", "Pair"))
gamma_long$Mid_date <- as.Date(gamma_long$Mid_date)
# Friendly facet labels for the tails.
gamma_long$Tail <- factor(gamma_long$Tail, levels = c("Negative", "Positive"),
                          labels = c("Lower (crash) tail", "Upper (rally) tail"))

p_gamma <- ggplot(gamma_long, aes(x = Mid_date, y = Gamma, group = Pair, colour = Cluster)) +
  geom_line(alpha = 0.7, linewidth = 0.9) + facet_wrap(~ Tail, ncol = 1) +
  geom_vline(data = EVENTS, aes(xintercept = date), linetype = "dashed", alpha = 0.5, colour = "grey55") +
  labs(title = t_if("Pairwise Gamma Evolution (K-Means Clustered)"),
       x = NULL, y = expression(Gamma), colour = "Cluster") +
  theme_pub(15) +
  theme(legend.position = "right") +
  guides(colour = guide_legend(override.aes = list(linewidth = 2.2, alpha = 1)))
save_fig("gamma_evolution_clustered.png", p_gamma, 12, 8)

# ============================================================================
# ---- 10. Export to Excel & Generate Movies ---------------------------------
# ============================================================================

message("\n=== [8/8] Exporting to Excel & Generating Movies ===")

wb <- createWorkbook()
addWorksheet(wb, "GlobalStats"); writeData(wb, "GlobalStats", global_df)
addWorksheet(wb, "Centrality"); writeData(wb, "Centrality", centrality_df)
addWorksheet(wb, "RhoSelection"); writeData(wb, "RhoSelection", rho_df)
addWorksheet(wb, "ResidualDiag"); writeData(wb, "ResidualDiag", diag_df)
addWorksheet(wb, "BTC_chi"); writeData(wb, "BTC_chi", btc_chi_long)
addWorksheet(wb, "ETH_chi"); writeData(wb, "ETH_chi", eth_chi_long)
addWorksheet(wb, "GammaClusters"); writeData(wb, "GammaClusters", cluster_assignments)

saveWorkbook(wb, file.path(OUT_DIR, "crypto_extreme_results.xlsx"), overwrite = TRUE)

for (suffix in c("pos", "neg", "bggm")) {
  img_files <- sort(list.files(FIG_DIR, pattern = sprintf("^net_W[0-9]+_%s\\.png$", suffix), full.names = TRUE))
  if (length(img_files) > 0) {
    out_file <- file.path(FIG_DIR, sprintf("movie_evolution_%s.mp4", suffix))
    tryCatch(av_encode_video(img_files, output = out_file, framerate = 3), error = function(e) message("Video encoding failed."))
  }
}

message("\n=== Pipeline complete ===")
message(sprintf("  Results -> %s/", OUT_DIR))
message(sprintf("  Figures -> %s/", FIG_DIR))


# ============================================================================
# ---- 11. Sector-level Extremal Dependence (Average Chi) --------------------
# ============================================================================

message("\n=== Computing Sector-Level Average Chi ===")

sector_chi_list <- list()

for (w in seq_len(n_windows)) {
  mid_date <- rho_df$Mid_date[w]
  
  for (tdir in c("Positive", "Negative")) {
    chi_mat <- if (tdir == "Positive") all_chi_pos[[w]] else all_chi_neg[[w]]
    if (is.null(chi_mat)) next
    
    # Convert matrix to long format
    chi_long <- as.data.frame(as.table(chi_mat), stringsAsFactors = FALSE)
    colnames(chi_long) <- c("Coin1", "Coin2", "Chi")
    
    # Remove self-pairs and duplicates (keep only upper triangle pairs)
    chi_long <- chi_long[chi_long$Coin1 < chi_long$Coin2, ]
    
    # Map to sectors
    chi_long$Sector1 <- SECTORS[chi_long$Coin1]
    chi_long$Sector2 <- SECTORS[chi_long$Coin2]
    
    # Create standardized sector pair label (alphabetical order)
    chi_long$Sector_Pair <- apply(chi_long[, c("Sector1", "Sector2")], 1, function(x) {
      paste(sort(x), collapse = " - ")
    })
    
    # Determine if within or between sectors
    chi_long$Type <- ifelse(chi_long$Sector1 == chi_long$Sector2, 
                            "Within-Sector", "Between-Sector")
    
    # Aggregate mean chi
    agg_chi <- chi_long %>%
      group_by(Sector_Pair, Type) %>%
      summarise(Mean_Chi = mean(Chi, na.rm = TRUE), .groups = "drop") %>%
      mutate(Window_n = w, Mid_date = mid_date, Tail = tdir)
    
    sector_chi_list[[length(sector_chi_list) + 1]] <- agg_chi
  }
}

sector_chi_df <- do.call(rbind, sector_chi_list)
sector_chi_df$Mid_date <- as.Date(sector_chi_df$Mid_date)
# Friendly facet labels.
sector_chi_df$Tail <- factor(sector_chi_df$Tail, levels = c("Negative", "Positive"),
                             labels = c("Lower (crash) tail", "Upper (rally) tail"))

# Plotting
p_sector_chi <- ggplot(sector_chi_df, aes(x = Mid_date, y = Mean_Chi, colour = Sector_Pair)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.6) +
  facet_grid(Type ~ Tail) +
  geom_vline(data = EVENTS, aes(xintercept = date), linetype = "dashed", alpha = 0.5, colour = "grey55") +
  labs(title = t_if("Average Extremal Dependence by Sector Pairs"),
       subtitle = t_if("Macro-level structural shifts across crypto sectors"),
       x = NULL, y = expression(bar(chi)), colour = "Sector pair") +
  theme_pub(14) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(linewidth = 2, size = 2.4)))

save_fig("sector_average_chi.png", p_sector_chi, 14, 9)
message("  Saved: sector_average_chi.png")