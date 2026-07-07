# Preprocess MSCI country index returns into covariance matrices for the
# low-volatility portfolio sampling pipeline.
#
# Reads ../../backtest-light/data/msci_country_indices.csv, applies a sliding
# window, estimates covariance per window via QIS shrinkage, and saves the
# results as .rds files consumable by run_real_data_small.R / run_real_data_large.R.
#
# Output:
#   msci_covariance_matrices_small.rds  — G7 countries (7 assets)
#   msci_covariance_matrices_large.rds  — all 24 countries

library(volesti)

window   <- 252
step     <- 63
m_levels <- 5

data_path <- "../../backtest-light/data/msci_country_indices.csv"
stopifnot(file.exists(data_path))
raw     <- read.csv(data_path)
returns <- as.matrix(raw[, -1])
T       <- nrow(returns)
all_cc  <- colnames(returns)

# ---- Small (G7) ----
cc  <- c("US", "GB", "JP", "DE", "FR", "CA", "IT")
idx <- match(cc, all_cc)
R   <- returns[, idx]

lCov_small          <- list()
lVola_targets_small <- list()

for (start in seq(1, T - window, by = step)) {
  end <- start + window - 1
  Rw  <- R[start:end, ]
  Rw  <- Rw[rowSums(abs(Rw)) > 0, ]
  if (any(colSums(abs(Rw)) < 1e-12)) next
  sigma <- cov.qis(Rw)
  Cs    <- get_sequence_of_volatilities_2(sigma, m_levels)
  lCov_small[[length(lCov_small) + 1]]                         <- sigma
  lVola_targets_small[[length(lVola_targets_small) + 1]]       <- Cs
}
saveRDS(list(lCov = lCov_small, lVola_targets = lVola_targets_small),
        "msci_covariance_matrices_small.rds")
cat(sprintf("Small: %d windows, %d assets\n", length(lCov_small), length(cc)))

# ---- Large (all countries) ----
R <- returns

lCov_large          <- list()
lVola_targets_large <- list()

for (start in seq(1, T - window, by = step)) {
  end <- start + window - 1
  Rw  <- R[start:end, ]
  Rw  <- Rw[rowSums(abs(Rw)) > 0, ]
  if (any(colSums(abs(Rw)) < 1e-12)) next
  sigma <- cov.qis(Rw)
  Cs    <- get_sequence_of_volatilities_2(sigma, m_levels)
  lCov_large[[length(lCov_large) + 1]]                         <- sigma
  lVola_targets_large[[length(lVola_targets_large) + 1]]       <- Cs
}
saveRDS(list(lCov = lCov_large, lVola_targets = lVola_targets_large),
        "msci_covariance_matrices_large.rds")
cat(sprintf("Large: %d windows, %d assets\n", length(lCov_large), length(all_cc)))
