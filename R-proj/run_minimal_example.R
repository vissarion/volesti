# Minimal reproducible example for volesti portfolio sampling
# Generates synthetic data and samples constant-volatility long-only portfolios

library(volesti)

set.seed(42)

# --- Step 1: Create a synthetic covariance matrix (5 assets) ---
n <- 5
A <- matrix(rnorm(n * 3), nrow = 3, ncol = n)
sigma <- t(A) %*% A + diag(0.1, n)  # positive definite
sigma <- cov2cor(sigma)              # correlation-like scale

# --- Step 2: Set a target volatility level ---
c <- 0.3  # target portfolio variance

# --- Step 3: Sample portfolios at constant volatility ---
M <- 2000  # number of portfolios
result <- sample_ptfs_constant_volatility(sigma, c, M, ignore_smallest_components = TRUE)

# --- Step 4: Inspect results ---
samples <- result$overall_samples[[1]]
cat("Sampled", ncol(samples), "portfolios of dimension", nrow(samples), "\n")

# Verify: each column sums to 1 (long-only simplex constraint)
col_sums <- colSums(samples)
cat("Column sums (should all be 1):", range(col_sums), "\n")

# Verify: empirical volatility matches target c
portfolio_returns <- t(samples) %*% sigma %*% samples
empirical_var <- mean(diag(portfolio_returns))
cat("Target volatility:", c, "\n")
cat("Empirical mean variance:", empirical_var, "\n")
cat("First 3 portfolios:\n")
print(samples[, 1:min(3, ncol(samples))])
