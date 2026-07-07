# Minimal reproducible example for volesti portfolio sampling
# Generates synthetic data and samples constant-volatility long-only portfolios

library(volesti)
set.seed(42)

# ---- Step 1: Create two synthetic covariance matrices (5 assets) ----
n <- 5
A1 <- matrix(rnorm(n * 3), 3, n)
sigma_in  <- cov2cor(t(A1) %*% A1 + diag(0.1, n))   # period t (in-sample)

A2 <- matrix(rnorm(n * 3, mean = 0.1), 3, n)
sigma_out <- cov2cor(t(A2) %*% A2 + diag(0.1, n))   # period t+1 (out-of-sample)

# ---- Step 2: Set a target volatility level ----
c <- 0.3  # target portfolio variance

# ---- Step 3: Sample portfolios at constant volatility ----
M <- 500  # number of portfolios
result <- sample_ptfs_constant_volatility(sigma_in, c, M)
samples <- result$overall_samples[[1]]
cat("Sampled", ncol(samples), "portfolios of dimension", nrow(samples), "\n\n")

# ---- Step 4: Verify constraints ----

# Long-only: each column sums to 1
col_sums <- colSums(samples)
cat("Column sums (should all be 1):", range(col_sums), "\n")

# In-sample variance: exactly c (hard geometric constraint)
var_in <- mean(diag(t(samples) %*% sigma_in %*% samples))
cat(sprintf("In-sample variance:  %.8f  (target: %.1f) ← exact by construction\n",
            var_in, c))

# Out-of-sample variance: will differ (covariance changed next period)
var_out <- mean(diag(t(samples) %*% sigma_out %*% samples))
cat(sprintf("Out-of-sample var:  %.6f  (target: %.1f) ← estimation error\n",
            var_out, c))

# Number of disconnected components found
cat(sprintf("Connected components: %d\n", length(result$overall_samples)))

# ---- Step 5: Show a few portfolios ----
cat("\nFirst 3 portfolios:\n")
print(round(samples[, 1:min(3, ncol(samples))], 4))
