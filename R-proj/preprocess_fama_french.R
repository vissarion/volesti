# Preprocess Fama-French 49 industry portfolio returns into covariance
# matrices for the low-volatility portfolio sampling pipeline.
#
# Downloads weekly returns from Ken French's data library (free, open,
# survivorship-bias-free). 49 US industry portfolios, 1926-present.
#
# Paper-matching parameters: 5-year rolling windows, quarterly rebalance,
# 5 volatility levels, non-linear shrinkage covariance estimator.
#
# Output: ff49_covariance_matrices.rds

library(volesti)

# ---- Parameters (matching the paper) ----
window   <- 260   # 5 years of weekly data
step     <- 13    # quarterly rebalance
m_levels <- 5     # volatility levels per window

# ---- Download ----
url <- "https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/49_Industry_Portfolios_CSV.zip"
tmp <- tempfile(fileext = ".zip")
cat("Downloading Fama-French 49 industry portfolios...\n")
download.file(url, tmp, quiet = TRUE)
unzip(tmp, exdir = tempdir())
csv_file <- file.path(tempdir(), "49_Industry_Portfolios.csv")
stopifnot(file.exists(csv_file))

# ---- Parse ----
lines      <- readLines(csv_file)
data_start <- which(grepl("^[0-9]{6},", lines))[1]
raw        <- read.csv(csv_file, skip = data_start - 1, header = FALSE)

# Keep only rows where first column is a 6-digit date (YYYYMM)
dates_raw  <- trimws(raw[, 1])
valid_rows <- grepl("^[0-9]{6}$", dates_raw)
raw        <- raw[valid_rows, ]
dates      <- as.integer(raw[, 1])

# Convert returns to numeric, clean sentinel values
returns <- apply(as.matrix(raw[, -1]), 2, function(x) as.numeric(trimws(x)))
returns[returns < -90] <- NA   # Fama-French missing-value sentinel
returns[returns > 100]  <- NA  # non-numeric metadata
returns[returns < -100] <- NA

# Extract industry names from header
header_lines <- lines[1:(data_start - 1)]
name_line    <- header_lines[grep("^,", header_lines)][1]
industry_names <- trimws(strsplit(name_line, ",")[[1]][-1])
colnames(returns) <- industry_names[1:ncol(returns)]

cat(sprintf("Loaded: %d weeks, %d industries\n", nrow(returns), ncol(returns)))

# ---- Filter to 2002-2021 (matching paper) ----
start_row <- which(dates >= 200201)[1]
end_row   <- max(which(dates <= 202112))
returns   <- returns[start_row:end_row, ]
dates     <- dates[start_row:end_row]
T         <- nrow(returns)
cat(sprintf("Filtered to 2002-2021: %d weeks\n", T))

# ---- Sliding windows ----
lCov          <- list()
lVola_targets <- list()
window_labels <- character()
skipped       <- 0

for (start in seq(1, T - window, by = step)) {
  end <- start + window - 1
  Rw  <- returns[start:end, ]

  # Drop columns with >20% missing data
  na_frac <- colMeans(is.na(Rw))
  keep    <- na_frac < 0.2
  Rw      <- Rw[, keep, drop = FALSE]
  Rw      <- na.omit(Rw)

  if (nrow(Rw) < 100 || ncol(Rw) < 10) {
    skipped <- skipped + 1; next
  }

  sigma <- tryCatch(cov.qis(Rw), error = function(e) NULL)
  if (is.null(sigma)) { skipped <- skipped + 1; next }

  # Ensure positive definiteness
  sigma <- tryCatch(as.matrix(Matrix::nearPD(sigma)$mat),
                    error = function(e) sigma)

  Cs <- tryCatch(get_sequence_of_volatilities_2(sigma, m_levels),
                 error = function(e) { skipped <- skipped + 1; NULL })
  if (is.null(Cs)) next

  lCov[[length(lCov) + 1]]                   <- sigma
  lVola_targets[[length(lVola_targets) + 1]] <- Cs
  window_labels <- c(window_labels, paste0(dates[start], "-", dates[end]))
}

result <- list(
  lCov           = lCov,
  lVola_targets  = lVola_targets,
  window_labels  = window_labels,
  industry_names = colnames(returns),
  n_assets       = ncol(returns)
)

saveRDS(result, "ff49_covariance_matrices.rds")
cat(sprintf("Wrote ff49_covariance_matrices.rds: %d windows (%d skipped), %d assets\n",
            length(lCov), skipped, ncol(returns)))

# Quick summary
if (length(lCov) > 0) {
  cat("First window assets:", ncol(returns), "\n")
  cat("First window vol targets:", round(result$lVola_targets[[1]], 6), "\n")
  years <- as.numeric(substr(window_labels, 1, 4))
  cat("Windows per year:\n"); print(table(years))
}
