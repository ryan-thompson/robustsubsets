# This file contains a number of auxiliary functions to support primary functions

# Hard thresholding operator
H <- \(x, keep) {
  x[!(seq_along(x) %in% order(abs(x), decreasing = TRUE)[0:keep])] <- 0
  return(x)
}

# Loss function
f <- \(x, y, beta, eta) 0.5 * norm(y - x %*% beta - eta, '2') ^ 2

# Robust scale
rob.scale <- \(x, center = TRUE, scale = TRUE) {
  if (center) x <- apply(x, 2, \(x) x - stats::median(x))
  if (scale) {
    if (center) {
      x <- apply(x, 2, \(x) x / stats::mad(x))
    } else {
      x <- apply(x, 2, \(x) x / stats::mad(x, center = 0))
    }
  }
  return(x)
}

# Non-robust scale
scale2 <- \(x, center = TRUE, scale = TRUE) {
  if (center) x <- apply(x, 2, \(x) x - base::mean(x))
  if (scale) {
    if (center) {
      x <- apply(x, 2, \(x) x / stats::sd(x))
    } else {
      x <- apply(x, 2, \(x) x / sqrt(sum(x ^ 2) / (length(x) - 1)))
    }
  }
  return(x)
}

# Trimmed mean square prediction error
tmspe <- \(x, alpha = 0.25) mean(utils::head(sort(x ^ 2), length(x) * (1 - alpha)))

# Mean square prediction error
mspe <- \(x) mean(x ^ 2)
