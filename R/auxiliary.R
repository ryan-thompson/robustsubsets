# This file contains a number of auxiliary functions to support primary functions

# Hard thresholding operator
H <- function(x, keep) {
  x[!(1:length(x) %in% order(abs(x), decreasing = T)[0:keep])] <- 0
  return(x)
}

# Loss function
f <- function(X, y, beta, eta) 0.5 * norm(y - X %*% beta - eta, '2') ^ 2

# Robust scale
rob.scale <- function(X, center = T, scale = T) {
  if (center) X <- apply(X, 2, function(x) x - stats::median(x))
  if (scale) {
    if (center) {
      X <- apply(X, 2, function(x) x / stats::mad(x))
    } else {
      X <- apply(X, 2, function(x) x / stats::mad(x, center = 0))
    }
  }
  return(X)
}

# Trimmed mean square prediction error
tmspe <- function(x, alpha = 0.25) mean(sort(x ^ 2)[1:round(length(x) * (1 - alpha))])

# Mean square prediction error
mspe <- function(x) mean(x ^ 2)
