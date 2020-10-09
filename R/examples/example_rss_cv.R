# Generate training data with mixture error
set.seed(1)
n <- 100
p <- 10
p0 <- 5
n.c <- 5
beta <- c(rep(1, p0), rep(0, p - p0))
X <- matrix(rnorm(n * p), n, p)
e <- rnorm(n, c(rep(10, n.c), rep(0, n - n.c)))
y <- X %*% beta + e

# Cross-validate robust subset selection models
cv <- rss.cv(X, y, k = 0:10, h = function(n) round(c(0.95, 1.00) * n), n.core = 1)

# Plot cross-validation results
plot(cv)
