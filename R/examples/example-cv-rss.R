# Generate training data with mixture error
set.seed(123)
n <- 100
p <- 10
p0 <- 5
ncontam <- 5
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n, c(rep(10, ncontam), rep(0, n - ncontam)))
y <- x %*% beta + e

# Robust subset selection with cross-validation
fit <- cv.rss(x, y)

# Extract model coefficients, generate predictions, and plot cross-validation results
coef(fit)
predict(fit, x)
plot(fit)
