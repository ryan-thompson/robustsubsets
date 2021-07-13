# Generate training data with mixture error
set.seed(0)
n <- 100
p <- 10
p0 <- 5
ncontam <- 10
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n, c(rep(10, ncontam), rep(0, n - ncontam)))
y <- x %*% beta + e

# Robust subset selection with cross-validation
cl <- parallel::makeCluster(2)
fit <- cv.rss(x, y, cluster = cl)
parallel::stopCluster(cl)

# Extract model coefficients, generate predictions, and plot cross-validation results
coef(fit)
predict(fit, x[1:3, ])
plot(fit)
