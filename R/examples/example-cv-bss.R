# Generate training data
set.seed(0)
n <- 100
p <- 10
p0 <- 5
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n)
y <- x %*% beta + e

# Best subset selection with cross-validation
cl <- parallel::makeCluster(2)
fit <- cv.bss(x, y, cluster = cl)
parallel::stopCluster(cl)

# Extract model coefficients, generate predictions, and plot cross-validation results
coef(fit)
predict(fit, x[1:3, ])
plot(fit)
