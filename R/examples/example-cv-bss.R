# Generate training data
set.seed(123)
n <- 100
p <- 10
p0 <- 5
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n)
y <- x %*% beta + e

# Best subset selection with cross-validation
fit <- cv.bss(x, y)

# Extract model coefficients, generate predictions, and plot cross-validation results
coef(fit)
predict(fit, x)
plot(fit)
