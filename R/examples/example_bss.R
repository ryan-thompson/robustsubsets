# Generate training data
set.seed(1)
n <- 100
p <- 10
p0 <- 5
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n)
y <- x %*% beta + e

# Fit best subset selection models
# Run the mixed-integer solver on the k that minimises the cv error
fit <- bss(x, y, k = 0:10, mio = 'min', n.core = 1)

# Extract model coefficients and generate predictions
coef(fit)
predict(fit, x)

# Plot coefficient profiles and cross-validation results
plot(fit, type = 'profile')
plot(fit, type = 'cv')
