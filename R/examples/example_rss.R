# Generate training data with mixture error
set.seed(1)
n <- 100
p <- 10
p0 <- 5
n.c <- 5
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
e <- rnorm(n, c(rep(10, n.c), rep(0, n - n.c)))
y <- x %*% beta + e

# Fit robust subset selection models
# Run the mixed-integer solver on the (k,h) that minimises the cv error
fit <- rss(x, y, k = 0:10, h = function(n) round(c(0.95, 1.00) * n), mio = 'min', n.core = 1)

# Extract model coefficients and generate predictions
coef(fit)
predict(fit, x)

# Plot coefficient profiles and cross-validation results
plot(fit, type = 'profile')
plot(fit, type = 'cv')
