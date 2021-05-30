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

# Fit robust subset selection models and run the mixed-integer solver
fit <- rss.fit(x, y, k = 0:p, h = n - n.c, k.mio = 0:p, h.mio = n - n.c)

# Extract model coefficients and generate predictions
coef(fit, k = p0, h = n - n.c)
predict(fit, x, k = p0, h = n - n.c)

# Plot coefficient profiles
plot(fit)
