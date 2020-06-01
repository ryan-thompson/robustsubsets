# Set simulation parameters
set.seed(1)
n <- 100
p <- 10
p0 <- 5
n.c <- 10

# Generate training data with mixture error
beta <- c(rep(1, p0), rep(0, p - p0))
X <- matrix(rnorm(n * p), n, p)
e <- rnorm(n, c(rep(10, n.c), rep(0, n - n.c)))
y <- X %*% beta + e

# Fit best/robust subset selection models
fit.bss <- bss(X, y, n.cores = 1)
fit.rss <- rss(X, y, n.cores = 1)

# Extract model coefficients
bss.beta <- coef(fit.bss)
rss.beta <- coef(fit.rss)

# Check estimation error
ee.bss <- norm(bss.beta - c(0, beta), '2')
ee.rss <- norm(rss.beta - c(0, beta), '2')
cat('Best subsets estimation error:', ee.bss, '\n')
cat('Robust subsets estimation error:', ee.rss, '\n')

# Plot coefficient profiles
plot(fit.rss, type = 'profile')
# Each facet corresponds to a different value of h

# Plot cross-validation results
plot(fit.rss, type = 'cv')
# Each line corresponds to a different value of h

# Generate test data
X.test <- matrix(rnorm(n * p), n, p)
e.test <- rnorm(n)
y.test <- X.test %*% beta + e.test

# Make model predictions (using best parameters from cv)
pred.bss <- predict(fit.bss, X.test)
pred.rss <- predict(fit.rss, X.test)

# Compute prediction error
pe.bss <- 1 / n * norm(y.test - pred.bss, '2') ^ 2
pe.rss <- 1 / n * norm(y.test - pred.rss, '2') ^ 2
cat('Best subsets prediction error:', pe.bss, '\n')
cat('Robust subsets prediction error:', pe.rss, '\n')

