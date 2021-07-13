


# robustsubsets

## Overview

An R implementation of [robust subset
selection](https://arxiv.org/abs/2005.08217).

Robust subset selection is a robust adaption of the classic best subset
selection estimator, and is defined by the constrained least squares
problem:

![](man/figures/README-equation.png)<!-- -->

Robust subsets seeks out the best subset of predictors and observations
and performs a least squares fit on this subset. The number of
predictors used in the fit is controlled by the parameter `k` and the
observations by the parameter `h`.

## Installation

You should install Gurobi and the associated R package gurobi before
installing robustsubsets. Gurobi is available for free under academic
license at <https://www.gurobi.com/>.

To install `robustsubsets` from GitHub, run the following code:

``` r
devtools::install_github('ryan-thompson/robustsubsets')
```

## Usage

The `rss()` function fits a robust subset regression model for a grid of
`k` and `h`. The `cv.rss()` function provides a convenient way to
automatically cross-validate these parameters.

``` r
library(robustsubsets)

# Generate training data with contaminated predictor matrix
set.seed(123)
n <- 100 # Number of observations
p <- 10 # Number of predictors
p0 <- 5 # Number of relevant predictors
ncontam <- 5 # Number of contaminated observations
beta <- c(rep(1, p0), rep(0, p - p0))
x <- matrix(rnorm(n * p), n, p)
y <- x %*% beta + rnorm(n)
x[1:ncontam, ] <- matrix(rnorm(ncontam * p, mean = 10), ncontam, p)

# Fit the robust subset selection regularisation path
fit <- rss(x, y)
coef(fit, k = p0, h = n - ncontam)
```

    ##  [1] 0.1719763 1.0787905 0.9532380 0.8548021 1.1739899 1.0471684 0.0000000
    ##  [8] 0.0000000 0.0000000 0.0000000 0.0000000

``` r
# Cross-validate the robust subset selection regularisation path
fit <- cv.rss(x, y)
coef(fit)
```

    ##  [1] 0.1328087 1.0782082 0.9546673 0.8761584 1.0892592 1.0549823 0.0000000
    ##  [8] 0.2406216 0.0000000 0.0000000 0.0000000

## Documentation

See the package [reference manual](robustsubsets_1.1.0.pdf).
