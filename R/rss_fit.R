#==================================================================================================#
# This is the master function for robust subset selection. It calls mio and pbgd to perform
# mixed-integer optimisation and projected block-coordinate gradient descent, respectively.
#==================================================================================================#

#' @title Fit a robust subset selection model
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Fits a sequence of robust subset selection models using a combination of heuristics
#' and mixed-integer optimisation (mio).
#'
#' @param X a matrix of predictors
#' @param y a vector of the response
#' @param k the number of predictors to minimise sum of squares over (i.e. the model sparsity); by
#' default a sequence from 0 to 20
#' @param h the number of observations to minimise sum of squares over; by default a sequence from
#' 75 to 100 percent of sample size (in increments of 5 percent)
#' @param int a logical indicating whether to include an intercept
#' @param k.mio the subset of \code{k} for which the mixed-integer solver should be run
#' @param h.mio the subset of \code{h} for which the mixed-integer solver should be run
#' @param time a time limit in seconds on each call to the mixed-integer solver
#' @param tau a positive number greater than 1 used to tighten variable bounds in the mixed-integer
#' formulation; small values give quicker run times but can also exclude the optimal solution
#' @param output a logical indicating whether to print status updates
#' @param params an optional list of additional Gurobi parameters (the parameters Time and OutputFlag
#' are controlled by \code{time} and \code{output})
#' @param robust a logical indicating whether to standardise the data robustly; median/mad for true
#' and mean/sd for false
#' @param max.iter.ns the maximum number of neighbourhood search iterations to perform; if
#' output is true then the number of iterations required for convergence will be printed
#' @param max.iter.gd the maximum number of gradient descent iterations to perform
#' @param tol a numerical tolerance parameter used to declare convergence
#' @param ... any other arguments
#'
#' @return An object of class \code{rss.fit}; a list with the following components:
#' \item{beta}{an array of estimated regression coefficients; each column of regression
#'  coefficients corresponds to fixed value of \code{k} and each matrix to fixed value of \code{h}}
#' \item{weights}{an array of binary weights; weights equal to one correspond to good observations
#' selected for inclusion in the least squares fit; each column of weights corresponds to fixed
#' value of \code{k} and each matrix to fixed value of \code{h}}
#' \item{objval}{a matrix with the objective function values; each row corresponds to a value for
#' different \code{k} and each column to a value for different \code{h}}
#' \item{k}{the value of \code{k} that was passed in}
#' \item{h}{the value of \code{h} that was passed in}
#' \item{int}{whether an intercept was included}
#'
#' @details The function first computes solutions over all combinations of \code{k} and \code{h}
#'  using heuristics. The solutions can then be refined further using the mixed-integer solver.
#' The parameters that the solver operates on are specified by the \code{k.mio} and \code{h.mio}
#' parameters, which must be subsets of \code{k} and \code{h}. \cr
#' If robust is set to true and the median of any predictor is zero, then the data cannot be
#' standardised (the median absolute deviation is undefined) and an error message will be returned.
#'
#' @example R/examples/example_rss_fit.R
#'
#' @export
#'
#' @importFrom Matrix "Matrix"
#' @importFrom Matrix "bdiag"
#' @importFrom Matrix "Diagonal"

rss.fit <- function(X, y,
                    k = 0:min(nrow(X) - int, ncol(X), 20), h = round(seq(0.75, 1, 0.05) * nrow(X)),
                    int = T, k.mio = NULL, h.mio = NULL, time = 60, tau = 1.25, output = F,
                    params = NULL, robust = T, max.iter.ns = 1e2, max.iter.gd = 1e5, tol = 1e-4,
                    ...) {

  # Preliminaries
  X <- as.matrix(X)
  y <- as.matrix(y)
  n <- nrow(X)
  p <- ncol(X)
  n.k <- length(k)
  n.h <- length(h)
  n.k.mio <- length(k.mio)
  n.h.mio <- length(h.mio)

  # Centre y and centre/scale X; median/mad for robust and mean/sd for non-robust
  X.o <- X
  y.o <- y
  if (robust) {
    X <- rob.scale(X, center = int, scale = T)
    if (any(is.na(X))) stop('At least one predictor has MAD = 0 causing standardisation to fail.')
    y <- rob.scale(y, center = int, scale = F)
  } else {
    X <- scale(X, center = int, scale = T)
    y <- scale(y, center = int, scale = F)
  }

  # Run neighbourhood search
  if (output) cat ('Running neighbourhood search... \n \n')
  step <- 1 / norm(X, '2') ^ 2
  results <- list(x = array(dim = c(p + n + p + n, n.k, n.h)), objval = array(dim = c(n.k, n.h)))
  for (j in 1:n.h) {
    for (i in 1:n.k) {
      result <- pbgd(X, y, rep(0, p), rep(0, n), k[i], h[j], step, max.iter.gd, tol)
      results$x[, i, j] <- result$x
      results$objval[i, j] <- result$objval
    }
  }
  if (output) cat('Iteration:', 0, 'Total objectives:', sum(results$objval), '\n')
  for (iter in 1:max.iter.ns) {
    old <- results$objval
    for (j in 1:n.h) {
      for (i in 1:n.k) {
        for (a in (i - 1):(i + 1)) {
          if (a < 1 | a > n.k) next
          for (b in (j - 1):(j + 1)) {
            if (abs(i - a) + abs(j - b) != 1 | b < 1 | b > n.h)  next
            beta.init <- H(results$x[1:p, a, b], k[i])
            eta.init <- H(results$x[(p + 1):(p + n), a, b], n - h[j])
            result <- pbgd(X, y, beta.init, eta.init, k[i], h[j], step, max.iter.gd, tol)
            if (result$objval < results$objval[i, j]) {
              results$x[, i, j] <- result$x
              results$objval[i, j] <- result$objval
            }
          }
        }
      }
    }
    if (output) cat('Iteration:', iter, 'Total objectives:', sum(results$objval), '\n')
    new <- results$objval
    if (sum(old - new) <= tol) break
  }

  # Run mio
  if (!any(is.null(k.mio)) & !any(is.null(h.mio))) {
    if (output) cat('\nRunning mixed-integer optimisation... \n \n')
    for (j in 1:n.h.mio) {
      for (i in 1:n.k.mio) {
        if (k.mio[i] == 0) next
        result <- mio(X, y, k.mio[i], h.mio[j],
                      results$x[, which(k == k.mio[i]), which(h == h.mio[j])],
                      time, tau, output, params)
        results$x[, which(k == k.mio[i]), which(h == h.mio[j])] <- result$x[1:(p + n + p + n)]
        results$objval[which(k == k.mio[i]), which(h == h.mio[j])] <- result$objval
      }
    }
  }

  # Refit the models on the original (unscaled) data
  results.final <- list(beta = array(dim = c(p + int, n.k, n.h)),
                        weights = array(dim = c(n, n.k, n.h)),
                        objval = array(dim = c(n.k, n.h)), k = k, h = h, int = int)
  for (j in 1:n.h) {
    for (i in 1:n.k) {
      if (k[i] != 0) {
        beta <- results$x[1:p, i, j]
        eta <- results$x[(p + 1):(p + n), i, j]
        beta[abs(beta) < 1e-5] <- 0
        eta[abs(eta) < 1e-5] <- 0
        id.beta <- 1:p %in% which(beta != 0)
        id.eta <- 1:n %in% which(eta != 0)
        if(int) {
          beta <- c(0, beta)
          beta[c(T, id.beta)] <- stats::lsfit(X.o[!id.eta, id.beta], y.o[!id.eta], int = T)$coef
        } else {
          beta[id.beta] <- stats::lsfit(X.o[!id.eta, id.beta], y.o[!id.eta], int = F)$coef
        }
      } else {
        if (int) {
          eta <- H(y.o, n - h[j])
          id.eta <- 1:n %in% which(eta != 0)
          beta <- matrix(0, p + 1)
          beta[1] <- mean(y.o[!id.eta])
        } else {
          id.eta <- rep(F, n)
          beta <- matrix(0, p)
        }
      }
      if (int) eta <- y.o - cbind(1, X.o) %*% beta else eta <- y.o - X.o %*% beta
      eta[!id.eta] <- 0
      if (int) objval <- f(cbind(1, X.o), y.o, beta, eta) else objval <- f(X.o, y.o, beta, eta)
      results.final$beta[, i, j] <- beta
      results.final$weights[, i, j] <- as.numeric(eta == 0)
      results.final$objval[i, j] <- objval
    }
  }

  class(results.final) <- 'rss.fit'
  return(results.final)

}

#==================================================================================================#
# This function is responsible for mixed-integer optimisation. It contains the Gurobi specifications
# for solving the robust subset selection problem. The mio solver is capable of solving the problem
# to global optimality.

# Based on Ryan Tibshirani's Gurobi implementation of best subsets at
# https://github.com/ryantibs/best-subset
#==================================================================================================#

mio <- function(X, y, k, h, init, time, tau, output, params) {

  # Preliminaries
  n <- nrow(X)
  p <- ncol(X)
  W <- cbind(X, diag(1, n))
  form <- ifelse(ncol(X) <= nrow(X) & h == nrow(X), 1, 2)

  # Set bounds for variables
  Mb <- tau * max(abs(init[1:p]))
  Me <- tau * max(abs(init[(p + 1):(p + n)]))

  # Set the solver problem formulation
  model <- list()

  # Problem formulation 1
  # -------------------------------------------------------------
  # Objective function: x ^ T Q x + c ^ T x + a
  # Problem variable: x = [beta, eta, s, z]
  # Constraints
  # C1:      sum(s) <= k     <==> ||beta||_0 <= k
  # C2:      sum(z) <= n - h <==>  ||eta||_0 <= n - h
  # C3:   beta - Mb <= 0     <==>       beta <= Mb
  # C4: - beta - Mb <= 0     <==>       beta >= - Mb
  # C5:    eta - Me <= 0     <==>        eta <= Me
  # C6:  - eta - Me <= 0     <==>        eta >= - Me

  if (form == 1) {
    model$vtype <- c(rep('C', p + n), rep('B', p + n)) # Problem variable x
    model$Q <- 0.5 * bdiag(t(W) %*% W, Matrix(0, p + n, p + n)) # Matrix Q in objective function
    model$obj <- c(- t(W) %*% y, rep(0, p + n)) # Vector c in objective function
    model$objcon <- 0.5 * t(y) %*% y # Constant a in objective function
    model$A <- # LHS of constraints
      rbind(
        cbind(Matrix(0, 1, p), Matrix(0, 1, n), Matrix(1, 1, p), Matrix(0, 1, n)), # C1
        cbind(Matrix(0, 1, p), Matrix(0, 1, n), Matrix(0, 1, p), Matrix(1, 1, n)), # C2
        cbind(Diagonal(p, 1), Matrix(0, p, n), Diagonal(p, - Mb), Matrix(0, p, n)), # C3
        cbind(Diagonal(p, - 1), Matrix(0, p, n), Diagonal(p, - Mb), Matrix(0, p, n)), # C4
        cbind(Matrix(0, n, p), Diagonal(n, 1), Matrix(0, n, p), Diagonal(n, - Me)), # C5
        cbind(Matrix(0, n, p), Diagonal(n, - 1), Matrix(0, n, p), Diagonal(n, - Me)) # C6
      )
    model$sense <- c('<=', '<=', rep('<=', 2 * p), rep('<=', 2 * n)) # Constraint types
    model$rhs <- c(k, n - h, rep(0, 2 * p), rep(0, 2 * n)) # RHS of constraints
    model$lb <- c(rep(- Mb, p), rep(- Me, n), rep(0, p + n)) # Lower bound on variables
    model$ub <- c(rep(Mb, p), rep(Me, n), rep(1, p + n)) # Upper bound on variables
    model$start <- init  # Warm start
  }

  # Problem formulation 2
  # -------------------------------------------------------------
  # Objective function: x ^ T Q x + c ^ T x + a
  # Problem variable: x = [beta, eta, s, z, xi]
  # Constraints
  # C1:            sum(s) <= k     <==>   ||beta||_0 <= k
  # C2:            sum(z) <= n - h <==>    ||eta||_0 <= n - h
  # C3:         beta - Mb <= 0     <==>         beta <= Mb
  # C4:       - beta - Mb <= 0     <==>         beta >= - Mb
  # C5:          eta - Me <= 0     <==>          eta <= Me
  # C6:        - eta - Me <= 0     <==>          eta >= - Me
  # C7: W [beta, eta] - xi = 0     <==> W [beta, eta] = xi

  if (form == 2) {
    Mx <- tau * max(abs(W %*% init[1:(p + n)]))
    model$vtype <- c(rep('C', p + n), rep('B', p + n), rep('C', n)) # Problem variable x
    model$Q <- 0.5 * bdiag(Matrix(0, p + n + p + n, p + n + p + n), Diagonal(n, 1)) # Matrix Q in objective function
    model$obj <- c(- t(W) %*% y, rep(0, p + n + n)) # Vector c in objective function
    model$objcon <- 0.5 * t(y) %*% y # Constant a in objective function
    model$A <- # LHS of constraints
      rbind(
        cbind(Matrix(0, 1, p), Matrix(0, 1, n), Matrix(1, 1, p), Matrix(0, 1, n), Matrix(0, 1, n)), # C1
        cbind(Matrix(0, 1, p), Matrix(0, 1, n), Matrix(0, 1, p), Matrix(1, 1, n), Matrix(0, 1, n)), # C2
        cbind(Diagonal(p, 1), Matrix(0, p, n), Diagonal(p, - Mb), Matrix(0, p, n), Matrix(0, p, n)), # C3
        cbind(Diagonal(p, - 1), Matrix(0, p, n), Diagonal(p, - Mb), Matrix(0, p, n), Matrix(0, p, n)), # C4
        cbind(Matrix(0, n, p), Diagonal(n, 1), Matrix(0, n, p), Diagonal(n, - Me), Matrix(0, n, n)), # C5
        cbind(Matrix(0, n, p), Diagonal(n, - 1), Matrix(0, n, p), Diagonal(n, - Me), Matrix(0, n, n)), # C6
        cbind(W, Matrix(0, n, p), Matrix(0, n, n), Diagonal(n, - 1)) # C7
      )
    model$sense <- c('<=', '<=', rep('<=', 2 * p), rep('<=', 2 * n), rep('=', n)) # Constraint types
    model$rhs <- c(k, n - h, rep(0, 2 * p), rep(0, 2 * n), rep(0, n)) # RHS of constraints
    model$lb <- c(rep(- Mb, p), rep(- Me, n), rep(0, p + n), rep(- Mx, n)) # Lower bound on variables
    model$ub <- c(rep(Mb, p), rep(Me, n), rep(1, p + n), rep(Mx, n)) # Upper bound on variables
    model$start <- c(init, W %*% init[1:(p + n)])  # Warm start
  }

  # Set the solver parameters
  if (is.null(params)) params <- list()
  params$TimeLimit <- time
  params$OutputFlag <- ifelse(output, 1, 0)

  # Solve the model using Gurobi
  result <- gurobi::gurobi(model, params)

  return(result)

}

#==================================================================================================#
# Coefficient function
#==================================================================================================#

#' @title Coefficient function for rss.fit object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Extracts coefficients for a given parameter pair \code{(k,h)}.
#'
#' @param object an object of class \code{rss.fit}
#' @param k the number of predictors indexing the desired fit
#' @param h the number of observations indexing the desired fit
#' @param ... any other arguments
#'
#' @return An array of coefficients.
#'
#' @method coef rss.fit
#'
#' @export
#'
#' @importFrom stats "coef"

coef.rss.fit <- function(object, k = NULL, h = NULL, ...) {

  if (!is.null(k)) index1 <- which(object$k == k) else index1 <- 1:length(object$k)
  if (!is.null(h)) index2 <- which(object$h == h) else index2 <- 1:length(object$h)
  object$beta[, index1, index2]

}

#==================================================================================================#
# Predict function
#==================================================================================================#

#' @title Predict function for rss.fit object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Generate predictions for new data using a given parameter pair \code{(k,h)}.
#'
#' @param object an object of class \code{rss.fit}
#' @param X.new a matrix of new values for the predictors
#' @param k the number of predictors indexing the desired fit
#' @param h the number of observations indexing the desired fit
#' @param ... any other arguments
#'
#' @return An array of predictions.
#'
#' @method predict rss.fit
#'
#' @export
#'
#' @importFrom stats "predict"

predict.rss.fit <- function(object, X.new, k = NULL, h = NULL, ...) {

  X.new <- as.matrix(X.new)
  if (object$int) X.new <- cbind(1, X.new)
  beta <- coef.rss.fit(object, k, h, ...)
  if (length(dim(beta)) < 3) X.new %*% beta else apply(beta, 2:3, function(beta) X.new %*% beta)

}

#==================================================================================================#
# Plot function
#==================================================================================================#

#' @title Plot function for rss.fit object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Plot the coefficient profiles from robust subset selection.
#'
#' @param x an object of class \code{rss.fit}
#' @param ... any other arguments
#'
#' @return A plot of the coefficient profiles.
#'
#' @method plot rss.fit
#'
#' @export
#'
#' @importFrom graphics "plot"

plot.rss.fit <- function(x, ...) {

  beta <- x$beta
  if (x$int) beta <- beta[- 1, , , drop = F]
  beta[beta == 0] <- NA
  beta <- data.frame(beta = as.vector(beta),
                     predictor = as.factor(rep(1:nrow(beta), length(x$k) * length(x$h))),
                     k = as.factor(rep(rep(x$k, each = nrow(beta), length(x$h)))),
                     h = as.factor(rep(x$h, each = nrow(beta) * length(x$k))))
  beta <- stats::na.omit(beta)
  ggplot2::ggplot(beta, ggplot2::aes_string('k', 'beta', col = 'predictor')) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(. ~ h, labeller = ggplot2::label_both)

}
