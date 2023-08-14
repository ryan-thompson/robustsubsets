#==================================================================================================#
# This is the master function for robust subset selection. It calls heuristics and mio.
#==================================================================================================#

#' @title Robust subset selection
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Fits a sequence of regression models using robust subset selection.
#'
#' @param x a predictor matrix
#' @param y a response vector
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param h the number of observations to minimise sum of squares over; by default a sequence from
#' 75 to 100 percent of sample size (in increments of 5 percent)
#' @param k.mio the subset of \code{k} for which the mixed-integer solver should be run
#' @param h.mio the subset of \code{h} for which the mixed-integer solver should be run
#' @param params a list of parameters (settings) to pass to the mixed-integer solver (Gurobi)
#' @param tau a positive number greater than or equal to 1 used to tighten coefficient bounds in the
#' mixed-integer solver; small values give quicker run times but can also exclude the optimal
#' solution; can be \code{Inf}
#' @param warm.start a logical indicating whether to warm start the mio solver using the heuristics
#' @param robust a logical indicating whether to standardise the data robustly; median/mad for
#' \code{TRUE} and mean/sd for \code{FALSE}
#' @param max.ns.iter the maximum number of neighbourhood search iterations allowed
#' @param max.gd.iter the maximum number of gradient descent iterations allowed per value of
#' \code{k} and \code{h}
#' @param eps a numerical tolerance parameter used to declare convergence
#'
#' @return An object of class \code{rss}; a list with the following components:
#' \item{beta}{an array of estimated regression coefficients; columns correspond to \code{k} and
#' matrices to \code{h}}
#' \item{weights}{an array of binary weights; weights equal to one correspond to good observations
#' selected for inclusion in the least squares fit; columns correspond to \code{k} and matrices to
#' \code{h}}
#' \item{objval}{a matrix with the objective function values; rows correspond to \code{k} and
#' columns to \code{h}}
#' \item{mipgap}{a matrix with the optimality gap values; rows correspond to \code{k} and columns
#' to \code{h}}
#' \item{k}{a vector containing the values of \code{k} used in the fit}
#' \item{h}{a vector containing the values of \code{h} used in the fit}
#'
#' @details The function first computes solutions over all combinations of \code{k} and \code{h}
#' using heuristics. The heuristics include projected block-coordinate gradient descent and
#' neighbourhood search (see \href{https://arxiv.org/abs/2005.08217}{arXiv}). The solutions produced
#' by the heuristics can be refined further using the mixed-integer solver. The tuning parameters
#' that the solver operates on are specified by the \code{k.mio} and \code{h.mio} parameters,
#' which must be subsets of \code{k} and \code{h}. \cr
#'
#' By default, the mixed-integer optimization problem is formulated with SOS constraints and
#' bound constraints. The bound constraints are estimated as \eqn{\tau\|\hat{\beta}\|_\infty}, where
#' \eqn{\hat{\beta}} is output from the heuristics. For finite values of \code{tau}, the
#' mixed-integer solver automatically converts the SOS constraints to Big-M constraints, which are
#' more numerically efficient to optimise.
#'
#' @references Thompson, R. (2021). 'Robust subset selection'. arXiv:
#' \href{https://arxiv.org/abs/2005.08217}{2005.08217}.
#'
#' @example R/examples/example-rss.R
#'
#' @export

rss <- \(x, y, k = 0:min(nrow(x) - 1, ncol(x), 20), h = round(seq(0.75, 1, 0.05) * nrow(x)),
         k.mio = NULL, h.mio = NULL, params = list(TimeLimit = 60, OutputFlag = 0), tau = 1.5,
         warm.start = TRUE, robust = TRUE, max.ns.iter = 1e2, max.gd.iter = 1e5, eps = 1e-4) {

  # Check data is valid
  if (!is.matrix(x)) x <- as.matrix(x)
  attributes(x)$dimnames <- NULL
  if (!is.matrix(y)) y <- as.matrix(y)
  attributes(y)$dimnames <- NULL

  # Check arguments
  if (!is.null(k.mio)) if (!all(k.mio %in% k)) stop('k.mio not a subset of k')
  if (!is.null(h.mio)) if (!all(h.mio %in% h)) stop('h.mio not a subset of h')

  # Preliminaries
  n <- nrow(x)
  p <- ncol(x)
  nk <- length(k)
  nh <- length(h)
  nk.mio <- length(k.mio)
  nh.mio <- length(h.mio)

  # Centre y and centre/scale x; median/mad for robust and mean/sd for nonrobust
  x.o <- x
  y.o <- y
  if (robust) {
    x <- rob.scale(x, center = TRUE, scale = TRUE)
    if (any(is.na(x))) stop('At least one predictor has MAD = 0 causing robust standardisation to
                            fail. Set robust = FALSE or remove offending predictors.')
    y <- rob.scale(y, center = TRUE, scale = FALSE)
  } else {
    x <- scale(x, center = TRUE, scale = TRUE)
    y <- scale(y, center = TRUE, scale = FALSE)
  }

  # Run neighbourhood search
  fits <- robustsubsets:::ns(x, y, k, h, max.ns.iter, max.gd.iter, eps)

  # Run mio
  fits$mipgap <- array(dim = c(nk, nh))
  if (!any(is.null(k.mio)) & !any(is.null(h.mio))) {
    for (j in 1:nh.mio) {
      for (i in 1:nk.mio) {
        if (k.mio[i] == 0) next
        k.ind <- which(k == k.mio[i])
        h.ind <- which(h == h.mio[j])
        fit <- mio(x, y, k.mio[i], h.mio[j], fits$beta[, k.ind, h.ind], fits$eta[, k.ind, h.ind],
                   tau, params, warm.start)
        fits$beta[, k.ind, h.ind] <- fit$beta
        fits$eta[, k.ind, h.ind] <- fit$eta
        fits$objval[k.ind, h.ind] <- fit$objval
        fits$mipgap[k.ind, h.ind] <- fit$mipgap
      }
    }
  }

  # Refit the models on the original (unscaled) data
  fits.final <- list(beta = array(dim = c(p + 1, nk, nh)), weights = array(dim = c(n, nk, nh)),
                     objval = array(dim = c(nk, nh)), mipgap = fits$mipgap, k = k, h = h)
  colnames(fits.final$objval) <- colnames(fits.final$mipgap) <-
    dimnames(fits.final$beta)[[3]] <- dimnames(fits.final$weights)[[3]] <- h
  rownames(fits.final$objval) <- rownames(fits.final$mipgap) <-
    dimnames(fits.final$beta)[[2]] <- dimnames(fits.final$weights)[[2]] <- k
  for (j in 1:nh) {
    for (i in 1:nk) {
      if (any(fits$beta[, i, j] != 0)) {
        beta <- fits$beta[, i, j]
        eta <- fits$eta[, i, j]
        beta[abs(beta) < 1e-5] <- 0
        eta[abs(eta) < 1e-5] <- 0
        id.beta <- 1:p %in% which(beta != 0)
        id.eta <- 1:n %in% which(eta != 0)
        beta <- c(0, beta)
        beta[c(TRUE, id.beta)] <- stats::lsfit(x.o[!id.eta, id.beta], y.o[!id.eta])$coef
      } else {
        eta <- H(y.o, n - h[j])
        id.eta <- 1:n %in% which(eta != 0)
        beta <- matrix(0, p + 1)
        beta[1] <- mean(y.o[!id.eta])
      }
      eta <- y.o - cbind(1, x.o) %*% beta
      eta[!id.eta] <- 0
      objval <- f(cbind(1, x.o), y.o, beta, eta)
      fits.final$beta[, i, j] <- beta
      fits.final$weights[, i, j] <- as.integer(eta == 0)
      fits.final$objval[i, j] <- objval
    }
  }

  # Return fit
  class(fits.final) <- 'rss'
  return(fits.final)

}

#==================================================================================================#
# This function is responsible for mixed-integer optimisation. It contains the Gurobi specifications
# for solving the robust subset selection problem. The mio solver is capable of solving the problem
# to global optimality.
#==================================================================================================#

mio <- \(x, y, k, h, beta, eta, tau, params, warm.start) {

  # Preliminaries
  n <- nrow(x)
  p <- ncol(x)
  w <- cbind(x, diag(1, n))
  form <- ifelse(ncol(x) <= nrow(x) & h == nrow(x), 1, 2)

  # Set bounds for variables
  Mb <- tau * max(abs(beta))
  Mb <- ifelse(is.nan(Mb), Inf, Mb)
  Me <- tau * max(abs(eta))
  Me <- ifelse(is.nan(Me), Inf, Me)

  # Set the solver problem formulation
  model <- list()

  # Problem formulation 1 (low-dimensional)
  # -------------------------------------------------------------
  # Objective function: x ^ T Q x + c ^ T x + a
  # Problem variable: x = [beta, eta, s, z]
  # Constraints
  # C1: sum(s) >= p - k <==> ||beta||_0 <= k
  # C2:     sum(z) >= h <==> ||eta||_0 <= n - h
  # NOTE: Gurobi will transform SOS constraints to Big-M
  # constraints when coef bounds are finite
  # -------------------------------------------------------------
  if (form == 1) {
    model$vtype <- c(rep('C', p + n), rep('B', p + n)) # Problem variable x
    model$Q <- 0.5 * rbind(
      cbind(t(w) %*% w, matrix(0, p + n, p + n)),
      matrix(0, p + n, p + n + p + n)
    ) # matrix Q in objective function
    model$obj <- c(- t(w) %*% y, rep(0, p + n)) # Vector c in objective function
    model$objcon <- 0.5 * t(y) %*% y # Constant a in objective function
    model$A <- # LHS of constraints
      rbind(
        cbind(matrix(0, 1, p), matrix(0, 1, n), matrix(1, 1, p), matrix(0, 1, n)), # C1
        cbind(matrix(0, 1, p), matrix(0, 1, n), matrix(0, 1, p), matrix(1, 1, n)) # C2
      )
    model$sense <- c('>=', '>=') # Constraint types
    model$rhs <- c(p - k, h) # RHS of constraints
    model$lb <- c(rep(- Mb, p), rep(- Me, n), rep(0, p + n)) # Lower bound on variables
    model$ub <- c(rep(Mb, p), rep(Me, n), rep(1, p + n)) # Upper bound on variables
    model$sos <- vector('list', p + n)
    for (j in 1:(p + n)) { # SOS constraints
      model$sos[[j]]$type <- 1
      model$sos[[j]]$index <- c(j, p + n + j)
      model$sos[[j]]$weight <- c(1, 2)
    }
    if (warm.start) model$start <- c(beta, eta, beta == 0, eta == 0) # Warm start
  }

  # Problem formulation 2 (high-dimensional)
  # -------------------------------------------------------------
  # Objective function: x ^ T Q x + c ^ T x + a
  # Problem variable: x = [beta, eta, s, z, xi]
  # Constraints
  # C1:        sum(s) >= p - k <==> ||beta||_0 <= k
  # C2:            sum(z) >= h <==> ||eta||_0 <= n - h
  # C3: w [beta, eta] - xi = 0 <==> w [beta, eta] = xi
  # NOTE: Gurobi will transform SOS constraints to Big-M
  # constraints when coef bounds are finite
  # -------------------------------------------------------------
  if (form == 2) {
    Mx <- tau * max(abs(w %*% c(beta, eta)))
    Mx <- ifelse(is.nan(Mx), Inf, Mx)
    model$vtype <- c(rep('C', p + n), rep('B', p + n), rep('C', n)) # Problem variable x
    model$Q <- 0.5 * rbind(
      matrix(0, p + n + p + n, p + n + p + n + n),
      cbind(matrix(0, n, p + n + p + n), diag(n))
    ) # matrix Q in objective function
    model$obj <- c(- t(w) %*% y, rep(0, p + n + n)) # Vector c in objective function
    model$objcon <- 0.5 * t(y) %*% y # Constant a in objective function
    model$A <- # LHS of constraints
      rbind(
        cbind(matrix(0, 1, p), matrix(0, 1, n), matrix(1, 1, p), matrix(0, 1, n), matrix(0, 1, n)), # C1
        cbind(matrix(0, 1, p), matrix(0, 1, n), matrix(0, 1, p), matrix(1, 1, n), matrix(0, 1, n)), # C2
        cbind(w, matrix(0, n, p), matrix(0, n, n), - diag(n)) # C3
      )
    model$sense <- c('>=', '>=', rep('=', n)) # Constraint types
    model$rhs <- c(p - k, h, rep(0, n)) # RHS of constraints
    model$lb <- c(rep(- Mb, p), rep(- Me, n), rep(0, p + n), rep(- Mx, n)) # Lower bound on coefs
    model$ub <- c(rep(Mb, p), rep(Me, n), rep(1, p + n), rep(Mx, n)) # Upper bound on coefs
    model$sos <- vector('list', p + n)
    for (j in 1:(p + n)) { # SOS constraints
      model$sos[[j]]$type <- 1
      model$sos[[j]]$index <- c(j, p + n + j)
      model$sos[[j]]$weight <- c(1, 2)
    }
    if (warm.start) model$start <- c(beta, eta, beta == 0, eta == 0, w %*% c(beta, eta))  # Warm start
  }
  # Solve the model using Gurobi
  fit <- gurobi::gurobi(model, params)
  if (is.null(fit$x)) x <- rep(0, p + n) # In case Gurobi fails to solve the relaxation
  fit$beta <- fit$x[1:p]
  fit$eta <- fit$x[(p + 1):(p + n)]

  return(fit)

}

#==================================================================================================#
# Coefficient function
#==================================================================================================#

#' @title Coefficient function for rss object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Extracts coefficients for a given parameter pair \code{(k,h)}.
#'
#' @param object an object of class \code{rss}
#' @param k the number of predictors indexing the desired fit
#' @param h the number of observations indexing the desired fit
#' @param ... any other arguments
#'
#' @return An array of coefficients.
#'
#' @method coef rss
#'
#' @export
#'
#' @importFrom stats "coef"

coef.rss <- \(object, k = NULL, h = NULL, ...) {

  if (!is.null(k)) index1 <- which(object$k == k) else index1 <- seq_along(object$k)
  if (!is.null(h)) index2 <- which(object$h == h) else index2 <- seq_along(object$h)
  object$beta[, index1, index2]

}

#==================================================================================================#
# Predict function
#==================================================================================================#

#' @title Predict function for rss object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Generate predictions for new data using a given parameter pair \code{(k,h)}.
#'
#' @param object an object of class \code{rss}
#' @param x.new a matrix of new values for the predictors
#' @param k the number of predictors indexing the desired fit
#' @param h the number of observations indexing the desired fit
#' @param ... any other arguments
#'
#' @return An array of predictions.
#'
#' @method predict rss
#'
#' @export
#'
#' @importFrom stats "predict"

predict.rss <- \(object, x.new, k = NULL, h = NULL, ...) {

  x.new <- as.matrix(x.new)
  x.new <- cbind(1, x.new)
  beta <- coef.rss(object, k, h, ...)
  if (length(dim(beta)) < 3) x.new %*% beta else apply(beta, 2:3, \(beta) x.new %*% beta)

}

#==================================================================================================#
# Plot function
#==================================================================================================#

#' @title Plot function for rss object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Plot the coefficient profiles from robust subset selection.
#'
#' @param x an object of class \code{rss}
#' @param ... any other arguments
#'
#' @return A plot of the coefficient profiles.
#'
#' @method plot rss
#'
#' @export
#'
#' @importFrom graphics "plot"

plot.rss <- \(x, ...) {

  beta <- x$beta
  beta <- beta[- 1, , , drop = FALSE]
  beta[beta == 0] <- NA
  beta <- data.frame(beta = as.vector(beta),
                     predictor = as.factor(rep(seq_len(nrow(beta)), length(x$k) * length(x$h))),
                     k = as.factor(rep(rep(x$k, each = nrow(beta), length(x$h)))),
                     h = as.factor(rep(x$h, each = nrow(beta) * length(x$k))))
  beta <- stats::na.omit(beta)
  ggplot2::ggplot(beta, ggplot2::aes_string('k', 'beta', col = 'predictor')) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(. ~ h, labeller = ggplot2::label_both)

}
