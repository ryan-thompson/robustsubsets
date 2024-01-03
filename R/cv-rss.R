#==================================================================================================#
# This function performs (repeated) K-fold cross-validation for tuning (k,h) in a robust subset
# selection model.
#==================================================================================================#

#' @title Cross-validated robust subset selection
#'
#' @author Ryan Thompson
#'
#' @description  Fits a sequence of regression models using robust subset selection and then
#' cross-validates these models.
#'
#' @param x a predictor matrix
#' @param y a response vector
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param h a function that takes the sample size that returns the number of observations to
#' minimise sum of squares over; by default produces a sequence from 75 to 100 percent of sample
#' size (in increments of 5 percent); a function is used here to facilitate varying sample sizes in
#' cross-validation
#' @param mio one of 'min', 'all', or 'none' indicating whether to run the mixed-integer solver on
#' the \code{k} and \code{h} that minimise the cv error, all \code{k} and \code{h}, or none at all
#' @param nfold the number of folds to use in cross-validation
#' @param cv.loss an optional cross-validation loss-function to use; should accept a vector of
#' errors; by default trimmed mean square prediction error with 25\% trimming
#' @param cluster an optional cluster for running cross-validation in parallel; must be set up using
#' \code{parallel::makeCluster}
#' @param ... any other arguments
#'
#' @return An object of class \code{cv.rss}; a list with the following components:
#' \item{cv}{a matrix with the cross-validated values of \code{cv.loss}; rows correspond to \code{k}
#' and columns to \code{h}}
#' \item{k}{a vector containing the values of \code{k} used in the fit}
#' \item{h}{a vector containing the values of \code{h} used in the fit}
#' \item{k.min}{the \code{k} yielding the lowest cross-validated \code{cv.loss}}
#' \item{h.min}{the \code{h} yielding the lowest cross-validated \code{cv.loss}}
#' \item{fit}{the fit from running \code{rss()} on the full data}
#'
#' @example R/examples/example-cv-rss.R
#'
#' @export
#'

cv.rss <- \(x, y, k = 0:min(nrow(x) - 1, ncol(x), 20), h = \(n) round(seq(0.75, 1, 0.05) * n),
            mio = 'min', nfold = 10, cv.loss = tmspe, cluster = NULL, ...) {

  # Check data is valid
  if (!is.matrix(x)) x <- as.matrix(x)
  attributes(x)$dimnames <- NULL
  if (!is.matrix(y)) y <- as.matrix(y)
  attributes(y)$dimnames <- NULL

  # Preliminaries
  n <- nrow(x)
  nk <- length(k)
  nh <- length(h(n))

  # Loop over folds
  folds <- sample(rep_len(1:nfold, n))
  cvf <- \(fold) {
    error <- array(dim = c(n, nk, nh))
    fold.ind <- which(folds == fold)
    x.train <- x[- fold.ind, , drop = FALSE]
    y.train <- y[- fold.ind, , drop = FALSE]
    x.valid <- x[fold.ind, , drop = FALSE]
    y.valid <- y[fold.ind, , drop = FALSE]
    fold.n <- n - length(fold.ind)
    fold.h <- h(fold.n)
    fit.fold <- rss(x.train, y.train, k, fold.h, ...)
    apply(fit.fold$beta, c(2, 3), \(beta) y.valid - cbind(1, x.valid) %*% beta)
  }
  if (is.null(cluster)) {
    error <- lapply(1:nfold, cvf)
  } else {
    parallel::clusterCall(cluster, \() library(robustsubsets))
    error <- parallel::clusterApply(cluster, 1:nfold, cvf)
  }
  error <- abind::abind(error, along = 1)

  # Compute cv metric and best (k,h) pair
  cv <- apply(error, c(2, 3), cv.loss)
  best <- which(cv == min(cv), arr.ind = TRUE)
  k.min <- k[best[1, 1]]
  h.min <- h(n)[best[1, 2]]

  # Fit the models
  k.mio <- NULL
  h.mio <- NULL
  if (mio == 'min') {
    k.mio <- k.min
    h.mio <- h.min
  } else if (mio == 'all') {
    k.mio <- k
    h.mio <- h(n)
  }
  fit <- rss(x, y, k, h(n), k.mio, h.mio, ...)

  # Return result
  result <- list(cv = cv, k = k, h = h(n), k.min = k.min, h.min = h.min, fit = fit)
  class(result) <- 'cv.rss'
  return(result)

}

#==================================================================================================#
# Coefficient function
#==================================================================================================#

#' @title Coefficient function for cv.rss object
#'
#' @author Ryan Thompson
#'
#' @description Extracts coefficients for a given parameter pair \code{(k,h)}.
#'
#' @param object an object of class \code{cv.rss}
#' @param k the number of predictors indexing the desired fit; 'k.min' uses best \code{k} from
#' cross-validation
#' @param h the number of observations indexing the desired fit; 'h.min' uses best \code{h} from
#' cross-validation
#' @param ... any other arguments
#'
#' @return An array of coefficients.
#'
#' @method coef cv.rss
#'
#' @export
#'
#' @importFrom stats "coef"

coef.cv.rss <- \(object, k = 'k.min', h = 'h.min', ...) {

  if (!is.null(k)) if (k == 'k.min') k <- object$k.min
  if (!is.null(h)) if (h == 'h.min') h <- object$h.min
  coef.rss(object$fit, k = k, h = h, ...)

}

#==================================================================================================#
# Predict function
#==================================================================================================#

#' @title Predict function for cv.rss object
#'
#' @author Ryan Thompson
#'
#' @description Generate predictions given new data using a given parameter pair \code{(k,h)}.
#'
#' @param object an object of class \code{cv.rss}
#' @param x.new a matrix of new values for the predictors
#' @param k the number of predictors indexing the desired fit; 'k.min' uses best \code{k} from
#' cross-validation
#' @param h the number of observations indexing the desired fit; 'h.min' uses best \code{h} from
#' cross-validation
#' @param ... any other arguments
#'
#' @return An array of predictions.
#'
#' @method predict cv.rss
#'
#' @export
#'
#' @importFrom stats "predict"

predict.cv.rss <- \(object, x.new, k = 'k.min', h = 'h.min', ...) {

  if (!is.null(k)) if (k == 'k.min') k <- object$k.min
  if (!is.null(h)) if (h == 'h.min') h <- object$h.min
  predict.rss(object$fit, x.new, k = k, h = h, ...)

}

#==================================================================================================#
# Plot function
#==================================================================================================#

#' @title Plot function for cv.rss object
#'
#' @author Ryan Thompson
#'
#' @description Plot the cross-validation results from robust subset selection.
#'
#' @param x an object of class \code{cv.rss}
#' @param ... any other arguments
#'
#' @return A plot of the cross-validation results.
#'
#' @method plot cv.rss
#'
#' @export
#'
#' @importFrom graphics "plot"

plot.cv.rss <- \(x, ...) {

  cv <- data.frame(cv = as.vector(x$cv), k = as.factor(rep(x$k, length(x$h))),
                   h = as.factor(rep(x$h, each = length(x$k))))
  ggplot2::ggplot(cv, ggplot2::aes_string('k', 'cv', col = 'h')) +
    ggplot2::geom_point() +
    ggplot2::ylab('cv error') +
    ggplot2::geom_vline(xintercept = as.factor(x$k.min), linetype = 'dotted', alpha = 0.5) +
    ggplot2::geom_hline(yintercept = x$mean.cv[which(x$k == x$k.min), which(x$h == x$h.min)],
                        linetype = 'dotted', alpha = 0.5)

}
