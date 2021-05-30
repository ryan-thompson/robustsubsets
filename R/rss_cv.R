#==================================================================================================#
# This function performs (repeated) K-fold cross-validation for tuning (k,h) in a robust subset
# selection model.
#==================================================================================================#

globalVariables('cv') # Used in parallel for loop, to pass build checks
globalVariables('fold') # Used in parallel for loop, to pass build checks

#' @title Cross-validation for robust subset selection
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Does (repeated) \code{K}-fold cross-validation for robust subset selection in
#'  parallel. In the interest of speed, uses heuristics without the mixed-integer solver.
#'
#' @param x a matrix of predictors
#' @param y a vector of the response
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param h a function that takes the sample size that returns the number of observations to
#' minimise sum of squares over; by default produces a sequence from 75 to 100 percent of sample
#' size (in increments of 5 percent); a function is used here to facilitate varying sample sizes in
#' cross-validation
#' @param int a logical indicating whether to include an intercept
#' @param nfold the number of folds to use in cross-validation
#' @param ncv the number of times to repeat cross-validation; the results are averaged
#' @param n.core the number of cores to use in cross-validation; by default \code{nfold} *
#' \code{ncv} cores are used (if available)
#' @param cv.objective the cross-validation objective function; by default trimmed mean square
#' prediction error with 25 percent trimming
#' @param ... any other arguments
#'
#' @return An object of class \code{rss.cv}; a list with the following components:
#' \item{mean.cv}{a matrix with the cross-validated values of \code{cv.objective}; each row
#' corresponds to a value of \code{k} and each column to a value of \code{h}}
#' \item{k.min}{the \code{k} yielding the lowest cross-validated \code{cv.objective}}
#' \item{h.min}{the \code{h} yielding the lowest cross-validated \code{cv.objective}}
#' \item{k}{the value of \code{k} that was passed in}
#' \item{h}{the value of \code{h} that was passed in}
#'
#' @references Thompson, R. (2021). 'Robust subset selection'. arXiv:
#' \href{https://arxiv.org/abs/2005.08217}{2005.08217}.
#'
#' @example R/examples/example_rss_cv.R
#'
#' @export
#'
#' @importFrom foreach "%dopar%"
#' @importFrom foreach "%:%"

rss.cv <- function(x, y,
                   k = 0:min(nrow(x) - int, ncol(x), 20),
                   h = function(n) round(seq(0.75, 1, 0.05) * n),
                   int = TRUE, nfold = 10, ncv = 1, n.core = nfold * ncv,
                   cv.objective = tmspe, ...) {

  # Preliminaries
  x <- as.matrix(x)
  y <- as.matrix(y)
  n <- nrow(x)
  n.k <- length(k)
  n.core <- min(parallel::detectCores(), n.core)

  # Compute errors via parallel ncv * nfold cross-validation
  folds <- matrix(rep_len(1:nfold, n), n, ncv)
  folds <- apply(folds, 2, sample)
  cb.fold <- function(...) abind::abind(..., along = 1)
  cb.cv <- function(...) abind::abind(..., along = 4)
  cl <- parallel::makeCluster(n.core)
  doParallel::registerDoParallel(cl)
  errors <- foreach::foreach(cv = 1:ncv, .combine = 'cb.cv', .multicombine = TRUE,
                             .inorder = FALSE) %:%
    foreach::foreach(fold = 1:nfold, .combine = 'cb.fold', .multicombine = TRUE,
                     .inorder = FALSE) %dopar% {
      fold.ind <- which(folds[, cv] == fold)
      x.train <- x[- fold.ind, ]
      y.train <- y[- fold.ind]
      x.valid <- x[fold.ind, , drop = FALSE]
      y.valid <- y[fold.ind]
      fold.n <- n - length(fold.ind)
      fold.h <- h(fold.n)
      n.h <- length(fold.h)
      fit <- rss.fit(x.train, y.train, k, fold.h, int, ...)
      if (int) {
        error <- apply(fit$beta, c(2, 3), function(x) y.valid - cbind(1, x.valid) %*% x)
      } else {
        error <- apply(fit$beta, c(2, 3), function(x) y.valid - x.valid %*% x)
      }
      dim(error) <- c(n - fold.n, n.k, n.h)
      error
    }
  parallel::stopCluster(cl)

  # Compute cv metric, mean of cv metric, and best (k,h) pair
  cv.est <- apply(errors, 2:length(dim(errors)), cv.objective)
  mean.cv <- apply(cv.est, 1:2, mean)
  best <- which(mean.cv == min(mean.cv), arr.ind = TRUE)

  # Save results
  result <- list()
  result$mean.cv <- mean.cv
  result$k.min <- k[best[1, 1]]
  result$h.min <- h(n)[best[1, 2]]
  result$k <- k
  result$h <- h(n)

  class(result) <- 'rss.cv'
  return(result)

}

#==================================================================================================#
# Plot function
#==================================================================================================#

#' @title Plot function for rss.cv object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Plot the cross-validation results from robust subset selection.
#'
#' @param x an object of class \code{rss.cv}
#' @param ... any other arguments
#'
#' @return A plot of the cross-validation results.
#'
#' @method plot rss.cv
#'
#' @export
#'
#' @importFrom graphics "plot"

plot.rss.cv <- function(x, ...) {

  cv <- data.frame(cv = as.vector(x$mean.cv),
                   k = as.factor(rep(x$k, length(x$h))),
                   h = as.factor(rep(x$h, each = length(x$k))))
  ggplot2::ggplot(cv, ggplot2::aes_string('k', 'cv', col = 'h')) +
    ggplot2::geom_point() +
    ggplot2::ylab('cv error') +
    ggplot2::scale_y_log10() +
    ggplot2::geom_vline(xintercept = as.factor(x$k.min), linetype = 'dotted', alpha = 0.5) +
    ggplot2::geom_hline(yintercept = x$mean.cv[which(x$k == x$k.min), which(x$h == x$h.min)],
                        linetype = 'dotted', alpha = 0.5)

}
