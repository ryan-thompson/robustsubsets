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
#'  parallel. To achieve good run time, only uses the heuristics (by default).
#'
#' @param X a matrix of predictors
#' @param y a vector of the response
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param h the number of observations to minimise sum of squares over; by default a sequence from
#' 75 to 100 percent of sample size (in increments of 5 percent)
#' @param int a logical indicating whether to include an intercept
#' @param n.fold the number of folds to use in cross-validation
#' @param n.cv the number of times to repeat cross-validation; the results are averaged
#' @param n.cores the number of cores to use in cross-validation; by default all cores are used
#' @param cv.objective the cross-validation objective function; by default trimmed mean square
#' prediction error with 25 percent trimming
#' @param ... any other arguments
#'
#' @return An object of class \code{rss.cv}; a list with the following components:
#' \item{mean.cv}{a matrix with the cross-validated values of \code{cv.objective}; each row
#' corresponds to a value of \code{k} and each column to a value of \code{h}}
#' \item{min.k}{the \code{k} yielding the lowest cross-validated \code{cv.objective}}
#' \item{min.h}{the \code{h} yielding the lowest cross-validated \code{cv.objective}}
#' \item{k}{the value of \code{k} that was passed in}
#' \item{h}{the value of \code{h} that was passed in}
#'
#' @export
#'
#' @importFrom foreach "%dopar%"
#' @importFrom foreach "%:%"

rss.cv <- function(X, y, k = (!int):min(nrow(X) - int, ncol(X), 20),
                   h = floor(seq(0.75, 1, 0.05) * nrow(X)), int = T, n.fold = 10, n.cv = 1,
                   n.cores = parallel::detectCores(), cv.objective = tmspe, ...) {

  # Preliminaries
  X <- as.matrix(X)
  y <- as.matrix(y)
  n <- nrow(X)
  n.k <- length(k)
  n.h <- length(h)

  # Compute errors via parallel n.cv * n.fold cross-validation
  folds <- matrix(rep_len(1:n.fold, n * n.cv), n, n.cv)
  folds <- apply(folds, 2, sample)
  cb.fold <- function(...) abind::abind(..., along = 1)
  cb.cv <- function(...) abind::abind(..., along = 4)
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)
  errors <- foreach::foreach(cv = 1:n.cv, .combine = 'cb.cv', .multicombine = T, .inorder = F) %:%
    foreach::foreach(fold = 1:n.fold, .combine = 'cb.fold', .multicombine = T, .inorder = F) %dopar% {
      fold.ind <- which(folds[, cv] == fold)
      X.train <- X[- fold.ind, ]
      y.train <- y[- fold.ind]
      X.valid <- X[fold.ind, , drop = F]
      y.valid <- y[fold.ind]
      fold.n <- n - length(fold.ind)
      fold.h <- round(h / n * fold.n)
      fit <- rss.fit(X.train, y.train, k, fold.h, int, ...)
      if (int) {
        error <- apply(fit$beta, c(2, 3), function(x) y.valid - cbind(1, X.valid) %*% x)
      } else {
        error <- apply(fit$beta, c(2, 3), function(x) y.valid - X.valid %*% x)
      }
      dim(error) <- c(n - fold.n, n.k, n.h)
      error
    }
  parallel::stopCluster(cl)

  # Compute cv metric, mean of cv metric, and best (k,h) pair
  cv.est <- apply(errors, 2:length(dim(errors)), function(x) cv.objective(x))
  mean.cv <- apply(cv.est, 1:2, function(x) mean(x))
  best <- which(mean.cv == min(mean.cv), arr.ind = T)

  # Save results
  result <- list()
  result$mean.cv <- mean.cv
  result$min.k <- k[best[1, 1]]
  result$min.h <- h[best[1, 2]]
  result$k <- k
  result$h <- h

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
    ggplot2::geom_vline(xintercept = as.factor(x$min.k), linetype = 'dotted', alpha = 0.5) +
    ggplot2::geom_hline(yintercept = x$mean.cv[which(x$k == x$min.k), which(x$h == x$min.h)],
                        linetype = 'dotted', alpha = 0.5)

}
