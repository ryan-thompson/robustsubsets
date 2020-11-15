#==================================================================================================#
# This function is intended for the end-user; it performs model fitting and cross-validation.
#==================================================================================================#

#' @title Robust subset selection
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Fits a sequence of robust subset selection models and cross-validates the prediction
#' error from these models.
#'
#' @param X a matrix of predictors
#' @param y a vector of the response
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param h a function that takes the sample size that returns the number of observations to
#' minimise sum of squares over; by default produces a sequence from 75 to 100 percent of sample
#' size (in increments of 5 percent); a function is used here to facilitate varying sample sizes in
#' cross-validation
#' @param int a logical indicating whether to include an intercept
#' @param mio one of 'min', 'all', or 'none' indicating whether to run the mixed-integer solver on
#' the \code{k} and \code{h} that minimise the cv error, all \code{k} and \code{h}, or none at all
#' @param ... any other arguments (see \code{rss.fit} and \code{rss.cv})
#'
#' @return An object of class \code{rss}; a list with the following components:
#' \item{cv}{the output from \code{rss.cv}; see documentation}
#' \item{fit}{the output from \code{rss.fit}; see documentation}
#'
#' @details The function is simply a wrapper that combines \code{rss.fit} and \code{rss.cv} to fit
#' a sequence of robust subset selection models (indexed by \code{k} and \code{h}) and perform
#' cross-validation for these models. Fitting is initially done using heuristics. Afterwards, the
#' mixed-integer solver can be called to search for globally optimal solutions using the parameter
#' \code{mio}. By default \code{mio='min'}, indicating that the mixed integer solver will be run
#' only on the tuning parameters that produced the minimum cross-validation error. It is possible
#' to run the mixed-integer solver on all tuning parameters using \code{mio='all'}, or not to call
#' the solver at all using \code{mio='none'}.\cr
#' See \code{rss.fit} and \code{rss.cv} for further options controlling the model fit and
#' cross-validation.
#'
#' @example R/examples/example_rss.R
#'
#' @export

rss <- function(X, y,
                k = 0:min(nrow(X) - int, ncol(X), 20),
                h = function(n) round(seq(0.75, 1, 0.05) * n),
                int = T, mio = 'min', ...) {

  # Run cross-validation
  cv <- rss.cv(X, y, k, h, int, ...)

  # Fit the models
  k.mio <- NULL
  h.mio <- NULL
  if (mio == 'min') {
    k.mio <- cv$k.min
    h.mio <- cv$h.min
  } else if (mio == 'all') {
    k.mio <- k
    h.mio <- h(nrow(X))
  }
  fit <- rss.fit(X, y, k, h(nrow(X)), int, k.mio, h.mio, ...)

  # Save results
  result <- list()
  result$cv <- cv
  result$fit <- fit

  class(result) <- 'rss'
  return(result)

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
#' @param k the number of predictors indexing the desired fit; 'k.min' uses best \code{k} from
#' cross-validation
#' @param h the number of observations indexing the desired fit; 'h.min' uses best \code{h} from
#' cross-validation
#' @param ... any other arguments
#'
#' @return An array of coefficients.
#'
#' @method coef rss
#'
#' @export
#'
#' @importFrom stats "coef"

coef.rss <- function(object, k = 'k.min', h = 'h.min', ...) {

  if (!is.null(k)) if (k == 'k.min') k <- object$cv$k.min
  if (!is.null(h)) if (h == 'h.min') h <- object$cv$h.min
  coef.rss.fit(object$fit, k = k, h = h, ...)

}

#==================================================================================================#
# Predict function
#==================================================================================================#

#' @title Predict function for rss object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Generate predictions given new data using a given parameter pair \code{(k,h)}.
#'
#' @param object an object of class \code{rss}
#' @param X.new a matrix of new values for the predictors
#' @param k the number of predictors indexing the desired fit; 'k.min' uses best \code{k} from
#' cross-validation
#' @param h the number of observations indexing the desired fit; 'h.min' uses best \code{h} from
#' cross-validation
#' @param ... any other arguments
#'
#' @return An array of predictions.
#'
#' @method predict rss
#'
#' @export
#'
#' @importFrom stats "predict"

predict.rss <- function(object, X.new, k = 'k.min', h = 'h.min', ...) {

  if (!is.null(k)) if (k == 'k.min') k <- object$cv$k.min
  if (!is.null(h)) if (h == 'h.min') h <- object$cv$h.min
  predict.rss.fit(object$fit, X.new, k = k, h = h, ...)

}

#==================================================================================================#
# Plot function
#==================================================================================================#

#' @title Plot function for rss object
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Plot the cross-validation results or coefficient profiles from robust subset
#' selection.
#'
#' @param x an object of class \code{rss}
#' @param type one of 'cv' or 'profile'
#' @param ... any other arguments
#'
#' @return A plot of the cross-validation results or coefficient profiles.
#'
#' @method plot rss
#'
#' @export
#'
#' @importFrom graphics "plot"

plot.rss <- function(x, type = 'cv', ...) {

  if (type == 'profile') {
    plot.rss.fit(x$fit)
  } else if (type == 'cv') {
    plot.rss.cv(x$cv)
  }

}
