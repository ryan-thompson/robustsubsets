#==================================================================================================#
# This function is intended for the end-user; it performs model fitting and cross-validation.
#==================================================================================================#

#' @title Robust subset selection
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Fits a sequence of robust subset selection model and cross-validates the prediction
#' error from these models.
#'
#' @param X a matrix of predictors
#' @param y a vector of the response
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param h the number of observations to minimise sum of squares over; by default a sequence from
#' 75 to 100 percent of sample size (in increments of 5 percent)
#' @param int a logical indicating whether to include an intercept
#' @param mio a logical indicating whether to run the mixed-integer solver
#' @param ... any other arguments (see \code{rss.fit} and \code{rss.cv})
#'
#' @return An object of class \code{rss}; a list with the following components:
#' \item{cv}{the output from \code{rss.cv}; see documentation}
#' \item{fit}{the output from \code{rss.fit}; see documentation}
#'
#' @details This function fits a sequence of models and cross-validates the prediction
#' error associated with these models. In the interest of speed, these steps are carried out using
#' heuristic optimisation methods. The parameters that produce the lowest cv error are run through
#' the mixed-integer solver which (given sufficient time) will converge to a global minimum. \cr
#' See \code{rss.fit} and \code{rss.cv} for further options controlling the model fit and
#' cross-validation.
#'
#' @example R/examples/example_rss.R
#'
#' @export

rss <- function(X, y,
                k = (!int):min(nrow(X) - int, ncol(X), 20), h = floor(seq(0.75, 1, 0.05) * nrow(X)),
                int = T, mio = T, ...) {

  # Run cross-validation
  cv <- rss.cv(X, y, k, h, int, ...)

  # Fit the models
  k.mio <- ifelse(mio, cv$min.k, NA)
  h.mio <- ifelse(mio, cv$min.h, NA)
  fit <- rss.fit(X, y, k, h, int, k.mio, h.mio, ...)

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
#' @param k the number of predictors indexing the desired fit; 'min.k' uses best \code{k} from
#' cross-validation
#' @param h the number of observations indexing the desired fit; 'min.h' uses best \code{h} from
#' cross-validation
#' @param ... any other arguments
#'
#' @return A vector of coefficients.
#'
#' @method coef rss
#'
#' @export
#'
#' @importFrom stats "coef"

coef.rss <- function(object, k = 'min.k', h = 'min.h',...) {

  if (k == 'min.k') k <- object$cv$min.k
  if (h == 'min.h') h <- object$cv$min.h
  coef(object$fit, k = k, h = h, ...)

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
#' @param k the number of predictors indexing the desired fit; 'min.k' uses best \code{k} from
#' cross-validation
#' @param h the number of observations indexing the desired fit; 'min.h' uses best \code{h} from
#' cross-validation
#' @param ... any other arguments
#'
#' @return A vector of predictions.
#'
#' @method predict rss
#'
#' @export
#'
#' @importFrom stats "predict"

predict.rss <- function(object, X.new, k = 'min.k', h = 'min.h', ...) {

  X.new <- as.matrix(X.new)
  if (object$fit$int) {
    cbind(1, X.new) %*% coef.rss(object, k, h, ...)
  } else {
    X.new %*% coef.rss(object, k, h, ...)
  }

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
    plot(x$fit)
  } else if (type == 'cv') {
    plot(x$cv)
  }

}
