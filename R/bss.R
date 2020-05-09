#==================================================================================================#
# This is a wrapper function for best subset selection problem via the rss functions (note: best
# subsets is a special case of robust subsets).
#==================================================================================================#

#' @title Best subset selection
#'
#' @author Ryan Thompson <ryan.thompson@monash.edu>
#'
#' @description Fits a sequence of best subset selection models. This function is just a wrapper for
#' the \code{rss} function. The function solves the robust subset selection problem with
#' \code{h}=\code{n}, using nonrobust measures of location and scale to standardise, as well as a
#' nonrobust measure of prediction error in cross-validation.
#'
#' @param X a matrix of predictors
#' @param y a vector of the response
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param int a logical indicating whether to include an intercept
#' @param mio a logical indicating whether to run the mixed-integer solver
#' @param ... any other arguments
#'
#' @return See documentation for the \code{rss} function.
#'
#' @export

bss <- function(X, y, k = (!int):min(nrow(X) - int, ncol(X), 20), int = T, mio = T, ...) {

  rss(X, y, k, nrow(X), int, mio, robust = F, cv.objective = mspe, ...)

}
