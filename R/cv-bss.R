#==================================================================================================#
# This function is a wrapper function for cross-validated best subset selection via the rss.cv
# function
#==================================================================================================#

#' @title Cross-validated best subset selection
#'
#' @author Ryan Thompson
#'
#' @description Fits a sequence of regression models using best subset selection and then
#' cross-validates these models. This function is just a wrapper for the \code{cv.rss} function.
#' The function solves the robust subset selection problem with \code{h}=\code{n}, using nonrobust
#' measures of location and scale to standardise, and a nonrobust measure of prediction error in
#' cross-validation.
#'
#' @param x a predictor matrix
#' @param y a response vector
#' @param k the number of predictors to minimise sum of squares over; by default a sequence from 0
#' to 20
#' @param mio one of 'min', 'all', or 'none' indicating whether to run the mixed-integer solver on
#' the \code{k} that minimises the cv error, all \code{k}, or none at all
#' @param nfold the number of folds to use in cross-validation
#' @param cv.loss an optional cross-validation loss-function to use; should accept a vector of
#' errors; by default mean square prediction error
#' @param ... any other arguments
#'
#' @return See documentation for the \code{cv.rss} function.
#'
#' @example R/examples/example-cv-bss.R
#'
#' @export

cv.bss <- \(x, y, k = 0:min(nrow(x) - 1, ncol(x), 20), mio = 'min', nfold = 10, cv.loss = mspe, ...) {

  cv.rss(x, y, k, \(n) n, mio, nfold, cv.loss, robust = FALSE,  ...)

}
