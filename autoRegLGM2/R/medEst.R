#' Getting the median of a chain
#'
#' This convenience function gets the median of a markov chain after tossing out the burn-in values
#' @param parm the chain
#' @param low The lowest value of the chain you want
#' @param high The highest value of the chain you want
#' @keywords
#' @export
#' @examples
#' @useDynLib autoregLGM
#' medEst()
# Returns the median estimate from a markov chain
medEst <- function(parm, low, high){
  apply(parm[ , low:high], 1, median)
}
