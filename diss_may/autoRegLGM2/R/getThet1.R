#' Obtaining Initial trait estimates
#'
#' Before an FFBS is run, we need an initial estimate of the trait.
#' This function uses values from getZs(), the HLM, and item characteristics to get this initial estimate.
#' It returns a list of precisions and means for each person.
#' @param zMat Latent Normally distributed outcome variables
#' @param alpha initial ability estimates from the HLM
#' @param nItems Number of items
#' @param ak A numeric vector of item-level discrimination values
#' @param bk A numeric vector of item-level difficulties
#' @param varL1
#' @keywords Initial initializing
#' @export
#' @examples
#' getThet1()
# Returns Initial theta values
getThet1 <- function(zMat, alpha, ak, bk, varL1) {
  v = 1/sum(ak*ak)
  nItems <- length(ak)
  zMat1 <- zMat[ , 1:nItems]
  for(k in 1:nItems){
    zMat1[ , k] <- ak[k]*(zMat[ , k] + bk[k])
  }
  thetVec1 <- apply(zMat1, MARGIN = 1, sum)*v
  thet1Precis <- 1/(1/v + 1/varL1)
  thet1Mean <- thetVec1/v + alpha/varL1
  ret <- list()
  ret$precis <- thet1Precis
  ret$mean <- thet1Mean
  return(ret)
}
