#' Sampling trait estimates
#'
#' After the forward filtering is complete, backwards sampling begins.
#' @param mMat N by t matrix of posterior ability estimates. Column gives the time point.
#' @param cMat N by t matrix of posterior variance estimates. Column gives the time point.
#' @param AMat N by t vector of prior ability estimates. Column gives the time point.
#' @param rMat N by t matrix of prior variance estimates. Column gives the time point.
#' @param beta2 Person-level autoregressive trends
#' @keywords Sampling filter backwards
#' @export
#' @examples
#' backFilt()
# Returns theta estimates

backFilt <- function(mMat, cMat, AMat, rMat,
                     beta2){
  nLGM <- length(beta2)
  M <- dim(mMat)[2]
  thet <- matrix(0, nLGM, M)
  thet[ , M] <- rnorm(nLGM, mean = mMat[ , M], sd = sqrt(cMat[ , M]))
  for(j in (M-1):1){
    h <- mMat[ , j] + beta2*cMat[ , j]/rMat[ , j+1]*(thet[ , j + 1]-AMat[ , j + 1])
    b <- cMat[ , j] - beta2^2 *cMat[ , j]^2/rMat[ , j + 1]
    thet[ , j] <- rnorm(nLGM, mean = h, sd = sqrt(b))
  }
  return(thet)
}
