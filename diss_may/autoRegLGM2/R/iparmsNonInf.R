#' Non-informative item sampling
#'
#' Without using any prior knowledge, this function gets estimates of item traits
#' @param thet an N by M matrix of ability estimates
#' @param zMat an N by I*M matrix of latent measures
#' @keywords Sampling item non-informative
#' @export
#' @examples
#' iparmsNonInf()
# Returns item discrimination and difficulty measures

iparmsNonInf <- function(thet, zMat){
  nTimes <- dim(thet)[2]
  nLGM <- dim(thet)[1]
  nItems <- dim(zMat)[2]/nTimes
  thetvec <- as.vector(thet)
  H <- cbind(thetvec, rep(1, nLGM*nTimes))
  covar <- solve(t(H)%*%H)
  preMean <- covar%*%t(H)
  ak <- vector(mode = "numeric"); bk <- vector(mode = "numeric"); #ak[1] = 1; bk[1] = 0;
  for (k in 1:nItems){
    out <- matrix(data = zMat[ , seq(from=k, to=k+nItems*(nTimes-1), by = nItems)], nrow = nTimes*nLGM)
    mean <- preMean%*%out
    x <- mvrnorm(mu=mean, Sigma=covar)
    ak[k] <- x[1]; bk[k] <- x[2]
  }
  ret <- list(); ret$ak <- ak; ret$bk <- -bk;
  return(ret);
}
