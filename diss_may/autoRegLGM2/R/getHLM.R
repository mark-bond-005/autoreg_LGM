#' Getting a maximum-likelihood estimate of the HLM
#'
#' This convenience function gets the maximum-likelihood estimate of the HLM.
#' @param thet An N by M matrix of ability estimates
#' @keywords
#' @export
#' @examples
#' medEst()
#' Maximum-likelihod HLM

getHLM <- function(thet){
  thetInits <- cbind(c(thet[ , 1], thet[ , 2], thet[ , 3], thet[ , 4]),
                     c(rep(0, nLGM), rep(1, nLGM), rep(2, nLGM), rep(3, nLGM)),
                     c(rep(0, nLGM), thet[ , 1], thet[ , 2], thet[ , 3]),
                     rep(1:nLGM, 4))
  thetTest <- as.data.frame(thetInits); names(thetTest) <- c("theta", "time", "lag", "unit")
  init.Model <- lmer(theta ~ time + lag + (1 + time + lag | unit), data = thetTest)
  return(init.Model)
}
