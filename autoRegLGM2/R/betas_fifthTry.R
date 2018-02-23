#' Sampling level-one HLM parameters
#'
#' This function gets person level initial abilities, person-level growth terms, and person-level autoregressive trends
#' @param thet an N by M matrix of ability estimates
#' @param TMat a 3 by 3 level-two variance-covariance matrix
#' @param varL1 Level-one variance
#' @param gammas Level-two means
#' @keywords Sampling HLM level one
#' @export
#' @examples
#' betas_fifthTry()
# Returns person-level HLM variables

betas_fifthTry <- function(thet, TMat, varL1, gammas){
  meanVarL1 <- .Call('autoRegLGM2_getMeanVar', PACKAGE = 'autoRegLGM2', thet, as.matrix(gammas, ncol=3), varL1, TMat)
  nLGM <- dim(meanVarL1$means)[1]
  alpha <- vector(mode="numeric"); beta1 <- vector(mode="numeric"); beta2 <- vector(mode="numeric")
  # using apply instead does not really speed this process up.
  for(indivI in 1:nLGM)
  {
    mean <- meanVarL1$means[indivI, ]
    covar <- matrix(meanVarL1$D[indivI, ], ncol = length(gammas))
    alphaBeta <- try(mvrnorm(n=1, mu=mean, Sigma=covar), silent=T)
    if (class(alphaBeta) == "try-error"){
      newSigma <- matrix(c(covar[1,1], 0,0, 0, covar[2,2], 0, 0,0, covar[3,3]), ncol=length(gammas))
      alphaBeta <- mvrnorm(n=1, mu=mean, Sigma=newSigma)
    }
    alpha[indivI] <- alphaBeta[1]
    beta1[indivI] <- alphaBeta[2]
    #beta2[indivI] <- alphaBeta[3]
    beta2[indivI] <- rtruncnorm(n=1, a=-1, b=1, mean=mean[3], sd=sqrt(covar[3,3]))
  }
  betas <- cbind(alpha, beta1, beta2)
  return(betas)
}
