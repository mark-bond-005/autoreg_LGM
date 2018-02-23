### We need data before we can estimate the model I am proposing.
# This file generates appropriate data, using the variables nLGM and phi from master_file.R

# The a and b parameters refer to the slope and difficulty of the items, respectively.
bparms <- parameters[seq(1, 59, 2), "b"]
aparms <- parameters[seq(1, 59, 2), "a"]
nItems <- length(bparms)
occVar <- 1/3

growParms <- rtmvnorm(nLGM, mean=c(0, 0.34, phi), 
                      sigma = rbind(c(1, 0, 0), c(0, 0.17, 0), c(0, 0, 0.273)),
                      lower= c(-Inf, -Inf, -1), upper = c(Inf, Inf, 1))

# Here I generate initial ability, growth, and autoregression for each simulee. 
theta_ <- growParms[ , 1]; beta1_ <- growParms[ , 2]; beta2_ <- growParms[ , 3]
m0 <- rnorm(n=nLGM, sd=sqrt(occVar))
#note m0 below
theta1 <- theta_ + 0*beta1_ + 0*beta2_ + rnorm(n = nLGM, sd = sqrt(1-occVar^2))
theta2 <- theta_ + 1*beta1_ + theta1*beta2_ + rnorm(n = nLGM, sd = sqrt(1-occVar^2))
theta3 <- theta_ + 2*beta1_ + theta2*beta2_ + rnorm(n = nLGM, sd = sqrt(1-occVar^2))
theta4 <- theta_ + 3*beta1_ + theta3*beta2_ + rnorm(n = nLGM, sd = sqrt(1-occVar^2))

# This function generates responses based on the model.
genResp <- function(nLGM, nItems, theta, timePt){
  responses <- matrix(0, nLGM, nItems)
  for(i in 1:nLGM){
    tmp <- rep(theta[i], nItems)
    prob <- exp(tmp*aparms - bparms)/(1+exp(tmp*aparms - bparms))
    resp <- rbinom(nItems, 1, prob)
    responses[i, ] <- resp
  }
  responses <- as.data.frame(responses)
  names(responses) <- c(paste("t", timePt, paste("q", 1:nItems, sep = ""), sep = ""))
  return(responses)
}

responses1 <- genResp(nLGM, nItems, theta1, "1"); responses1$ID <- 1:nLGM
responses2 <- genResp(nLGM, nItems, theta2, "2"); responses2$ID <- 1:nLGM
responses3 <- genResp(nLGM, nItems, theta3, "3"); responses3$ID <- 1:nLGM
responses4 <- genResp(nLGM, nItems, theta4, "4"); responses4$ID <- 1:nLGM

responses <- merge(responses1, responses2, by="ID") 
responses <- merge(responses, responses3, by="ID") 
responses <- merge(responses, responses4, by="ID")

responses <- responses[, c(-1)]
