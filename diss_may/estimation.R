### Generates summary statistics for the Gibbs sampler.

start_ <- burnIn
N_ <- nIter

M <- 3

# The convenience functions below handle most of the summarizing I'll need to do. 
rHat <- function(chain, start=start_, N=N_){
  M <- dim(chain)[2]
  chain <- chain[seq(from=start, to=N),]
  
  thetM <- apply(chain, 2, mean)
  thetVar <- apply(chain, 2, var)
  thetMean <- mean(chain)
  
  B <- N/(M-1)*sum((thetM-thetMean)^2)
  W <- 1/M*sum(thetVar)
  
  vHat <- (N-1)*W/N + (M+1)*B/(M*N)
  return(vHat/W)
}

getMeans <- function(chain, start=start_, N=N_){
  chain <- chain[,start:N,]
  apply(chain, 1, mean)
}

getMedians <- function(chain, start=start_, N=N_){
  chain <- chain[,start:N,]
  apply(chain, 1, median)
}


relBias <- function(est, parm, start=start_, N=N_){
  (est[start:N, ]-parm)/est[start:N, ]
}

MSE <- function(est, parm, start=start_, N=N_){
 mean((est[start:N,]-parm)^2) 
}

TVec <- as.vector(chain_TMatDraw)
TMat2 <- array(chain_TMatDraw, dim=c(9,N_,M))

thet_Rhat <- apply(chain_thetDraw, 1, rHat)
print(summary(thet_Rhat))
beta_Rhat<- apply(chain_betasDraw, 1, rHat)
#print(summary(beta_Rhat))
TMat_Rhat <- apply(TMat2, 1, rHat)
print(TMat_Rhat)
a_Rhat <- apply(chain_aDraw, 1, rHat)
b_Rhat <- apply(chain_bDraw, 1, rHat)
varL1_Rhat <- rHat(chain_varL1Draw)
print(varL1_Rhat)
gamma_Rhat <- apply(chain_gammaDraw, 1, rHat)
print(gamma_Rhat)

thetTrue <- c(theta1, theta2, theta3, theta4)
TMatTrue <- c(1, 0, 0, 0, 0.17, 0, 0, 0, 0.273)
betaTrue <- c(theta_, beta1_, beta2_)
varL1True <- 1/3
aTrue <- aparms
bTrue <- bparms

medThet <- getMedians(chain_thetDraw)
medBeta <- getMedians(chain_betasDraw)
medTMat <- getMedians(TMat2)
medA <- getMedians(chain_aDraw)
medB <- getMedians(chain_bDraw)
medvarL1 <- median(chain_varL1Draw[start_:N_, ])
medGamma <- getMedians(chain_gammaDraw)
