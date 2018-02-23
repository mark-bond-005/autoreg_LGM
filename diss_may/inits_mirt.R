
### Grab initial theta values

allButConstrain <- "
  F1 = 1-30
F2 = 31-60
F3 = 61-90
F4 = 91-120
MEAN = F2, F3, F4
COV = F1*F2*F3*F4
CONSTRAIN = "

allQuestions <- function(x){
  x1 <- paste("(", paste(x, x+30, x+60, x+90, sep = ", "), ", d)", sep="")
  x2 <- paste("(", paste(x, x+30, x+60, x+90, sep = ", "), ", a1, a2, a3, a4)", sep="")
  paste(x1, x2, sep=", ")
}

query <- allButConstrain
for(i in 1:30){
  query <- ifelse(i == 1, paste(query, allQuestions(i)), paste(query, allQuestions(i), sep=", "))
}
mod <- mirt.model(query)
init.IRT<- mirt(data=responses, mod, verbose=F, method="QMCEM")
#test_ <- ltm(longResp ~ z1, constraint = rbind(c(1, 1, 0), c(1, 2, 1)), IRT.param = F)
ak <- coef(init.IRT, simplify=T)$items[,1][1:30]
bk <- coef(init.IRT, simplify=T)$items[,5][1:30]
ak;bk
thet <- fscores(init.IRT, method="MAP")

### BEGIN CHAINS
nChain = 3
chain_ak <- array(data=NA, dim=c(length(ak),nChain))
chain_bk <- array(data=NA, dim=c(length(ak),nChain))
chain_thet <- array(data=NA, dim=c(dim(thet), nChain))
chain_ak[ , 1] <- ak
chain_bk[ , 1] <- bk
chain_thet[ , , 1] <- thet
chain_ak[ , 2:3] <- ak + rnorm(n= 2*length(ak), sd = 0.1)
chain_bk[ , 2:3] <- bk + rnorm(n= 2*length(bk), sd = 0.1)
chain_thet[ , , 2] <- thet + rnorm(n=dim(thet)[1]*dim(thet)[2], sd=0.1)
chain_thet[ , , 3] <- thet + rnorm(n=dim(thet)[1]*dim(thet)[2], sd=0.1)

# Use thetas to get initial estimates of HLM terms.
getHLM <- function(thet){
  thetInits <- cbind(c(thet[ , 1], thet[ , 2], thet[ , 3], thet[ , 4]), 
                  c(rep(0, nLGM), rep(1, nLGM), rep(2, nLGM), rep(3, nLGM)),
                  c(rep(0, nLGM), thet[ , 1], thet[ , 2], thet[ , 3]), 
                  rep(1:nLGM, 4))
  thetTest <- as.data.frame(thetInits); names(thetTest) <- c("theta", "time", "lag", "unit")
  init.Model <- lmer(theta ~ time + lag + (1 + time + lag | unit), data = thetTest)
  return(init.Model)
}
#1000 -0.826476589 -0.0703562280 0.62037894

init.Model<- apply(chain_thet, 3, getHLM)
#alpha <- coef(init.Model)[ , 1]; beta1 <- coef(init.Model)[ , 2]; beta2 <- coef(init.Model)[ , 3]

chain_varL1 <- sapply(init.Model, FUN=sigma)^2
chain_gammas <- sapply(init.Model, lme4::fixef)
#chain_betas <- lapply(init.Model, coef)

getTMat <- function(x){
  xDat <- as.data.frame(VarCorr(x))
  matrix(data = c(xDat[1,4], xDat[4,4], xDat[5,4], xDat[4,4], xDat[2,4],xDat[5,4],xDat[5,4],xDat[6,4], xDat[3,4]), nrow = 3)
} 

prechain_TMat <- lapply(init.Model, getTMat)
chain_TMat <- array(data = c(prechain_TMat[[1]], prechain_TMat[[2]], prechain_TMat[[3]]), dim = c(3,3,3))
chain_betas <- array(data=0, dim=c(nLGM, 3, 3))
chain_betas[,1,1] <- as.vector(coef(init.Model[[1]])[[1]])[,1]
chain_betas[,2,1] <- as.vector(coef(init.Model[[1]])[[1]])[,2]
chain_betas[,3,1] <- as.vector(coef(init.Model[[1]])[[1]])[,3]
chain_betas[,1,2] <- as.vector(coef(init.Model[[2]])[[1]])[,1]
chain_betas[,2,2] <- as.vector(coef(init.Model[[2]])[[1]])[,2]
chain_betas[,3,2] <- as.vector(coef(init.Model[[2]])[[1]])[,3]
chain_betas[,1,3] <- as.vector(coef(init.Model[[3]])[[1]])[,1]
chain_betas[,2,3] <- as.vector(coef(init.Model[[3]])[[1]])[,2]
chain_betas[,3,3] <- as.vector(coef(init.Model[[3]])[[1]])[,3]

chain_TMat[1, ,1] <- chain_TMat[1,,1]*1/sqrt(chain_TMat[1,1,1])
chain_TMat[,1,1] <- chain_TMat[,1,1]*1/chain_TMat[1,1,1]

chain_TMat[1,,2] <- chain_TMat[1,,2]*1/sqrt(chain_TMat[1,1,2])
chain_TMat[,1,2] <- chain_TMat[,1,2]*1/chain_TMat[1,1,2]

chain_TMat[1,,3] <- chain_TMat[1,,3]*1/sqrt(chain_TMat[1,1,3])
chain_TMat[,1,3] <- chain_TMat[,1,3]*1/chain_TMat[1,1,3]