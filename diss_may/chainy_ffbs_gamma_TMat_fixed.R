# Each iteration takes about 4 seconds
# The forward filter takes about 0.9 seconds
# HLM level 1 takes 0.3 seconds. 

# Setup draw matrices and M matrix
  M <- nIter
  chain_thetDraw <- array(NA, dim=c(nLGM*4, M, nChain)) 
  chain_aDraw <- array(NA, dim=c(nItems, M, nChain)) 
  chain_bDraw <- array(NA, dim=c(nItems, M, nChain)) 
  chain_betasDraw <- array(NA, dim=c(nLGM*3, M, nChain)) 
  chain_varL1Draw <- array(NA, dim=c(M, nChain))
  chain_TMatDraw <- array(NA, dim=c(3,3 , M, nChain))
  chain_gammaDraw <- array(NA, dim=c(3, M, nChain))
  
  gamSum2 <- list();   gamSum <- list()
  zMat <- array(NA, dim = c(nLGM, nItems*4, nChain))

# Most of these 
for (i in 1:nIter){
  # Get Zs
  #Sys.time()
  zMat[ ,  ,1] <- getZs(chain_thet[ , , 1], responses, chain_ak[ , 1], chain_bk[ , 1])
  zMat[ ,  ,2] <- getZs(chain_thet[ , , 2], responses, chain_ak[ , 2], chain_bk[ , 2])
  zMat[ ,  ,3] <- getZs(chain_thet[ , , 3], responses, chain_ak[ , 3], chain_bk[ , 3])
  
  # Use Zs to estimate the first measure of theta
  thet1_1 <- getThet1(zMat[,,1], chain_betas[,1,1], chain_ak[,1], chain_bk[,1], chain_varL1[1])
  thet1_2 <- getThet1(zMat[,,2], chain_betas[,1,2], chain_ak[,2], chain_bk[,2], chain_varL1[2])
  thet1_3 <- getThet1(zMat[,,3], chain_betas[,1,3], chain_ak[,3], chain_bk[,3], chain_varL1[3])
  
  # Now FFBS
  # Filter forward to start
  tmp_1 <- .Call('autoRegLGM2_forFilt2', PACKAGE = 'autoRegLGM2', 
                 as.matrix(chain_thet[,,1]), as.matrix(zMat[,,1]), chain_betas[,1,1], chain_betas[,2,1], rep(chain_gammas[3,1],nLGM),
                 chain_ak[,1], chain_bk[,1], thet1_1$mean*thet1_1$precis, chain_varL1[1], thet1_1$precis)
  mMat_1 <- tmp_1$mMat; cMat_1 <- tmp_1$cMat;
  AMat_1 <- tmp_1$AMat; rMat_1 <- tmp_1$rMat
  
  tmp_2 <- .Call('autoRegLGM2_forFilt2', PACKAGE = 'autoRegLGM2', 
                 as.matrix(chain_thet[,,2]), as.matrix(zMat[,,2]), chain_betas[,1,2], chain_betas[,2,2], rep(chain_gammas[3,2],nLGM),  
                 chain_ak[,2], chain_bk[,2],  thet1_2$mean*thet1_2$precis, chain_varL1[2], thet1_2$precis)
  mMat_2 <- tmp_2$mMat; cMat_2 <- tmp_2$cMat;
  AMat_2 <- tmp_2$AMat; rMat_2 <- tmp_2$rMat
  
  tmp_3 <- .Call('autoRegLGM2_forFilt2', PACKAGE = 'autoRegLGM2', 
                 as.matrix(chain_thet[,,3]), as.matrix(zMat[,,3]), 
                 chain_betas[,1,3], chain_betas[,2,3], rep(chain_gammas[3,3],nLGM),
                 chain_ak[,3], chain_bk[,3],
                 thet1_3$mean*thet1_3$precis, chain_varL1[3], thet1_3$precis)
  mMat_3 <- tmp_3$mMat; cMat_3 <- tmp_3$cMat;
  AMat_3 <- tmp_3$AMat; rMat_3 <- tmp_3$rMat
  # Forward filter complete
  # Now backwards filter
  chain_thet[,,1] <- backFilt(mMat_1, cMat_1, AMat_1, rMat_1, 
                   rep(chain_gammas[3,1],nLGM))
  
  chain_thet[,,2] <- backFilt(mMat_2, cMat_2, AMat_2, rMat_2, 
                   rep(chain_gammas[3,2],nLGM))
  
  chain_thet[,,3] <- backFilt(mMat_3, cMat_3, AMat_3, rMat_3, 
                   rep(chain_gammas[3,3],nLGM))
  # Get item parameters
  iparms_1 <- iparmsNonInf(chain_thet[,,1], zMat[,,1])
  chain_ak[,1] <- iparms_1$ak; chain_bk[,1] <- iparms_1$bk
  
  iparms_2 <- iparmsNonInf(chain_thet[,,2], zMat[,,2])
  chain_ak[,2] <- iparms_2$ak; chain_bk[,2] <- iparms_2$bk
  
  iparms_3 <- iparmsNonInf(chain_thet[,,3], zMat[,,3])
  chain_ak[,3] <- iparms_3$ak; chain_bk[,3] <- iparms_3$bk
  
  # Now get the betas
  # Often run into non-positive definite issues - if that's the case, just keep the old betas. Check chains later.
  betaTry1 <- try(betas_fifthTry(chain_thet[,,1], as.matrix(chain_TMat[,,1]), chain_varL1[1], chain_gammas[,1]))
  betaTry2 <- try(betas_fifthTry(chain_thet[,,2], as.matrix(chain_TMat[,,2]), chain_varL1[2], chain_gammas[,2]))
  betaTry3 <- try(betas_fifthTry(chain_thet[,,3], as.matrix(chain_TMat[,,3]), chain_varL1[3], chain_gammas[,3]))

  chain_betas[,,1] <- ifelse(matrix(class(betaTry1) == "try-error", nrow = nLGM, ncol = 3), 
                             chain_betas[,,1], betaTry1)
  chain_betas[,,2] <- ifelse(matrix(class(betaTry2) == "try-error", nrow = nLGM, ncol = 3), 
                             chain_betas[,,2], betaTry2)
  chain_betas[,,3] <- ifelse(matrix(class(betaTry3) == "try-error", nrow = nLGM, ncol = 3), 
                             chain_betas[,,3], betaTry3)
  
  # Grab gammas
  for(i2 in 1:nLGM){
    gamSum[[i2]] <- solve(as.matrix(chain_TMat[,,1]))%*%chain_betas[i2,,1] 
  }
  means <- Reduce('+', gamSum)
  precis <- solve(nLGM*solve(as.matrix(chain_TMat[,,1])))
  gam_tmp <- mvrnormArma(n = 1, mu = (precis%*%means)[2:3], sigma = precis[2:3, 2:3])
  chain_gammas[,1] <- c(0, gam_tmp)
  
  for(i2 in 1:nLGM){
    gamSum[[i2]] <- solve(as.matrix(chain_TMat[,,2]))%*%chain_betas[i2,,2] 
  }
  means <- Reduce('+', gamSum)
  precis <- solve(nLGM*solve(as.matrix(chain_TMat[,,2])))
  gam_tmp <- mvrnormArma(n = 1, mu = (precis%*%means)[2:3], sigma = precis[2:3, 2:3])
  chain_gammas[,2] <- c(0, gam_tmp)
  
  for(i2 in 1:nLGM){
    gamSum[[i2]] <- solve(as.matrix(chain_TMat[,,3]))%*%chain_betas[i2,,3] 
  }
  means <- Reduce('+', gamSum)
  precis <- solve(nLGM*solve(as.matrix(chain_TMat[,,3])))
  gam_tmp <- mvrnormArma(n = 1, mu = (precis%*%means)[2:3], sigma = precis[2:3, 2:3])
  chain_gammas[,3] <- c(0, gam_tmp)
  
  
  # get varL1
  s1 <- rbind(chain_thet[,1,1]-chain_betas[,1,1])%*%cbind((chain_thet[,1,1]-chain_betas[,1,1]))
  s2 <- rbind(chain_thet[,2,1]-chain_betas[,1,1]-chain_betas[,3,1]*chain_thet[,1,1] - chain_betas[,2,1]*1)%*%cbind(chain_thet[,2,1]-chain_betas[,1,1]-chain_betas[,3,1]*chain_thet[,1,1] - chain_betas[,2,1]*1)
  s3 <- rbind(chain_thet[,3,1]-chain_betas[,1,1]-chain_betas[,3,1]*chain_thet[,2,1] - chain_betas[,2,1]*2)%*%cbind(chain_thet[,3,1]-chain_betas[,1,1]-chain_betas[,3,1]*chain_thet[,2,1] - chain_betas[,2,1]*2)
  s4 <- rbind(chain_thet[,4,1]-chain_betas[,1,1]-chain_betas[,3,1]*chain_thet[,3,1] - chain_betas[,2,1]*3)%*%cbind(chain_thet[,4,1]-chain_betas[,1,1]-chain_betas[,3,1]*chain_thet[,3,1] - chain_betas[,2,1]*3)
  s_sq <- (s1+s2+s3+s4)/(4*nLGM)
  #sigL1 <- 1/rgamma(1, 4*nLGM/2, s_sq/2)
  chain_varL1[1] <- rinvchisq(n = 1, df = 4*nLGM, scale = s_sq)
  
  s1 <- rbind(chain_thet[,1,2]-chain_betas[,1,2])%*%cbind((chain_thet[,1,2]-chain_betas[,1,2]))
  s2 <- rbind(chain_thet[,2,2]-chain_betas[,1,2]-chain_betas[,3,2]*chain_thet[,1,2] - chain_betas[,2,2]*1)%*%cbind(chain_thet[,2,2]-chain_betas[,1,2]-chain_betas[,3,2]*chain_thet[,1,2] - chain_betas[,2,2]*1)
  s3 <- rbind(chain_thet[,3,2]-chain_betas[,1,2]-chain_betas[,3,2]*chain_thet[,2,2] - chain_betas[,2,2]*2)%*%cbind(chain_thet[,3,2]-chain_betas[,1,2]-chain_betas[,3,2]*chain_thet[,2,2] - chain_betas[,2,2]*2)
  s4 <- rbind(chain_thet[,4,2]-chain_betas[,1,2]-chain_betas[,3,2]*chain_thet[,3,2] - chain_betas[,2,2]*3)%*%cbind(chain_thet[,4,2]-chain_betas[,1,2]-chain_betas[,3,2]*chain_thet[,3,2] - chain_betas[,2,2]*3)
  s_sq <- (s1+s2+s3+s4)/(4*nLGM)
  #sigL1 <- 1/rgamma(1, 4*nLGM/2, s_sq/2)
  chain_varL1[2] <- rinvchisq(n = 1, df = 4*nLGM, scale = s_sq)
  
  
  s1 <- rbind(chain_thet[,1,3]-chain_betas[,1,3])%*%cbind((chain_thet[,1,3]-chain_betas[,1,3]))
  s2 <- rbind(chain_thet[,2,3]-chain_betas[,1,3]-chain_betas[,3,3]*chain_thet[,1,3] - chain_betas[,2,3]*1)%*%cbind(chain_thet[,2,3]-chain_betas[,1,3]-chain_betas[,3,3]*chain_thet[,1,3] - chain_betas[,2,3]*1)
  s3 <- rbind(chain_thet[,3,3]-chain_betas[,1,3]-chain_betas[,3,3]*chain_thet[,2,3] - chain_betas[,2,3]*2)%*%cbind(chain_thet[,3,3]-chain_betas[,1,3]-chain_betas[,3,3]*chain_thet[,2,3] - chain_betas[,2,3]*2)
  s4 <- rbind(chain_thet[,4,3]-chain_betas[,1,3]-chain_betas[,3,3]*chain_thet[,3,3] - chain_betas[,2,3]*3)%*%cbind(chain_thet[,4,3]-chain_betas[,1,3]-chain_betas[,3,3]*chain_thet[,3,3] - chain_betas[,2,3]*3)
  s_sq <- (s1+s2+s3+s4)/(4*nLGM)
  #sigL1 <- 1/rgamma(1, 4*nLGM/2, s_sq/2)
  chain_varL1[3] <- rinvchisq(n = 1, df = 4*nLGM, scale = s_sq)
  
  
  #Get level 2 var-covar last.
  # about 37 ms
  
  matSum2 <- matrix(0, nrow = 3, ncol = 3)
  for(i2 in 1:nLGM){
    betVec <- matrix(chain_betas[i2,,1], nrow = 3)
    matSum2 <- matSum2 + (betVec-chain_gammas[,1])%*%(t(betVec-chain_gammas[,1]))
  }
  # Not solve(matSUm2). Nice!
  chain_TMat[,,1] <- riwish(nLGM, matSum2)
  
  matSum2 <- matrix(0, nrow = 3, ncol = 3)
  for(i2 in 1:nLGM){
    betVec <- matrix(chain_betas[i2,,2], nrow = 3)
    matSum2 <- matSum2 + (betVec-chain_gammas[,2])%*%(t(betVec-chain_gammas[,2]))
  }
  chain_TMat[,,2] <- riwish(nLGM, matSum2)

  matSum2 <- matrix(0, nrow = 3, ncol = 3)
  for(i2 in 1:nLGM){
    betVec <- matrix(chain_betas[i2,,3], nrow = 3)
    matSum2 <- matSum2 + (betVec-chain_gammas[,3])%*%(t(betVec-chain_gammas[,3]))
  }
  chain_TMat[,,3] <- riwish(nLGM, matSum2)
  
  #Fixing TMat[1,1] = 1 here.
  chain_TMat[1, ,1] <- chain_TMat[1,,1]*1/sqrt(chain_TMat[1,1,1])
  chain_TMat[,1,1] <- chain_TMat[,1,1]*1/(chain_TMat[1,1,1])
  
  chain_TMat[1,,2] <- chain_TMat[1,,2]*1/sqrt(chain_TMat[1,1,2])
  chain_TMat[,1,2] <- chain_TMat[,1,2]*1/(chain_TMat[1,1,2])
  
  chain_TMat[1,,3] <- chain_TMat[1,,3]*1/sqrt(chain_TMat[1,1,3])
  chain_TMat[,1,3] <- chain_TMat[,1,3]*1/(chain_TMat[1,1,3])
   
  chain_thetDraw[,i,] <- chain_thet
  chain_aDraw[,i,] <- chain_ak
  chain_bDraw[,i,] <- chain_bk
  chain_betasDraw[,i,] <- chain_betas
  chain_varL1Draw[i,] <- chain_varL1
  chain_TMatDraw[,,i,] <- chain_TMat
  chain_gammaDraw[,i,] <- chain_gammas
  
  #Sys.time()
  if(i%%20 == 0) print(i)
}
