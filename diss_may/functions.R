library(Rcpp); library(RcppArmadillo); library(truncnorm);library(MASS); library(geoR); library(mirt); library(lme4); library(tmvtnorm); library(devtools)
library(autoRegLGM2); 

# Taken from MCMCpack
# Sadly cant figure out how to install it on stampede.
riwish <- function (v, S) 
{
  return(solve(rwish(v, solve(S))))
}

rwish <- function (v, S) 
{
  if (!is.matrix(S)) 
    S <- matrix(S)
  if (nrow(S) != ncol(S)) {
    stop(message = "S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop(message = "v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z[rep(p * pseq, pseq) + unlist(lapply(pseq, seq))] <- rnorm(p * 
                                                                  (p - 1)/2)
  }
  return(crossprod(Z %*% CC))
}