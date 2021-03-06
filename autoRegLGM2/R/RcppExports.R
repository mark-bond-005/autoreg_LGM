# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

mvrnormArma <- function(n, mu, sigma) {
    .Call('autoRegLGM2_mvrnormArma', PACKAGE = 'autoRegLGM2', n, mu, sigma)
}

getMeanVar <- function(thet, gammas, varL1, TMat) {
    .Call('autoRegLGM2_getMeanVar', PACKAGE = 'autoRegLGM2', thet, gammas, varL1, TMat)
}

forFilt <- function(thet, zMat, alpha, beta1, beta2, ak, bk, varL1, m0, c0) {
    .Call('autoRegLGM2_forFilt', PACKAGE = 'autoRegLGM2', thet, zMat, alpha, beta1, beta2, ak, bk, varL1, m0, c0)
}

forFilt2 <- function(thet, zMat, alpha, beta1, beta2, ak, bk, m0, varL1, c0) {
    .Call('autoRegLGM2_forFilt2', PACKAGE = 'autoRegLGM2', thet, zMat, alpha, beta1, beta2, ak, bk, m0, varL1, c0)
}

