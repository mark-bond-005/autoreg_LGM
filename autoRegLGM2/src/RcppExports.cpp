// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;


// mvrnormArma
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma);
RcppExport SEXP autoRegLGM2_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// getMeanVar
List getMeanVar(arma::mat thet, arma::vec gammas, const float varL1, arma::mat TMat);
RcppExport SEXP autoRegLGM2_getMeanVar(SEXP thetSEXP, SEXP gammasSEXP, SEXP varL1SEXP, SEXP TMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type thet(thetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gammas(gammasSEXP);
    Rcpp::traits::input_parameter< const float >::type varL1(varL1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type TMat(TMatSEXP);
    rcpp_result_gen = Rcpp::wrap(getMeanVar(thet, gammas, varL1, TMat));
    return rcpp_result_gen;
END_RCPP
}
// forFilt
List forFilt(const arma::mat thet, const arma::mat zMat, const arma::vec alpha, const arma::vec beta1, const arma::vec beta2, const arma::vec ak, const arma::vec bk, const arma::vec m0, const float varL1, const float c0);
RcppExport SEXP autoRegLGM2_forFilt2(SEXP thetSEXP, SEXP zMatSEXP, SEXP alphaSEXP, SEXP beta1SEXP, SEXP beta2SEXP, SEXP akSEXP, SEXP bkSEXP, SEXP m0SEXP, SEXP varL1SEXP, SEXP c0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type thet(thetSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type zMat(zMatSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type ak(akSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type bk(bkSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type m0(m0SEXP);
    Rcpp::traits::input_parameter< const float >::type varL1(varL1SEXP);
    Rcpp::traits::input_parameter< const float >::type c0(c0SEXP);
    rcpp_result_gen = Rcpp::wrap(forFilt2(thet, zMat, alpha, beta1, beta2, ak, bk, m0, varL1, c0));
    return rcpp_result_gen;
END_RCPP
}
