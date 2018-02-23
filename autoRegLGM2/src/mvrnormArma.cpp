# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// !!! draw multivariate samples
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
	int ncols = sigma.n_cols;
	arma::mat Y = arma::randn(n, ncols);
	return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
