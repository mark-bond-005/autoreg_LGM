# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// !!! Do matrix algebra for the level one equations
// each person in the LGM gets his own intercept, slope, and AR mean & covar
// [[Rcpp::export]]
List getMeanVar(arma::mat thet, arma::vec gammas, 
			const float varL1, arma::mat TMat){
	// Variable declaration
	const int nLGM = thet.n_rows;
	const int nTimes = thet.n_cols;
	arma::mat meanOut(nLGM, 3);
	arma::mat DOut(nLGM, 9);
	for (int iLGM = 0; iLGM < nLGM; iLGM++)
	{
		arma::vec thetJ(nTimes);
		arma::mat xTi(3, nTimes);
		// The first column is a pain, so code it manually.
		xTi(0, 0) = 1;
		xTi(1, 0) = 0;
		xTi(2, 0) = 0;
		thetJ(0) = thet(iLGM, 0);
		// Building IV and DV matrices
		for (int iTime = 1; iTime < nTimes; iTime++)
		{
			xTi(0, iTime) = 1;
			xTi(1, iTime) = iTime;
			xTi(2, iTime) = thet(iLGM, iTime-1);
			thetJ(iTime)  = thet(iLGM, iTime);
		}
		arma::mat sumSquares;
		sumSquares = xTi*xTi.t();
		arma::mat betaHat;
		betaHat = sumSquares.i()*xTi*thetJ;
		arma::mat bigSig;
		bigSig = varL1 * sumSquares.i();
		
		arma::vec d;
		d = bigSig.i()*betaHat + TMat.i()*gammas;
		arma::mat D;
		D = (bigSig.i() + TMat.i()).i();
		arma::vec mean;
		mean = D*d;
		
		meanOut(iLGM, 0) = mean(0);
		meanOut(iLGM, 1) = mean(1);
		meanOut(iLGM, 2) = mean(2);

		DOut(iLGM, 0) = D(0);
		DOut(iLGM, 1) = D(1);
		DOut(iLGM, 2) = D(2);
		DOut(iLGM, 3) = D(3);
		DOut(iLGM, 4) = D(4);
		DOut(iLGM, 5) = D(5);
		DOut(iLGM, 6) = D(6);
		DOut(iLGM, 7) = D(7);
		DOut(iLGM, 8) = D(8);
	}
	List ret;
	ret["means"] = meanOut;
	ret["D"] = DOut;
	return(ret);	
}
