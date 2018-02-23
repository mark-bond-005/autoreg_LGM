# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// !!! Implement the forward filter (Kalman filter)
// !!! The model is:
// !!! thetVec_t = [1, thet[t]]
// !!! F = [ak, bk]
// !!! Then Zs = F*thetVec + vt   vt the k*k identity
// !!! And if G = [[1, 0][alpha+beta1*t, beta2]]
// !!! Then thetVec_t = G*thetVec_t-1 + wt   wt =[[0, 0][0, varL1]]
// !!! Takes person-level data and parameters as input and returns
// !!! prior and posterior means and variances for their theta values.

// [[Rcpp::export()]]
List forFilt(const arma::mat thet, const arma::mat zMat,
	const arma::vec alpha, const arma::vec beta1, const arma::vec beta2,
	const arma::vec ak, const arma::vec bk,
	const arma::vec m0, const float varL1, const float c0){
	// !!! Declare and initialize variables we will loop through
	// Get indices
	const int nLGM = thet.n_rows;
	const int nTimes = thet.n_cols;
	const int nItems = ak.n_elem;
	// Initialize the forward-filter matrices at zero
	// m and c posterior mean and variance
	// A and r prior mean and varaince
	arma::mat mMat(nLGM, nTimes);
	arma::mat cMat(nLGM, nTimes);
	arma::mat rMat(nLGM, nTimes);
	arma::mat AMat(nLGM, nTimes);
	mMat.fill(0);
	cMat.fill(0);
	rMat.fill(0);
	AMat.fill(0);

	// !!! Declare and initialize variables that do not change 
	// !!! during the filtering process. Item information and 
	// !!! variance within the filter does not change (homoscedasticity)
	// Set up the item vector matrix F 
	arma::mat F(nItems, 2);
	for (int iItem = 0; iItem < nItems; iItem++)
	{
		F(iItem, 0) = -1 * bk(iItem);
		F(iItem, 1) = ak(iItem);
	}
	// Set vt. Item-level variance assumed to be the identity.
	arma::mat vt(nItems, nItems);
	vt.eye();
	// Set wt. Theta variance assumed to be varL1. 
	// The constant 1 in thetvec never changes
	arma::mat wt;
	wt << 0 << 0 << arma::endr
	   << 0 << varL1 << arma::endr;

	// !!! Begin the Kalman Filter
	for (int iLGM = 0; iLGM < nLGM; iLGM++)
	{
		//Populate the first column of the matrices only
		mMat(iLGM, 0) = m0(iLGM);
		cMat(iLGM, 0) = c0;
		rMat(iLGM, 0) = beta2(iLGM)*beta2(iLGM)*c0 + varL1;
		AMat(iLGM, 0) = alpha(iLGM) + beta1(iLGM) + beta2(iLGM)*m0(iLGM);
		arma::mat zRow(1, nItems);
		zRow = zMat.row(iLGM);
		for (int iTime = 1; iTime < nTimes; iTime++)
		{
			// Populate the outcome vector for each person and time pt
			arma::mat zij(1, nItems);
			int zLow = (iTime)*nItems;
			int zHigh = (iTime+1)*nItems - 1;
			zij = zRow.cols(zLow, zHigh);
			// Set G depending on the item parameters;
			arma::mat G(2, 2);
			G << 1 << 0 << arma::endr
				<< alpha(iLGM) + beta1(iLGM)*(iTime) << beta2(iLGM) << arma::endr;
			// Get posteriors for t-1
			arma::vec mOld;
			mOld << 1 << mMat(iLGM, iTime - 1);

			arma::mat cOld;
			cOld << 0 << 0 << arma::endr
				<< 0 << cMat(iLGM, iTime - 1) << arma::endr;
			// Get the predictive mean & var for t
			arma::mat at;
			at = G*mOld;
			AMat(iLGM, iTime) = at(1, 0);

			arma::mat Rt;
			Rt = G*cOld*G.t() + wt;
			rMat(iLGM, iTime) = Rt(1, 1);
			// Get predictive mean and var, calc error
			arma::mat ft;
			ft = F*at;
			arma::mat et;
			et = zij.t() - ft;

			arma::mat Qt;
			Qt = F*Rt*F.t() + vt;

			// Get posterior mean and var
			arma::mat mt;
			mt = at + Rt*F.t()*Qt.i()*et;
			mMat(iLGM, iTime) = mt(1, 0);

			arma::mat Ct;
			Ct = Rt - Rt*F.t()*Qt.i()*F*Rt;
			cMat(iLGM, iTime) = Ct(1, 1);
		}
	}

	// !!! Start returning values
	List ret;
	ret["nLGM"] = nLGM;
	ret["nTimes"] = nTimes;
	ret["mMat"] = mMat;
	ret["cMat"] = cMat;
	ret["rMat"] = rMat;
	ret["AMat"] = AMat;
	return(ret);
}

