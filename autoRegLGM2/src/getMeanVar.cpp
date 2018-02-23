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
		// The first column is a pain in the ass, so do it here.
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
		
		// DONT JUDGE ME SENPAI!!!
		// ;.;
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
	const float varL1, const float m0, const float c0){
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
		F(iItem, 0) = -1*bk(iItem);
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
		mMat(iLGM, 0) = m0;
		cMat(iLGM, 0) = c0;
		rMat(iLGM, 0) = beta2(iLGM)*beta2(iLGM)*c0 + varL1;
		AMat(iLGM, 0) = alpha(iLGM) + beta2(iLGM)*m0;
		arma::mat zRow(1, nItems);
		zRow = zMat.row(iLGM);
		for (int iTime = 1; iTime < nTimes; iTime++)
		{
			// Populate the outcome vector for each person and time pt
			arma::mat zij(1, nItems);
			int zLow = (iTime - 1)*nItems;
			int zHigh = iTime*nItems-1;
			zij = zRow.cols(zLow, zHigh);
			// Set G depending on the item parameters;
			arma::mat G(2, 2);
			G << 1 << 0 << arma::endr
				<< alpha(iLGM) + beta1(iLGM)*(iTime - 1) << beta2(iLGM) << arma::endr;
			// Get posteriors for t-1
			arma::vec mOld; 
			mOld << 1 << mMat(iLGM, iTime - 1);

			arma::mat cOld;
			cOld << 0 << 0 << arma::endr
				<< 0 << cMat(iLGM, iTime - 1) << arma::endr;
			// Get the prior mean & var for t
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
			et = ft - zij.t();

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

// [[Rcpp::export()]]
List forFilt2(const arma::mat thet, const arma::mat zMat,
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
