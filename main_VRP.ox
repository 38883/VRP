/*
  Model: 5 vector equations (3 states, 2 observed; return include jump; kappa tau << 1)
  BV_t = 0.5 kappa_i theta_i tau + (1 - 0.5kappa_i tau) sig2_(t-tau) + epsilon,
  VIX^2_it = a_it(gamma_it) + b_t(gamma_it) sigma_it + c_it * EJ_t - omega_it + epsilon_it
    a_it = 1 - exp(-(kappa+gamma_it)T) / (T*(gamma_it + kappa)) 
    b_it = theta (1 - a_it)
    EJ_t = (RV_it - BV_it) + epsilon_t
  sig2_t = sig2_(t-tau) + Kappa(theta  - sig2_(t-tau)) tau + nu sig2_(t-tau) epsilon,
  dgamma_t = chol(Sig_gam) * dW_t^gam
  domega_t = hocl(Sig_ome) * dW_t^ome
*/

#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import  <maximize>
#include "Hawkes_jump.ox"
#include "LinearApproxMH.ox"
#include "Heston_vol.ox"

fKalmanFilter(const mY, const amX, const vBeta, const amZ, const amG, const amW,
			const amT, const amF, const mF0, const iRnd);
fTsvar(const vx, const iBm);
fCD(const vx, const iBm);
fMyWishart(const n,const mS);
fvG(const vGamt, const vKappa, const vTheta, const sTau, const sTerm, const vSigt2);


decl giLinear = 0, giKTVIX = 0, giNuWishart = 1,  giSigmaKalman = 1, giSelectRP = 0;
// giLinear = 1 linear for VIX model, 0 otherwise (generally select 0)
// giKTVIX: include VIX equation for the estimation of kappa and theta if 1. generally 0
// giNuWishart = 1 if Nu is Wishart dist, 0 for gamma process (diagonal)
// giSelectRP: select VRP or JRP to caluculate confidence interval: 0 for save gamma, 1 for save omega
decl giJumpSepCalc = 0;
// giJumpSepCalc = 1 if jump parameters are separately calculated from this code.

main(){

	decl mKappa, vTheta, mSigBV2inv, vKappa, sThetaG, mSigBV2, mSigGam2, mSigOme2, mSigGam2inv, mSigOme2inv;
	decl mData, mRV, mRV0, mVIX, mBV, mBV0, vSigVIX2, mOut, mOme, mGam, msigt, mG;	
	decl sT, sD, sPara, sSim, sBurnin, sFactor, i, j, k, t, vsigtProp, sTau, sTerm;
	decl mA0, vA0, vS, mS, sS, mSS, vSS, mB0, vB0, sC0, mC0, sD0, vD0, sE0, mE0, sF0, vF0, sH0, mH0;
	decl amSp, amSp2, vsDate, time, ans, an0, sBM, vCI, maxAm, vAccept, mSigVIX2inv;
	decl mHhat, mHprop, amSigstar, amGHhat0, v0, sTrue, sProp, sTrueOld, sPropOld, vNu;
	decl amX, vBeta, amZ, amG, amW, amT, amF, mF0, mG0, vG0, mI0, vI0;
	decl vKapProp, vThetaProp, vS1, mS1, vS2, mSigNu2, mSigNu2inv, vSigBV2;
	decl amSpAll, vCI2, amSig2Jinv, mLam, vEta, vLamAlp, vLamBet, vLamAve;
	decl mJParams, mJLatent, mEJ, vOmeBarProp, mOmeBar, vOmeProp, aP, sNnods;
	decl vLower, vUpper, mRange, mGreal, mEJreal, mZ, mRVBV, sThres, mRet;
	decl vPhiAR1V, vMuAR1V, vPhiAR1J, vMuAR1J;

	// parameters for sampling
	time = timer();	
	ans = zeros(0,0);
	sBM = 500;
	sSim = 10000;	                 // number of samples
	sBurnin = sSim * 0.1;            // number of burn-in
 	vCI = <.025, .975>;	 			 // confident intervals
	vCI2 = <0.05, 0.95>;
	sNnods = 25;  // 5;                     // number of nods in nonlinear AR sampling
	sThres = 0;

	// data inputs and adjustments
	mData = loadmat(sprint("./input.xlsx"));
	sT = rows(mData);				 // sample period of input
	sD = (columns(mData) - 1 ) / 4;	 // dimenstion (#countries) of RV, BV, IV 
	sPara = 8 + 6;	                 // number of parameters
	maxAm = 8;	                     //  number of latent variables (jump variables are regarded one)

	// extraction of holidays
	for (k = 0; k < sT; k++ ) {
		if (prodr(mData[k][]) == 0 || isnan(prodr(mData[k][]))) {
			if (k == 0)	mData = mData[1:sT-1][];
			else if (k == sT-1) mData = mData[0:k-1][];
			else mData = mData[0:k-1][] | mData[k+1:sT-1][];
			sT -= 1;
			k = k - 1;
		}
	}

	// set variables
	vsDate = mData[1:sT-1][0];
	mRV  = mData[1:sT-1][1:sD] * 252;
	mRV0 = mData[0:sT-2][1:sD] * 252;
	mBV  = mData[1:sT-1][sD+1:sD*2] * 252;
	mBV0 = mData[0:sT-2][sD+1:sD*2] * 252;
	mVIX = (mData[1:sT-1][sD*2+1:sD*3] / 100) .^ 2;
	mRet = mData[1:sT-1][sD*3+1:sD*4];
	mRet = (mRet .> 0) .? 1 .: -1;
	sT = sT - 1;	
	mRVBV = mRV - mBV;

	sTerm = 1/12;   // one month 
	sTau = 1/252;   // one-day annualization
	
	mHprop = zeros(sT, sD);
	vAccept = zeros(3,1);	
	
	// initial value of parameters
	vKappa = 5 * ones(sD, 1);
	mKappa = diag(vKappa);
	vTheta = 0.1 * ones(sD, 1);
	vSigBV2 = ones(sD, 1);
	mSigBV2inv = invert(diag(vSigBV2));
	vSigVIX2 = 0.1 * ones(sD, 1);
	mSigVIX2inv = invert(diag(vSigVIX2));
	mSigGam2 = mSigGam2inv = 0.1 * unit(sD);
	mSigOme2 = 0.1 * unit(sD);
	mSigOme2inv = invert(mSigOme2);
	
	vNu = 0.1 * ones(sD, 1);
	mSigNu2 = diag(vNu);
	mSigNu2inv = invert(mSigNu2);
	if (giJumpSepCalc == 1)	{
		mJParams = loadmat("out_jump_params.xlsx");
		mJLatent = loadmat("out_jump_latent.xlsx");
	} else {
		mJParams = ones(sD, 8);
		mJParams[][0] = 0.01 * ones(sD, 1);
		mJParams[][2] += 8 * ones(sD, 1); 
   		mJLatent = 0.5 * ones(sT, sD*4);
	}
	
	// initial value of latent variables
	mGam = -10 * ones(sT, sD);			  // DRP
	mOme = -20 * ones(sT, sD);			  // JRP
	msigt = 0.01 * ones(sT, sD);		  // sigma_t
	mG = mGreal = zeros(sT, sD);
	mEJ = mEJreal = zeros(sT, sD);        // EJ

	vMuAR1V = vMuAR1J = 0.1 * ones(sD, 1);
	vPhiAR1V = vPhiAR1J = ones(sD, 1);
	
	
	//priors
	vA0 = 3 * ones(sD, 1);  // mean for kappa
	mA0 = 1 * unit(sD);   	  // covar for kappa
	vB0 = 0.1 * ones(sD, 1); // mean for theta
	mB0 = 0.01 * unit(sD);		  // covar for theta
	sC0 = 10;				  // dof for SigmaBV
	mC0 = 1 * unit(sD);	  	  // scale for SigmaBV
	sD0 = sD+1;				  // dof for sigma_VIX
	vD0 = 0.01 * ones(sD, 1);    // vol for sigma_VIX 
	sE0 = sD+1;				  // dof for sig_gam
	mE0 = 0.1 * unit(sD);		  // scale for Sig_gam
	sH0 = sD+1;					  // dof for sig_ome
	mH0 = 0.1 * unit(sD);		  // scale for sig_ome 
	sF0 = 10;
	vF0 = 1 * ones(sD, 1);
	mG0 = 0.0001 * unit(sD);	  // covar for gamma
	vG0 = -1 * ones(sD, 1);	  // mean for gamma
	mI0	= -0.01 * unit(sD);	   // prior for omega
	vI0 = zeros(sD, 1);		   // prior for omega

	// variables for filtering
	amX = amZ = amG = amW = amT = amF =  new array[sT];

	// storing parameters initial setting
	amSp = new array[maxAm];
	for (i = 0; i < 3; i++) amSp[i]  = zeros(sT, sD);
	amSp[3]	= zeros(sT, sD*4);  // jump related latent variable (mJ ~ mX ~ mLam)
	amSp[4] = amSp[5] = zeros(sT, sD);
	amSp[6] = amSp[7] = zeros(sT, sD*2);
	
	
	amSp2 = new array[sPara];
	for (i = 0; i < sPara;     i++) amSp2[i] = zeros(sSim, sD);
	
	amSpAll = new array[sD * 2];
	for (i = 0; i < sD; i++) amSpAll[i] = zeros(sSim, sT);

	// saving inputs
	savemat(sprint("out_input.xlsx"), vsDate ~ mRV ~ mBV ~ mVIX);


  // Gibbs sampling starts
  for (k = -sBurnin; k < sSim; k++) {
	 
  // sampling Kappa	
  	mS = invert(mA0);
	vS = mS * vA0;
	for (t = 1; t < sT; t++) {
		mSS = diag(msigt[t-1][]' - vTheta) * 0.5 * sTau;
		mS += mSS ' mSigBV2inv * mSS;
		vS += mSS ' mSigBV2inv * (msigt[t-1][]' - mBV[t][]');
		mS += 4 * mSS ' mSigNu2inv * mSS;
		vS += 2 * mSS ' mSigNu2inv * (msigt[t-1][]' - msigt[t][]');
		if (giKTVIX == 1) {
			mSS = diag(msigt[t][]' - vTheta) * 0.5 * sTerm;
			mS += mSS ' mSigVIX2inv * mSS;
			vS += mSS ' mSigVIX2inv * (msigt[t][]' - 0.5 * mGam[t][]' .* msigt[t][]' * sTerm
					+ mEJ[t][]' - mVIX[t][]');		
		}
	}
	mS = invert(mS);
	vS = mS * vS;
	vKapProp = vS + choleski(mS) * rann(sD, 1);
	
	if (min(vKapProp) <= 0)	for (i = 0; i < sD; i++) if(vKapProp[i] < 0) vKapProp[i] = 0;	
	if (giKTVIX == 1) {
		sTrue  = sTrueOld = sProp = sPropOld = 0;
		for (t = 0; t < sT; t++) {
			vS = mVIX[t][]' - fvG(mGam[t][]', vKapProp, vTheta, sTau, sTerm, msigt[t][]')
				 - mEJ[t][]';	  
			sTrue    += -0.5 * vS ' mSigVIX2inv * vS;
			vS = mVIX[t][]' - fvG(mGam[t][]', vKappa  , vTheta, sTau, sTerm, msigt[t][]')
				 - mEJ[t][]';	  
			sTrueOld += -0.5 * vS ' mSigVIX2inv * vS;
		
			vS = mVIX[t][]' - (ones(sD, 1) - 0.5 * (vKapProp + mGam[t][]') * sTerm) .* msigt[t][]'
					- 0.5 * vKapProp .* vTheta * sTerm - mEJ[t][]';	  
			sProp    += -0.5 * vS ' mSigVIX2inv * vS;
			vS = mVIX[t][]' - (ones(sD, 1) - 0.5 * (vKappa   + mGam[t][]') * sTerm) .* msigt[t][]'
					- 0.5 * vKappa 	 .* vTheta * sTerm - mEJ[t][]';	 		
			sPropOld += -0.5 * vS ' mSigVIX2inv * vS;	   
		}

		if (ranu(1,1) < exp(sTrue - sProp - sTrueOld + sPropOld)) {
			vKappa = vKapProp;
			if (k >= 0) vAccept[0]++;
			println("accept Kappa", k, vKappa');
		}	
	} else {
		if (min(vKapProp) <= 0)	for (i = 0; i < sD; i++) if(vKapProp[i] < 0) vKapProp[i] = 0;
		vKappa = vKapProp;
	}
	mKappa = diag(vKappa);

	// sampling Theta
	mS = invert(mB0);
	vS = mS	* vB0;
	for (t = 1; t < sT; t++) {
		mS += (sTau * 0.5) ^ 2 * mKappa ' mSigBV2inv * mKappa;
		vS += 0.5 * sTau * mKappa ' mSigBV2inv
			  * (mBV[t][]' - (unit(sD) - 0.5 * sTau * mKappa) * msigt[t-1][]');
		mS += sTau * mKappa ' mSigNu2inv * mKappa * sTau;
		vS += sTau * mKappa ' mSigNu2inv
			  * (msigt[t][]' - (unit(sD) - mKappa * sTau) * msigt[t-1][]');
		if (giKTVIX == 1) {
			mS += 0.5 * sTerm * mKappa ' mSigVIX2inv * mKappa * sTerm * 0.5;
			vS += 0.5 * sTerm * mKappa ' mSigVIX2inv
					* (mVIX[t][]' - mEJ[t][]'
					   - (ones(sD, 1) - 0.5 * (vKappa + mGam[t][]') * sTerm) .* msigt[t][]');		
		}
	}
	mS = invert(mS);
	vS = mS * vS;
	vThetaProp = vS + choleski(mS) * rann(sD, 1);
	if (min(vThetaProp) < 0) for (i = 0; i < sD; i++) if(vThetaProp[i] < 0) vThetaProp[i] = 0;
	if (giKTVIX == 1) {
		sTrue = sTrueOld = 0;
		for (t = 0; t < sT; t++) {
			vS = mVIX[t][]' - fvG(mGam[t][]', vKappa, vThetaProp, sTau, sTerm, msigt[t][]')
				 -  mEJ[t][]';	  
			sTrue    += -0.5 * vS ' mSigVIX2inv * vS;
			vS = mVIX[t][]' - fvG(mGam[t][]', vKappa, vTheta    , sTau, sTerm, msigt[t][]')
				 -  mEJ[t][]';	  
			sTrueOld += -0.5 * vS ' mSigVIX2inv * vS;
		
			vS = mVIX[t][]' - (ones(sD, 1) - 0.5 * (vKappa + mGam[t][]') * sTerm) .* msigt[t][]'
				 -0.5 * vKapProp	.* vThetaProp * sTerm - mEJ[t][]';	  
			sTrue    += +0.5 * vS ' mSigVIX2inv * vS;
			vS = mVIX[t][]' - (ones(sD, 1) - 0.5 * (vKappa + mGam[t][]') * sTerm) .* msigt[t][]'
				 -0.5 * vKappa   	.* vTheta     * sTerm - mEJ[t][]';	 		
			sTrueOld += +0.5 * vS ' mSigVIX2inv * vS;
		}
		if (ranu(1,1) < exp(sTrue - sTrueOld)) {
			vTheta = vThetaProp;
			if (k >= 0) vAccept[1]++;
		}
	} else {
		if (min(vThetaProp) < 0) for (i = 0; i < sD; i++) if(vThetaProp[i] < 0) vThetaProp[i] = 0;
		vTheta = vThetaProp;
	}
 	
	// sampling SigmaBV;
	 for (i = 0; i < sD; i++) {
		sS = sF0 + sT - 1;
		mS = vF0[i] + 0.5 * sumsqrc((mBV[1:sT-1][i]
					- (1 - mKappa[i][i] * sTau * 0.5) * msigt[0:sT-2][i]
					- mKappa[i][i] * vTheta[i] * sTau * 0.5));
			vSigBV2[i] = 1 / rangamma(1, 1, sS / 2, mS);
	}
	mSigBV2 = diag(vSigBV2);
	mSigBV2inv = invert(mSigBV2);	
	
	// sampling vSigma_it
	for (t = 0; t < sT; t++) {
		amX[t] = mKappa * vTheta * sTau ~ zeros(sD, 1);
		amZ[t] = unit(sD) - mKappa * sTau * 0.5;
		amG[t] = choleski(mSigBV2) ~ zeros(sD, sD);
		amW[t] = zeros(sD, 1) ~ mKappa * vTheta * sTau;
		amT[t] = unit(sD) - mKappa * sTau;
		amF[t] = zeros(sD, sD) ~ choleski(mSigNu2);
	 }
	vBeta = <0.5; 1>;
	mF0 = amF[0];
	mOut = fKalmanFilter(mBV, amX, vBeta, amZ, amG, amW, amT, amF, mF0, 1)[0:sT][];
	msigt[0:sT-2][] = mOut[2:sT][];
	msigt[sT-1][] = ((unit(sD)-mKappa*sTau)*msigt[sT-2][]'+mKappa*vTheta*sTau
					+choleski(mSigNu2)*rann(sD, 1))'; 
	msigt = (msigt .> 0) .? msigt .: 0;
					
	// sampling nu
	if (giNuWishart == 1) {
		sS  = sC0 + sT - 1;
		mS = invert(mC0);
		for (t = 1; t < sT; t++) {
			vSS = msigt[t][]' - (unit(sD) - mKappa * sTau) * msigt[t-1][]' 
				   - mKappa * vTheta * sTau; 
			mS += vSS * vSS';
		}
		mS = invert(mS);
		mSigNu2inv = fMyWishart(sS, mS);	
		mSigNu2 = invert(mSigNu2inv);

		for (i = 0; i < sD; i++) {
			sS = sF0 + sT - 1;
			mS = vF0[i] + 0.5 * sumsqrc((msigt[1:sT-1][i]
						- (1 - mKappa[i][i] * sTau) * msigt[0:sT-2][i]
						- mKappa[i][i] * vTheta[i] * sTau));
			vNu[i] = 1 / rangamma(1, 1, sS / 2, mS);
		}
		mSigNu2 = diag(vNu);
		mSigNu2inv = invert(mSigNu2);
	}
	
	// sampling vSigVIX
	
	sS = sD0 + sT;
	mS = vD0 + sumsqrc(mVIX - mG - mEJ)';
	for (i = 0; i < sD; i++) vSigVIX2[i] = 1 / rangamma(1, 1, sS/2, mS[i]/2);
	
	// sampling DRP (gamma)
	
	vLower = -150 * ones(sD,1);
	vUpper = zeros(sD, 1);
	mRange = vLower ~ vUpper;
	for (i = 0; i <sD; i++) {
		mGreal[][i] = fv_G(zeros(sT, 1), fvSetaP_G(zeros(sT, 1), msigt, vKappa, vTheta, sTerm, i, zeros(sT, sD))); 
	}	
	aP = fvSetaP_G(mGam, msigt, vKappa, vTheta, sTerm, 0, mGreal);

	mGam = fLinearApproxARMH(mVIX - mEJ, mGam, fvGt, fmLAIntercG, fmLACoeffG,
				aP, &vMuAR1V, &vPhiAR1V, diag(vSigVIX2), &mSigGam2, mRange, 20);
				
	for (i = 0; i <sD; i++) 
		mGreal[][i] = fv_G(zeros(sT, 1), fvSetaP_G(zeros(sT, 1), msigt, vKappa, vTheta, sTerm, i, zeros(sT, sD))); 
	
	mG = mGreal + fmLACoeffG(aP) .* mGam;
	
	
	// sampling jump
	if (giJumpSepCalc != 1) {
		[mJLatent, mJParams] = faHawkesMertonJump(mRVBV, mJLatent, mJParams);
	}

	// sampling JRP	(omega)
	vUpper = zeros(sD, 1);
	vLower = -500 * ones(sD, 1);
	mRange = vLower ~ vUpper;
 
	mEJreal = fmHawkesEJ(zeros(sT, sD), mJLatent, mJParams, sTerm);
 	aP = faSetParams(mJLatent, mJParams, sTerm, 0, mEJreal);
	
	mOme = fLinearApproxARMH(mVIX - mG, mOme, fvEJt, fmLAIntercEJ, fmLACoeffEJ, aP,
					&vMuAR1J, &vPhiAR1J, diag(vSigVIX2), &mSigOme2, mRange, 20);
   
	mEJ = mEJreal + fmLACoeffEJ(aP) .* mOme;

	// storing 
	 if (k >= 0) {
	 
		amSp[0] += mGam / sSim;
		amSp[1] += mOme / sSim;
		amSp[2] += msigt / sSim;
		amSp[3] += mJLatent ./ sSim;
		amSp[4] +=  mG ./sSim;
		amSp[5] +=  mEJ ./ sSim;
		amSp[6] += (mGreal - mG) ./sSim
				~ (mEJreal - mEJ) ./sSim;
		amSp[7] +=  mGam .* msigt /sSim
				~ mOme .* mJLatent[][sD*2:sD*3-1] ./ mJLatent[][sD*3:sD*4-1] / sSim;
				
		amSp2[0][k][] = diagonal(mKappa);
		amSp2[1][k][] = vTheta';
		amSp2[2][k][] = vSigBV2';
		amSp2[3][k][] = vSigVIX2';
		for (i = 1; i < 8; i++) amSp2[i+3][k][] = mJParams[][i]';
		amSp2[11][k][] = diagonal(mSigNu2);
		amSp2[12][k][] = diagonal(mSigGam2);
		amSp2[13][k][] = diagonal(mSigOme2);
	
														 	  
	}
	
	 if (!fmod(k, sSim/100)) println("Iteration:", "%5.0f", k, "  Time:", "%4.1f", timespan(time), "\n");

   }  // sampling end						 
 
	// save outputs	: DRP(Gamma) ~ JRP(Omega) ~ vol(sigma2) ~ JumpIndicator ~ jumpSize ~ Intensity
	savemat(sprint("out_LatentVariables.xlsx"), amSp[0] ~ amSp[1] ~ amSp[2] ~ amSp[3]);
	savemat(sprint("out_VRPs.xlsx"), amSp[6] ~ amSp[4] ~ amSp[5]);
	savemat(sprint("residuals.xlsx"), mVIX - amSp[4] - amSp[5] ~ amSp[7] ~ amSp[4] + amSp[7][][0:sD-1] ~ amSp[5] + amSp[7][][sD:2*sD-1]);
											
	// stats of samples 					      
	for (i = 0; i < sPara; i++)	{	    	
		println(meanc(amSp2[i]));	        
		an0 = zeros(0,0);
		if (i < 12) {
			for (j = 0; j < sD; j++)
				an0 = an0 ~ fCD(amSp2[i][][j], sBM);
		} else {
			for (j = 0; j < sD; j++) 
				an0 = an0 ~ fCD(amSp2[i][][j], sBM);
		}
			
		ans = ans ~
		(meanc(amSp2[i])
		 | sqrt(varc(amSp2[i]))
		 | quantilec(amSp2[i], vCI[0])
		 | quantilec(amSp2[i], vCI[1])
		 | an0);	
	}

	savemat(sprint("Out_StatsParams.xlsx"), ans);
}

fTsvar(const vx, const iBm)
{
	decl ns, vsp;

	ns = rows(vx);
	vsp = periodogram(vx, iBm, 2, 1) / ns;
	
	return M_2PI * vsp[0];
}


fCD(const vx, const iBm)
{
	decl ns, n1, n2, vx1, vx2, dx1bar, dx2bar;
	decl dvar1, dvar2, dz;

	ns = rows(vx);
	n1 = floor(0.1*ns);
	n2 = floor(0.5*ns);
	vx1 = vx[:n1-1];
	vx2 = vx[ns-n2:];
	dx1bar = meanc(vx1);
	dx2bar = meanc(vx2);
	dvar1 = fTsvar(vx1, iBm);
	dvar2 = fTsvar(vx2, iBm);
	dz = (dx1bar - dx2bar) / sqrt(dvar1 / n1 + dvar2 / n2);
	
	return 2 * tailn(fabs(dz));
}

fKalmanFilter(const mY, const amX, const vBeta, const amZ, const amG,
			const amW, const amT, const amF, const mF0, const iRnd)
{
	// yt = Xt beta + Zt alp_t + Gt u_t,
	// alp_t+1 = Wt beta + Tt alp_t + Ft u_t 
	// input:  mY:T*C, amX:T*C*K, vBeta:K*1, amZ:T*C*S, amG:T*C*L, amW:T*S*K, amT:T*S*S, amF:T*S*L, mF0:S*L
	// return:  (T+1)*S matrix 
	// T:#observations, C:#variables, S:#states, L:#stochastic Factors, K:#constants
	// iRnd: 1 include randomness in simulation smoother; 0 estimate expectation of the smoother
	
	decl i, sT, sC, sS, ma, me, amDinv, va, mP, mK, amL, mJ, mC, mPhi, vep, vr, mU, mV, meta, veta0, vout;

	sT = rows(mY);
	sC = rows(amZ[0]);
	sS = columns(amZ[0]);
	me = zeros(sT, sC);
	amL = amDinv = new array[sT];
	meta = zeros(sT, sS);
	vout = zeros(sT+1, sS);
	
	//Kalman filter	
	va = amW[0] * vBeta;
	mP = mF0 * mF0';
	for (i = 0; i < sT; i++){
		me[i][] = (mY[i][]' - amX[i] * vBeta - amZ[i] * va)';
		amDinv[i] = invert(amZ[i]*mP*amZ[i]' + amG[i]*amG[i]');
		mK = (amT[i] * mP * amZ[i]' + amF[i] * amG[i]') * amDinv[i];
		va = amW[i] * vBeta + amT[i] * va + mK * me[i][]';
		amL[i] = amT[i] - mK * amZ[i];
		mJ = amF[i] - mK * amG[i];		  
		mP = amT[i] * mP * amL[i]' + amF[i] * mJ';
	}

	// simulation smoother
	vr = zeros(sS, 1);
	mU = zeros(sS, sS);
	for (i = sT-1; i >= 0; i--) {
		mPhi = amF[i] * amF[i]';
		mC = mPhi - mPhi * mU * mPhi;
		if (iRnd == 0) vep = zeros(sS, 1);
		else vep = choleski(mC) * rann(sS, 1);
		meta[i][] = (mPhi * vr + vep)';
		mV = mPhi * mU * amL[i];
		vr = amZ[i] ' amDinv[i] * me[i][]' + amL[i] ' vr - mV ' invert(mC) * vep;
		mU = amZ[i] ' amDinv[i] * amZ[i] + amL[i] ' mU * amL[i] + mV ' invert(mC) * mV;
	
	}
	mPhi = mF0 * mF0';

	if (iRnd == 0) {
		veta0 = mPhi * vr;
	} else {
		mC = mPhi - mPhi * mU * mPhi;
		vep = choleski(mC) * rann(sS, 1);
		veta0 = mPhi * vr + vep;
	}
	
	vout[0][] = amW[0] * vBeta + veta0;
	for (i = 0; i < sT; i++)
		vout[i+1][] = (amW[i] * vBeta + amT[i] * vout[i][]' + meta[i][]' )';
	
	return vout;	
	
}
	
fMyWishart(const n,const mS){
// X~W(n,S), X:(pxp) symmetrx matrix
// f(X|n,S) = const *|X|^{(n-p-1)/2}*exp-0.5*tr(S^{-1}X)
// n: degrees of freedom 
// S: (p x p) symmetric parameter matrix
// E(X_{ij})=n*s_{ij}, s_{ij}: (i,j) element of S.
	decl ci, cj, cp, mA, mL, mT, mX;
	cp = rows(mS);
	mT = zeros(cp,cp);
	if(n < cp)
		print("d.f = ", n, " size = ", cp, " d.f. too small.\n");
	for(ci=0; ci<cp; ++ci) {
		mT[ci][ci] = sqrt(ranchi(1,1,n-ci));
		for(cj=0; cj<ci; ++cj) 
			mT[ci][cj]=rann(1,1);
	}
	mA = mT*mT'; // A=TT' ~W(n,I)
	mL = choleski(mS); // S=LL'
	if (mL==0)
		println(mS);
	mX = mL*mA*mL'; // X=LAL'~W(n,S)
	return(mX);
}
