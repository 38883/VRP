#import  <maximize>
extern fMyWishart(const n,const mS);
extern fKalmanFilter(const mY, const amX, const vBeta, const amZ, const amG,
			const amW, const amT, const amF, const mF0, const iRnd);

faHawkesMertonJump(const mRVBV, const mLatent, const mParams);
fvC(const vAlp, const vBeta, const sTau, const vEz2);
fvD(const vAlp, const vBeta, const vLamAve, const sTau, const vEz2);
fvEJ(const mOme, const aP);
fvd1EJ(const mOme, const aP);
fvd2EJ(const mOme, const aP);
faGetParams(const mOme, const aP);
faSetParams(const mLatent, const mParams, const sTerm, const sDim, const mG);
fvSetaP(const aPfull, const ir, const ird);
fmHawkesEJ(const mOme, const mLatent, const mParams, const sTerm);

faHawkesMertonJump(const mRVBV, const mLatent, const mParams) {

	decl sTau, sTerm;
	decl mOut, mParamOut;	
	decl sT, sD, i, j, k, t;
	decl vS, vSS, mS, mSS, mS0, vS0, sS, sSS, sF0, mF0, sG0, mG0;
	decl sJumpThres, vERVBV, vVRVBV, vERVBVall, vVRVBVall, mEz2, valp0, vbet0; 
	decl mX, vMuJ, mj, SigG2, SigG2inv;
	decl mPrJ, mJast, vPr, VF0, mLam, mEta;
	decl sNnod, sNnodEta, sDatesBlock, vNod, vNumSp, vNumJp, sNumNod;
	decl vLamAlp, vLamBeta, vLamAve, mSigLam, mSigLaminv;
	decl vLamAlpA0, vLamAlpB0, vLamBetaA0, vLamBetaB0, vLamAveA0, vLamAveB0;
	decl vLamA0, vXi, mZeta2, vXi0, mXi0, mZeta20, vZeta20, mSigZ2, mSigZ2inv, mZeta2inv;
	decl amX, vBeta, amZ, amG, amW, amT, amF, mQ0, mKalmanOut;
	decl mSigRVBV2, mSigRVBV2inv, sH0, mH0, sWindowSize, sWindow, sWStt;
	
	ranseed(today());
	
	sT = rows(mRVBV);
	sD = columns(mRVBV);	 // dimenstion

 	amX = amZ = amG = amW = amT = amF =  new array[sT];	

	mj   = mLatent[][0   :sD  -1];
	mX   = mLatent[][sD  :sD*2-1];
	mLam = mLatent[][sD*2:sD*3-1];
	mEta = mLatent[][sD*3:sD*4-1];
	
	
	vMuJ 	 = zeros(sD, 1);
	vLamAlp	 =      mParams[][2];
	vLamBeta =      mParams[][3];
	vLamAve  =      mParams[][4];
	mSigZ2   = diag(mParams[][5]);
	mSigLam  = diag(mParams[][6]);
	mSigRVBV2 = diag(mParams[][7]);
	mSigZ2inv = invert(mSigZ2);
	mSigRVBV2inv = invert(mSigRVBV2);
	mSigLaminv = invert(mSigLam);

	sJumpThres = 0.15^2;

	mF0 = unit(sD);
	sF0 = sD+1;
	sG0 = sH0 = 10;
	mG0 = 0.01 * unit(sD);
	mH0 = 0.01 * unit(sD);	 
	vS0 = 0.01 * ones(sD, 1);
	vLamAlpA0 = 5 * ones(sD, 1);
	vLamAlpB0 = 5 * ones(sD, 1);
	vLamBetaA0 = 1 * ones(sD, 1);
	vLamBetaB0 = 1 * ones(sD, 1);
	vLamAveA0 = 0.01 * ones(sD, 1);
	vLamAveB0 = ones(sD, 1);
	vLamA0 = 5/252 * ones(sD, 1);	

	vERVBVall = sumc(mRVBV      .* (mRVBV .> sJumpThres))' ./ sumc(mRVBV .> sJumpThres)';
	vVRVBVall = sumc(mRVBV .^ 2 .* (mRVBV .> sJumpThres))' ./ sumc(mRVBV .> sJumpThres)'
				- vERVBVall .^ 2;	
				
	sDatesBlock = 21*6;
	sTerm = 1/252;
	sTau = 1/252;

	sNnod = floor(sT / sDatesBlock);

	sWindowSize = 21; // average window size for the parameters in eta prior

	// sampling X (jump size)	
	for (t = 0; t < sT; t++) {
		for (i = 0; i < sD; i++) {
			if (mj[t][i] == 0) {
				mX[t][i] = max(mRVBV[t][i],0);
			} else {
				sS =  mRVBV[t][i] - mSigZ2[i][i] * mEta[t][i];
				mX[t][i] = max(sS + sqrt(mSigZ2[i][i]) * rann(1,1),0);
			}
		}
	}


	// sampling Eta
	for (t = 0; t < sT; t++) {

		sWindow = 1 + ranu(1, 1) * (2 * sWindowSize - 2); // stochastic window size
		if (t < 1) sWStt = t;
		else sWStt = max(0, t-sWindow+1);
		
		vERVBV = sumc(mRVBV[sWStt:t][] .* (mRVBV[sWStt:t][] .> sJumpThres))'
				 ./ sumc(mRVBV[sWStt:t][] .> sJumpThres)';
		vVRVBV =  sumc(mRVBV[sWStt:t][] .^ 2 .* (mRVBV[sWStt:t][] .> sJumpThres))'
				  ./ sumc(mRVBV[sWStt:t][] .> sJumpThres)'
				  - vERVBV .^ 2;
		
		vERVBV = (vERVBV .< 10 ^ (-5) .|| isnan(vERVBV)) .? vERVBVall .: vERVBVall;
		vVRVBV = (vVRVBV .< 10 ^ (-10) .|| isnan(vVRVBV)) .? vVRVBVall .: vVRVBVall;
	
		valp0 = (vERVBV.^2) ./ vVRVBV + 2;
		vbet0 = vERVBV .* (valp0 - 1) ;

		mEta[t][] = rangamma(1, sD, valp0', sumc(mj[sWStt:t][].* mX[sWStt:t][]) + vbet0');
		
	}
	mEz2 = ones(sT, sD) ./ mEta;

	// variance of error for RVBV: sigma_Z 
	for (i = 0; i < sD; i++) {
		sS = sG0 + sT;
		sSS = mG0[i][i] * sG0 + sumsqrc(mRVBV[][i]  - mj[][i] .* mX[][i]);
		mSigZ2inv[i][i] = rangamma(1, 1, sS / 2, sSS / 2);
		mSigZ2[i][i] = 1 / mSigZ2inv[i][i];
	}
	
	// sampling jt
	vPr = zeros(sD, 2);
	for (i = 0; i < sD; i++) {
		for (t = 0; t < sT; t++) {
			vPr[i][0] = log(1 - mLam[t][i] * sTerm) -0.5 * (mRVBV[t][i]) ^ 2 * mSigZ2inv[i][i];
			vPr[i][1] = log(mLam[t][i] * sTerm) -0.5 * (mRVBV[t][i] - mX[t][i]) ^ 2 * mSigZ2inv[i][i];
			mj[t][i] = ranu(1, 1) <= exp(vPr[i][1]) / (exp(vPr[i][0]) + exp(vPr[i][1])) ? 1 : 0;
		}
	}


	// sampling sigmaRVBV and lambda
	sS = sH0 + sT;
	sSS = mH0 * ones(sD, 1) * sH0;
	for  (t = 0; t < sT; t++) {
		vS =  fvD(vLamAlp, vLamBeta, vLamAve, sTau, mEz2[t][]');
		vSS = fvC(vLamAlp, vLamBeta, sTau, mEz2[t][]');
		sSS += (mRVBV[t][]' - vS - diag(vSS) * mLam[t][]') .^ 2;
	}
	mSigRVBV2inv = diag(rangamma(sD, 1, sS * ones(sD, 1) / 2, sSS / 2));
	mSigRVBV2 = invert(mSigRVBV2inv);	

	mS = choleski(mSigRVBV2);
	mSS = choleski(mSigLam);

	do {
		vNod = 0
	 	~ floor((sT-1)*(range(1,sNnod)+ranu(1,sNnod))/(sNnod+2))
   	    ~ (sT-1);
	} while (prodc(diff0(vNod',1)[1:]) == 0);
	
	for (j = 0; j <= sNnod; j++) {
		sNumNod = vNod[j+1]-vNod[j]; 
		amX = amZ = amG = amW = amT = amF = new array[sNumNod];
		for (t = 0; t < sNumNod; t++) {
			amX[t] = fvD(vLamAlp, vLamBeta, vLamAve, sTau, mEz2[t][]') ~ zeros(sD, 1);
			amZ[t] = diag(fvC(vLamAlp, vLamBeta, sTau, mEz2[t][]'));
			amG[t] = mS ~ zeros(sD, sD);
			amW[t] = zeros(sD, 1) ~ 0.5 * vLamAlp .* vLamAve * sTau + vLamBeta .* mj[vNod[j]+t][]';
			amT[t] = unit(sD) - 0.5 * diag(vLamAlp) * sTau;
			amF[t] = zeros(sD, sD) ~ mSS;
		}
		vBeta = <1; 1>;
		mQ0 = amF[0];
		mKalmanOut = fKalmanFilter(mRVBV[vNod[j]:vNod[j+1]-1][], amX, vBeta, amZ, amG, amW, amT, amF, mQ0, 1)[0:sNumNod][];
		mLam[vNod[j]:vNod[j+1]-1][] = mKalmanOut[1:sNumNod][];
	}
	mLam[sT-1][] = (amT[0]*mLam[sT-2][]' + 0.5 * vLamAlp .* vLamAve * sTau + vLamBeta .* mj[sT-1][]' + choleski(mSigLam)*rann(sD, 1) )';

	mLam = (mLam .< 0) .? 0 .: mLam;
	
	// sampling LamAlpha
	mS = invert(diag(vLamAlpB0));
	vS = mS * vLamAlpA0;
	for (t = 1; t < sT; t++) {
		mSS = diag(vLamAve - mLam[t-1][]') * 0.5 * sTau;
		mS += mSS ' mSigLaminv * mSS;
		vS += mSS ' mSigLaminv * (mLam[t][]' - mLam[t-1][]' - diag(vLamBeta) * mj[t][]');
	}
	mS = invert(mS);
	vS = mS * vS;
	vLamAlp = vS + choleski(mS) * rann(sD, 1);
	
	// sampling LamBeta
	mS = invert(diag(vLamBetaB0));
	vS = mS * vLamBetaA0;
	for (t = 1; t < sT; t++) {
		mSS = diag(mj[t][]);
		mS += mSS ' mSigLaminv * mSS;
		vS += mSS ' mSigLaminv * (mLam[t][]' - mLam[t-1][]' - 0.5 * vLamAlp .* (vLamAve - mLam[t-1][]') * sTau);
	}
	mS = invert(mS);
	vS = mS * vS;
	vLamBeta = maxc((vS + choleski(mS) * rann(sD, 1))' | zeros(1, sD))';


	// sampling LamAve
	mS = invert(diag(vLamAveB0));
	vS = mS * vLamAveA0;
	mSS = 0.5 * diag(vLamAlp) * sTau;	
	for (t = 1; t < sT; t++) {
		vS += mSS ' mSigLaminv * (mLam[t][]' - mLam[t-1][]' + 0.5 * vLamAlp .* mLam[t-1][]' * sTau - vLamBeta .* mj[t][]');
	}
	mS += (sT - 1) * mSS ' mSigLaminv * mSS;
	mS = invert(mS);
	vS = mS * vS;
	vLamAve = maxc((vS + choleski(mS) * rann(sD, 1))' | zeros(1, sD))';

	
	// sampling mSigLam
	sS  = sF0 + sT - 1;
	mS = invert(mF0);
	for (t = 1; t < sT; t++) {
		vS = mLam[t][]' - (unit(sD) - 0.5 * diag(vLamAlp) * sTau) * mLam[t-1][]' 
				   - 0.5 * diag(vLamAlp) * vLamAve * sTau - diag(vLamBeta) * mj[t][]'; 
		mS += vS * vS';
	}
	mS = invert(mS);
	mSigLaminv = fMyWishart(sS, mS);
	mSigLaminv = diag(diagonal(mSigLaminv));
	mSigLam = invert(mSigLaminv);


	mOut = mj ~ mX ~ mLam ~ mEta;
	mParamOut = vMuJ ~ meanc(mEta)' ~ vLamAlp ~ vLamBeta ~ vLamAve
			 ~ diagonal(mSigZ2)' ~ diagonal(mSigLam)' ~ diagonal(mSigRVBV2)';

	return {mOut, mParamOut};
}

fvC(const vAlp, const vBeta, const sTau, const vEz2) {
	decl vOut, vA, vB;
	vA = (vBeta-vAlp) * sTau;
	vB = exp(vA) - 1;
	vOut = vB ./ vA .* vEz2;
	return vOut;
}

fvD(const vAlp, const vBeta, const vLamAve, const sTau, const vEz2)	{
	decl vA, vB;
	vA = vAlp .* vLamAve ./ (vBeta - vAlp);
	vB = fvC(vAlp, vBeta, sTau, vEz2) - vEz2;
	return vA .* vB;
}

faSetParams(const mLatent, const mParams, const sTerm, const sDim, const mG) {
	decl aP;
	aP = new array[5];
//	aP[0] = sDim;
	aP[1] = sTerm;
	aP[2] = mParams;
	aP[3] = mLatent;
	aP[4] = mG;
	return aP;
}


fvEJt(const vOme, const aP, const t) {
	decl sT, sD, sTerm, vEta, vAlp, vBeta, vLamAve, vLam, vOut, vp;
	sTerm = aP[1];
//	vEta  = aP[2][][1];
	vAlp  = aP[2][][2];
	vBeta =	aP[2][][3];
	sD = rows(vAlp);
	vLamAve = aP[2][][4];
	vLam    = aP[3][t][sD*2:sD*3-1]';
	vEta	= aP[3][t][sD*3:sD*4-1]';
	sT = rows(vLam);

	vp = (exp((vBeta-vAlp-vOme)*sTerm)-ones(sD, 1)) ./ (vBeta-vAlp-vOme);
	vOut = vp .* vLam ./ (vEta * sTerm)
		+ vAlp.*vLamAve ./ (vBeta-vAlp-vOme) .* (vp-sTerm) ./ (vEta * sTerm);
	
	return vOut;
}

fmHawkesEJ(const mOme, const mLatent, const mParams, const sTerm)
{
	decl vEta, vAlp, vBeta, sT, vLamAve, mLam;
	decl EJ, sD, p, q, mXt;
	decl mEta;
	
	sT = rows(mLatent);
	sD = columns(mOme);
//	vEta    = mParams[][1];
	vAlp    = mParams[][2];
	vBeta   = mParams[][3];
	vLamAve = mParams[][4];
	mXt = ones(sT, 1) * (vBeta - vAlp)' - mOme;
	mLam = mLatent[][sD*2:sD*3-1];
	mEta = mLatent[][sD*3:sD*4-1];
	
	p = (exp(mXt * sTerm) - ones(sT, sD)) ./ mXt;
	q = (mLam .* p + ones(sT, 1) * (vAlp .* vLamAve)'./ mXt .* (p - sTerm)) ./ (sTerm * mEta);
	
	return q;
	
}


fmLACoeffEJ(const aP) {

	decl sD, sT, sTerm, vAlp, vBeta, vLamAve, mLam, mEta, mOut, va;
	sTerm	= aP[1];
	vAlp    = aP[2][][2];
	vBeta   = aP[2][][3];
	vLamAve = aP[2][][4];
	sD = rows(vAlp);
	mLam    = aP[3][][sD*2:sD*3-1];
	mEta    = aP[3][][sD*3:sD*4-1];
	sT = rows(aP[3]);
	
	va = exp((vBeta - vAlp) * sTerm) - ones(sD, 1) - (vBeta - vAlp) * sTerm ;
	mOut = ones(sT, 1) * (va ./ ((vBeta - vAlp) .^ 2 * sTerm))';
	mOut = (mLam + ones(sT, 1) * (vAlp .* vLamAve ./ (vBeta - vAlp))') .* mOut;
	mOut = (mOut - ones(sT, 1) * (vAlp .* vLamAve ./ (vBeta - vAlp))' * sTerm * 0.5) ./ mEta;
	
	return -1 * mOut;
}

fmLAIntercEJ(const aP) {
	return aP[4];
}


