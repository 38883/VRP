extern fKalmanFilter(const mY, const amX, const vBeta, const amZ, const amG,
			const amW, const amT, const amF, const mF0, const iRnd);
extern fMyWishart(const n,const mS);
extern fKalman1(const vy, const va, const vb, const ssigu, const some, const sphi, const ssige);
			
fvSetaP(const aPfull, const ir, const ird) {
	decl aP, mLatent, mGorJ, i;
	aP = new array[5];
	for (i = 0; i < 3; i++) aP[i] = aPfull[i];
	mLatent = aPfull[3];
	mGorJ   = aPfull[4];
	aP[3]   = mLatent[ir:ird][];
	aP[4]   = mGorJ[ir:ird][];
	return aP;
}

fLinearApproxARMH(const mY, const mX,
			const fv, const fmLAInterc, const fmLACoeff,  const aPfull,
			const vMuAR1, const vPhiAR1,
			const mSigU2, const  mSigV2, const mRange, const sNnods) {
// nonlinear state space form : Yt = f(Xt, p) + ut, Xt+1 = Xt + vt
// input vY: sample of Y (one dim); fvF = f(Xt,p), fvd1F = f'(Xt,p), fvd2F = f''(Xt,p); aP = p (array of parameters in f)
// sNnods: # blocks	, sSigU2inv: Inverse of VARIANCE for ut

	decl sT, sD, vk, sNFindMode, i, j, t, k, aP, scount, sMaxAR, sAlp, mSigU, mSigAR;
	decl id, ird, ir, mYi, mXs, mXn, mXo, sPtrue, sPprop, sPtruePrev, sPpropPrev;
	decl mSigU2inv, mSigV2inv, mOut, mIntc, mCoef, mSig2AR, mSig2ARinv, vV, vM;
	decl amX, amZ, amG, amW, amT, amF, amF0, vBeta, me, mF0, amDinv, mK, amL, iRnd, va0;
	decl mMu0, vMu0, mPhi0, vPhi0, vMu1, mMu1, vPhi1, mPhi1, vPhiAR1cand, sH0, vH0, sS, vS, vSS; 
	decl vSig2V;
	
	sT = rows(mY);
	sD = columns(mY);
	
	sNFindMode = 1;
	sMaxAR = 7;
	sAlp = 2.5;

	
	mSigU2inv = diag(ones(1, sD) ./ diagonal(mSigU2));
	mSigV2inv = invert(mSigV2[0]);

	mXs = mX;
	vMu0 = 0.1 * ones(sD, 1);
	mMu0 = unit(sD);
	vPhi0 = 0.01 * ones(sD, 1);
	mPhi0 = 0.01 * unit(sD);
	sH0 = 10;
	vH0 = 0.01 * ones(sD, 1);
	
	
// Let VRP coefficients Brownian motions
	vMuAR1[0] = zeros(sD, 1);
	vPhiAR1[0] = ones(sD, 1);	

	// sampling mSigU2
	sS = sH0 + sT - 1;
	vS = sS * vH0;
	for (t = 1; t < sT; t++) {
		vSS = mX[t][]' - vMuAR1[0] - diag(vPhiAR1[0]) * mX[t-1][]'; 
		vS += vSS .^ 2;
	}
	vSig2V = ones(sD, 1) ./ rangamma(sD, 1, sS * 0.5, vS * 0.5);

	mSigV2[0] = diag(vSig2V);
	mSig2AR = mSigV2[0];
	mSigAR = choleski(mSig2AR);
	mSig2ARinv = diag(ones(sD, 1) ./ vSig2V);
	
	mSigU = diag(sqrt(diagonal(mSigU2)));	

	vk = 0 ~ floor(sT * (range(1, sNnods) + ranu(1, sNnods)) / (sNnods + 2)) ~ sT;

	
	for (i = 0; i <= sNnods; i++) {
		ir = vk[i];
		id = vk[i+1] - vk[i];
		ird = vk[i+1] - 1;

		if (!id) break;

		mYi = mY[ir:ird][];
		mXo = mX[ir:ird][];
		mXn = zeros(id, sD);

 		aP  = fvSetaP(aPfull, ir, ird);
		mIntc = fmLAInterc(aP);
		mCoef = fmLACoeff(aP);


		for (j = 0; j < sD; j++) {
			mOut = fKalman1(mYi[][j], mIntc[][j], mCoef[][j], mSigU[j][j], vMuAR1[0][j], vPhiAR1[0][j], mSigAR[j][j]);
			mXn[][j] = (mOut     .<  mRange[j][1]) .? mOut     .: mRange[j][1];
			mXn[][j] = (mXn[][j] .>  mRange[j][0]) .? mXn[][j] .: mRange[j][0];
		}

		mXs[ir:ird][] = mXn;

	}

	
	return mXs;
}


fKalman1(const vy, const va, const vb, const ssigu, const some, const sphi, const ssige) {
	decl t, sT, vnu,vF, vP, vP0, vPT, vx, vx0, vxT, Pstar, vout;
	sT = rows(vy);

	vnu = vx = vF = vP = vout = zeros(sT, 1);
	vx0 = vP0 = vxT = vPT = zeros(sT+1, 1);
	
	vx0[0] = (sphi == 1) ? 0 : some / (1 - sphi);
	vP0[0] = (sphi == 1) ? ssige^2 : ssige^2 / (1 - sphi^2);
	for (t = 0; t < sT; t++) {
		vnu[t] = vy[t] - va[t] - vb[t] * vx0[t];
		vF[t]  = vb[t]^2 * vP0[t] + ssigu^2;
		vx[t]  = vx0[t] + vb[t] * vP0[t] * vnu[t] / vF[t];
		
		vP[t]  = vP0[t] - vb[t]^2 * vP0[t]^2 / vF[t];
		vx0[t+1] = some + sphi * vx[t];
		vP0[t+1] = sphi^2 * vP[t] + ssige^2;
	}
	vxT[sT] = vx0[sT];
	vPT[sT] = vP0[sT];
	for (t = sT-1; t >= 0; t--) {
		Pstar = sphi * vP[t] / vP0[t+1];
		vxT[t] = vx[t] + Pstar * (vxT[t+1] - vx0[t+1]);
		vPT[t] = vP[t] + Pstar^2 * (vPT[t+1] - vP0[t+1]);
		vout[t] = vxT[t] + sqrt(vPT[t]) * rann(1,1);
	}
	return vout;
}

