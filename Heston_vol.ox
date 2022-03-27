
fvG(const vGamt, const vKappa, const vTheta, const sTau, const sTerm, const vSigt2) {

	decl vd0, ve0, vd, sD = rows(vGamt);

	vd0 = (ones(sD, 1) - exp(-(vKappa + vGamt) * sTerm)) ./ ((vKappa + vGamt) * sTerm);
	ve0 = vKappa .* vTheta ./ (vKappa + vGamt);

	vd = vd0 .* vSigt2 + ve0 .* (ones(sD, 1) - vd0);
	return vd;
}

fvGt(const vGamt, const aP, const t) {

	decl vKappa, vTheta, sTerm;
	decl vd0, ve0, vd, sD, vSigt2;
	sTerm  = aP[1];
	vKappa = aP[2][][0];
	vTheta = aP[2][][1];
	sD = rows(vKappa);
	vSigt2  = aP[3][t][0:sD-1]';
	
	vd0 = (ones(sD, 1) - exp(-(vKappa + vGamt) * sTerm)) ./ ((vKappa + vGamt) * sTerm);
	ve0 = vKappa .* vTheta ./ (vKappa + vGamt);

	vd = vd0 .* vSigt2 + ve0 .* (ones(sD, 1) - vd0);
	return vd;
	
}

fvSetaP_G(const mGam, const msigt, const vKappa, const vTheta, const sTerm, const sDim, const mJ) {
	decl aP = new array[5];
	aP[0] = sDim;
	aP[1] = sTerm;
	aP[2] = vKappa ~ vTheta;
	aP[3] = msigt ~ mGam;
	aP[4] = mJ;
	return aP;
}

faGetG(const aP) {
	decl dim, vsigt, sKappa, sTheta, sTerm;
	dim    = aP[0]; 
	sTerm  = aP[1];
	sKappa = aP[2][dim][0];
	sTheta = aP[2][dim][1];
	vsigt  = aP[3][][dim];
	return {vsigt, sKappa, sTheta, sTerm};
}

fv_G(const vGam, const aP) {
	decl vsigt, sKappa, sTheta, sTerm, vOut, vp, vpdot, v1;
	[vsigt, sKappa, sTheta, sTerm] = faGetG(aP);
	v1 = ones(rows(vsigt), 1);	
	vp    = (v1- exp(-(sKappa*v1 + vGam) * sTerm)) ./ (sTerm * (sKappa*v1 + vGam));
	vOut = vp .* vsigt + sKappa * sTheta * v1 ./ (sKappa * v1 + vGam)	.* (v1 - vp);
	return vOut;
}

fmLACoeffG(const aP) {
// this function returns coefficient of linear apploximated G
	decl i, sT, sD, mOut, mLow, mX;
	decl msigt, vKappa, vTheta, sTerm;
	decl mAveSigt, sWindow = 21;
	
	sD = rows(aP[2]);
	sT = rows(aP[3]);
	mAveSigt = zeros(sT, sD);
	sTerm  = aP[1];
	vKappa = aP[2][][0];
	vTheta = aP[2][][1];
	msigt  = aP[3][][0:sD-1];
	mAveSigt = msigt;
	
	mOut = mAveSigt - ones(sT, 1) * vTheta';
	mOut = ones(sT, 1) * ((vKappa * sTerm + exp(-vKappa * sTerm) - ones(sD, 1))
			./ (vKappa .^ 2 * sTerm))' .* mOut
			+ ones(sT, 1) * vTheta' * 0.5 * sTerm;
	
	return -1 * mOut;
}

fmLAIntercG(const aP){
// this function returns intercept of linear approximated G
	decl t, mAveV, sT, sD, sWindow = 21;
	sT = rows(aP[4]);
	sD = columns(aP[4]);
	mAveV = zeros(sT, sD);
	mAveV = aP[4];
	return mAveV;
}

