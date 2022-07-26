double Process::matrixMSSMVirt_uu_suLsuR(double alphaS, double S, double T, 
   const double FiniteGs, const double Dminus4, int divergence, double mu) {   
	//ltini(); // for LoopTools
	setmudim(pow(mu,2));
	setlambda(divergence);   
    double U = pow(MassSq, 2) + pow(MassSq, 2) - S - T;  
	int Divergence;    // UV divergence from gauge coupling
	if(divergence == 0 || divergence == -2) {		    
	    Divergence = 0;     
	}
	else if(divergence == -1) {
		Divergence = 1;
	}   
    std::complex<double> matrix = (402.1238596594935*pow(alphaS,3)*(-1.*pow(MassSq,4) + T*U)*(0.21875*(C0i(cc0,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + C0i(cc0,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.)) + ((0.5625 + 0.09375*Dminus4)*B0i(bb0,T,0.,MassGlu*MassGlu) - 0.1875*(-1.*(MassGlu*MassGlu) - 1.*(MassSq*MassSq))*C0i(cc0,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + 0.1875*((MassGlu*MassGlu + MassSq*MassSq)*C0i(cc0,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + (-1.*(MassSq*MassSq) + T)*C0i(cc0,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.)) + 0.020833333333333332*(2.*(-1.*(MassSq*MassSq) + T)*(C0i(cc0,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq) + C0i(cc1,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq)) + MassGlu*MassGlu*(C0i(cc0,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq))) - 1.*(0.09375*Dminus4*(MassGlu*MassGlu - 1.*(MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu - 1.*(MassSq*MassSq) - 1.*T))*C0i(cc1,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.010416666666666666*(-2.*B0i(bb0,T,0.,MassSq*MassSq) + 2.*(-3.*(MassSq*MassSq) + T)*C0i(cc2,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu + MassSq*MassSq)*(C0i(cc1,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + C0i(cc2,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq)) + 0.1875*(-1.*(MassSq*MassSq) + T)*C0i(cc2,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.))/(-1.*(MassGlu*MassGlu) + T) - 0.14583333333333334*(2.*(-1.*(MassSq*MassSq) + U)*D0i(dd0,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd00,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-2.*(MassSq*MassSq) + U)*D0i(dd1,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) + 0.07291666666666667*(C0i(cc1,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + C0i(cc1,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + C0i(cc2,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) - 2.*(MassSq*MassSq)*D0i(dd11,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + U)*(D0i(dd12,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd2,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd23,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq))) + 0.14583333333333334*((MassSq*MassSq - 1.*U)*D0i(dd2,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + MassSq*MassSq*D0i(dd3,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) - 0.07291666666666667*((2.*(MassSq*MassSq) + 2.*U)*D0i(dd1,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (3.*(MassSq*MassSq) + U)*D0i(dd13,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassSq*MassSq + 3.*U)*D0i(dd3,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassSq*MassSq + U)*D0i(dd33,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq))))/(-1.*(MassGlu*MassGlu) + T) + (402.1238596594935*pow(alphaS,3)*(-1.*pow(MassSq,4) + T*U)*(0.21875*C0i(cc0,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq) + 0.14583333333333334*C0i(cc0,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + ((0.5625 + 0.09375*Dminus4)*B0i(bb0,U,0.,MassGlu*MassGlu) - 0.1875*(-1.*(MassGlu*MassGlu) - 1.*(MassSq*MassSq))*C0i(cc0,U,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + 0.1875*((MassGlu*MassGlu + MassSq*MassSq)*C0i(cc0,U,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + (-1.*(MassSq*MassSq) + U)*C0i(cc0,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.)) + 0.020833333333333332*(2.*(-1.*(MassSq*MassSq) + U)*(C0i(cc0,0.,U,MassSq*MassSq,0.,0.,MassSq*MassSq) + C0i(cc1,0.,U,MassSq*MassSq,0.,0.,MassSq*MassSq)) + MassGlu*MassGlu*(C0i(cc0,MassSq*MassSq,U,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,U,0.,MassGlu*MassGlu,0.,MassSq*MassSq))) - 1.*(0.09375*Dminus4*(MassGlu*MassGlu - 1.*(MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu - 1.*(MassSq*MassSq) - 1.*U))*C0i(cc1,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.010416666666666666*(-2.*B0i(bb0,U,0.,MassSq*MassSq) + 2.*(-3.*(MassSq*MassSq) + U)*C0i(cc2,0.,U,MassSq*MassSq,0.,0.,MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu + MassSq*MassSq)*(C0i(cc1,U,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + C0i(cc2,U,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq)) + 0.1875*(-1.*(MassSq*MassSq) + U)*C0i(cc2,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.))/(-1.*(MassGlu*MassGlu) + U) + (0.08333333333333333 + 0.041666666666666664*Dminus4)*(C0i(cc0,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) - 1.*D0i(dd00,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*U*D0i(dd11,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + (0.041666666666666664 + 0.020833333333333332*Dminus4)*(C0i(cc1,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) + C0i(cc2,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) - 1.*(-1.*(MassSq*MassSq) + U)*D0i(dd12,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*(-1.*(MassSq*MassSq) + U)*D0i(dd13,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) - 0.14583333333333334*((-1.*(MassSq*MassSq) + T)*D0i(dd0,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd00,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + T)*D0i(dd2,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) + 0.041666666666666664*(-1.*S*(D0i(dd0,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) + D0i(dd1,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) - 1.*S*D0i(dd2,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*S*D0i(dd3,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + 0.07291666666666667*(C0i(cc2,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq) + (3.*(MassGlu*MassGlu) - 2.*(MassSq*MassSq) + 2.*S - 1.*U)*D0i(dd0,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassGlu*MassGlu - 3.*(MassSq*MassSq) + 2.*S - 2.*U)*D0i(dd1,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-3.*(MassSq*MassSq) + 2.*S + 5.*U)*D0i(dd1,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + U)*D0i(dd11,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-3.*(MassSq*MassSq) + S - 1.*U)*D0i(dd13,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-2.*(MassSq*MassSq) + S + 2.*U)*D0i(dd13,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassSq*MassSq + 2.*S - 1.*U)*D0i(dd2,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassGlu*MassGlu + 2.*(MassSq*MassSq) - 3.*T - 4.*U)*D0i(dd3,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-2.*(MassSq*MassSq) + S)*D0i(dd33,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) - 0.07291666666666667*(C0i(cc1,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassSq*MassSq + U)*D0i(dd11,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + U)*(D0i(dd12,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd23,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) + (-5.*(MassSq*MassSq) + 3.*T)*D0i(dd3,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + T)*D0i(dd33,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq))))/(-1.*(MassGlu*MassGlu) + U) - (25.132741228718345*alphaS*(-1.*pow(MassSq,4) + T*U)*(-16.*(alphaS*alphaS)*(((0.1875*Dminus4*(MassGlu*MassGlu - 1.*T) + 0.375*(3.*(MassGlu*MassGlu) - 1.*T))*B0i(bb0,T,0.,MassGlu*MassGlu) + (MassGlu*MassGlu + T)*(-1.*(0.375 + 0.1875*Dminus4)*B0i(bb1,T,0.,MassGlu*MassGlu) + 0.0625*(10.*B0i(bb1,T,0.,MassSq*MassSq) + 2.*B0i(bb1,T,MassTop*MassTop,MassSq*MassSq))))/pow(-1.*(MassGlu*MassGlu) + T,2) + ((0.375 + 0.09375*Dminus4)*B0i(bb0,T,0.,MassGlu*MassGlu) + 0.1875*(-1.*(MassSq*MassSq) + T)*C0i(cc0,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.020833333333333332*(MassGlu*MassGlu)*(C0i(cc0,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) - 1.*(0.09375*Dminus4*(MassGlu*MassGlu - 1.*(MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu - 1.*(MassSq*MassSq) - 1.*T))*C0i(cc1,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.1875*(-1.*(MassSq*MassSq) + T)*C0i(cc2,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.))/(-1.*(MassGlu*MassGlu) + T) + 0.041666666666666664*(MassSq*MassSq)*D0i(dd11,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (0.08333333333333333 + 0.041666666666666664*Dminus4)*(C0i(cc0,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) - 1.*D0i(dd00,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*T*D0i(dd11,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + (0.041666666666666664 + 0.020833333333333332*Dminus4)*(C0i(cc1,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) + C0i(cc2,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) - 1.*(-1.*(MassSq*MassSq) + T)*D0i(dd12,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*(-1.*(MassSq*MassSq) + T)*D0i(dd13,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + 0.041666666666666664*(D0i(dd00,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) - 1.*S*(D0i(dd0,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) + D0i(dd1,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) - 1.*S*D0i(dd2,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*S*D0i(dd3,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + 0.020833333333333332*(C0i(cc0,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + C0i(cc2,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + T + 2.*U)*D0i(dd0,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + S + 2.*T + 3.*U)*D0i(dd1,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (3.*(MassSq*MassSq) + U)*D0i(dd12,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (3.*(MassSq*MassSq) + U)*D0i(dd13,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + MassSq*MassSq + T + 3.*U)*D0i(dd2,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassSq*MassSq + U)*D0i(dd22,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (2.*(MassSq*MassSq) + 2.*U)*D0i(dd23,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + MassSq*MassSq + T + 3.*U)*D0i(dd3,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassSq*MassSq + U)*D0i(dd33,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq))) - 12.566370614359172*alphaS*((-1.*(-1.*MassGlu*(-0.954929658551372*alphaS*MassGlu*Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*MassGlu*Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + 0.15915494309189535*alphaS*MassGlu*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))) + MassGlu*(-0.6366197723675814*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 0.954929658551372*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.15915494309189535*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))))) + T*(-0.6366197723675814*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 0.954929658551372*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.15915494309189535*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))))))/pow(-1.*(MassGlu*MassGlu) + T,2) + (0.01688686394038963*(12.566370614359172*alphaS*FiniteGs + 29.608813203268074*(4.*(-0.1193662073189215*alphaS*Divergence - 0.07957747154594767*alphaS*FiniteGs*(log((MassGlu*MassGlu)/(mu*mu)) + log((MassSq*MassSq)/(mu*mu)) + 0.3333333333333333*log((MassTop*MassTop)/(mu*mu)))) - 2.5464790894703255*alphaS*((0.16666666666666666 + 0.08333333333333333*Dminus4)*(Re(B0i(bb0,0.,0.,0.)) + Re(B0i(bb1,0.,0.,0.))) - 0.16666666666666666*Re(B0i(bb1,0.,MassGlu*MassGlu,MassSq*MassSq))) - 1.2732395447351628*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 1.909859317102744*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 2.5464790894703255*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.3183098861837907*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))) + 2.5464790894703255*alphaS*(-0.3333333333333333*Re(B0i(bb0,MassSq*MassSq,0.,MassGlu*MassGlu)) + 0.25*Re(B0i(bb0,MassSq*MassSq,0.,MassSq*MassSq)) + 0.16666666666666666*Re(B0i(bb1,MassSq*MassSq,0.,MassSq*MassSq)) - 0.3333333333333333*Re(B0i(bb1,MassSq*MassSq,MassGlu*MassGlu,0.)) - 0.3333333333333333*(MassSq*MassSq)*Re(B0i(dbb0,MassSq*MassSq,0.,MassGlu*MassGlu)) + 0.3333333333333333*(MassSq*MassSq)*Re(B0i(dbb0,MassSq*MassSq,0.,MassSq*MassSq)) + MassSq*MassSq*(0.16666666666666666*Re(B0i(dbb1,MassSq*MassSq,0.,MassSq*MassSq)) - 0.3333333333333333*Re(B0i(dbb1,MassSq*MassSq,MassGlu*MassGlu,0.)))))))/(-1.*(MassGlu*MassGlu) + T))))/(-1.*(MassGlu*MassGlu) + T) + (25.132741228718345*alphaS*(-1.*pow(MassSq,4) + T*U)*(16.*(alphaS*alphaS)*(((0.1875*Dminus4*(MassGlu*MassGlu - 1.*U) + 0.375*(3.*(MassGlu*MassGlu) - 1.*U))*B0i(bb0,U,0.,MassGlu*MassGlu) + (MassGlu*MassGlu + U)*(-1.*(0.375 + 0.1875*Dminus4)*B0i(bb1,U,0.,MassGlu*MassGlu) + 0.0625*(10.*B0i(bb1,U,0.,MassSq*MassSq) + 2.*B0i(bb1,U,MassTop*MassTop,MassSq*MassSq))))/pow(-1.*(MassGlu*MassGlu) + U,2) + ((0.375 + 0.09375*Dminus4)*B0i(bb0,U,0.,MassGlu*MassGlu) + 0.1875*(-1.*(MassSq*MassSq) + U)*C0i(cc0,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.020833333333333332*(MassGlu*MassGlu)*(C0i(cc0,MassSq*MassSq,U,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,U,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) - 1.*(0.09375*Dminus4*(MassGlu*MassGlu - 1.*(MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu - 1.*(MassSq*MassSq) - 1.*U))*C0i(cc1,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.1875*(-1.*(MassSq*MassSq) + U)*C0i(cc2,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.))/(-1.*(MassGlu*MassGlu) + U) + 0.041666666666666664*D0i(dd00,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + 0.041666666666666664*(MassSq*MassSq)*D0i(dd11,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + 0.020833333333333332*(C0i(cc0,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + C0i(cc2,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + 2.*T + U)*D0i(dd0,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + S + 3.*T + 2.*U)*D0i(dd1,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (3.*(MassSq*MassSq) + T)*D0i(dd12,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (3.*(MassSq*MassSq) + T)*D0i(dd13,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + MassSq*MassSq + 3.*T + U)*D0i(dd2,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassSq*MassSq + T)*D0i(dd22,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (2.*(MassSq*MassSq) + 2.*T)*D0i(dd23,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + MassSq*MassSq + 3.*T + U)*D0i(dd3,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassSq*MassSq + T)*D0i(dd33,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq))) + 12.566370614359172*alphaS*((-1.*(-1.*MassGlu*(-0.954929658551372*alphaS*MassGlu*Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*MassGlu*Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + 0.15915494309189535*alphaS*MassGlu*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))) + MassGlu*(-0.6366197723675814*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 0.954929658551372*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.15915494309189535*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))))) + U*(-0.6366197723675814*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 0.954929658551372*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.15915494309189535*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))))))/pow(-1.*(MassGlu*MassGlu) + U,2) + (0.01688686394038963*(12.566370614359172*alphaS*FiniteGs + 29.608813203268074*(4.*(-0.1193662073189215*alphaS*Divergence - 0.07957747154594767*alphaS*FiniteGs*(log((MassGlu*MassGlu)/(mu*mu)) + log((MassSq*MassSq)/(mu*mu)) + 0.3333333333333333*log((MassTop*MassTop)/(mu*mu)))) - 2.5464790894703255*alphaS*((0.16666666666666666 + 0.08333333333333333*Dminus4)*(Re(B0i(bb0,0.,0.,0.)) + Re(B0i(bb1,0.,0.,0.))) - 0.16666666666666666*Re(B0i(bb1,0.,MassGlu*MassGlu,MassSq*MassSq))) - 1.2732395447351628*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 1.909859317102744*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 2.5464790894703255*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.3183098861837907*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))) + 2.5464790894703255*alphaS*(-0.3333333333333333*Re(B0i(bb0,MassSq*MassSq,0.,MassGlu*MassGlu)) + 0.25*Re(B0i(bb0,MassSq*MassSq,0.,MassSq*MassSq)) + 0.16666666666666666*Re(B0i(bb1,MassSq*MassSq,0.,MassSq*MassSq)) - 0.3333333333333333*Re(B0i(bb1,MassSq*MassSq,MassGlu*MassGlu,0.)) - 0.3333333333333333*(MassSq*MassSq)*Re(B0i(dbb0,MassSq*MassSq,0.,MassGlu*MassGlu)) + 0.3333333333333333*(MassSq*MassSq)*Re(B0i(dbb0,MassSq*MassSq,0.,MassSq*MassSq)) + MassSq*MassSq*(0.16666666666666666*Re(B0i(dbb1,MassSq*MassSq,0.,MassSq*MassSq)) - 0.3333333333333333*Re(B0i(dbb1,MassSq*MassSq,MassGlu*MassGlu,0.)))))))/(-1.*(MassGlu*MassGlu) + U))))/(-1.*(MassGlu*MassGlu) + U);
//    ltexi();    // error message of LoopTools 
    double matrixReal = matrix.real();
    return matrixReal;
}
