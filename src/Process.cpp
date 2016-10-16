#include "Process.hpp"

Process::Process(std::string processID, boost::property_tree::ptree pt) {
   
   MassTop = pt.get<double>("masses.top");
   MassGlu = pt.get<double>("masses.gluino");
   MasssigmaO = pt.get<double>("masses.pseudoscalar_sgluon");
   MassphiO  = sqrt( pow(MasssigmaO,2) + 4.0 * pow(MassGlu, 2) );
   MassSuL = pt.get<double>("masses.suL");
   MassSuR = pt.get<double>("masses.suR");
   MassSdL = pt.get<double>("masses.sdL");
   MassSdR = pt.get<double>("masses.sdR");
   MassSsL = pt.get<double>("masses.ssL");
   MassSsR = pt.get<double>("masses.ssR");
   MassScL = pt.get<double>("masses.scL");
   MassScR = pt.get<double>("masses.scR");
   MassSbL = pt.get<double>("masses.sbL");
   MassSbR = pt.get<double>("masses.sbR");
   MassStL = pt.get<double>("masses.stL");
   MassStR = pt.get<double>("masses.stR");
   mu_r = pt.get<double>("collider setup.mu_r");
   mu_f = pt.get<double>("collider setup.mu_f");
   dS = pt.get<double>("technical parameters.dS");
   pdf = LHAPDF::mkPDF( pt.get<std::string>("collider setup.pdf") , 0);
   
   // @todo remove 
   MassSq = MassSuL;
   
    /* squark production, MRSSM */
	if(processID == "MRSSM-uu_suLsuR") {
      // tree-level is the same as in the MSSM
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuR;
      //matrixelementVirt = &Process::f;
      matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2;
      f2 = 2;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,ud_suLsdR") {
      //matrixelementTree = &Process::matrixMRSSMTree_ud_suLsdR;
      //matrixelementVirt = &matrixMRSSMVirt_ud_suLsdR;
      f1 = 2.;
      f2 = 1.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
    /* squark production, MSSM */
   else if(processID == "MSSM,uu_suLsuR") {
      //matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      //matrixelementVirt = &Process::matrixMSSMVirt_uu_suLsuR;
      //matrixelementReal_SC = &Process::matrixMRSSMSoft_uu_suLsuRg;
      matrixelementReal_SC = &Process::g;
      m1 = MassSuL;
      m2 = MassSuR;
      f1 = 2.;
      f2 = 2.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ud_suLsdR") {
      //matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      //matrixelementVirt = &matrixMSSMVirt_ud_suLsdR;
      f1 = 2.;
      f2 = 1.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
//    /* squark anti-squark production, MRSSM */
   else if(processID == "MRSSM,uubar_suLsuLdagger") {
      //matrixelementTree = &Process::matrixMRSSMTree_uubar_suLsuLdagger;
      //matrixelementVirt = &matrixMRSSMVirt_uubar_suLsuLdagger;
      f1 = 2.;
      f2 = -2.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,GG_suLsuLdagger") {
      //matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger;
      //matrixelementVirt = &matrixMRSSMVirt_GG_suLsuLdagger;
      f1 = 0.;
      f2 = 0.;
      k = 2.*2*8*8;
      h = 1.;
   }
   else if( processID == "sgluons-qqbar_OO" ) {
      matrixelementTree = &Process::matrixSgluonTree_qqbar_OO;
      matrixelementVirt = &Process::f;
      m1 =  MasssigmaO;
      m2 =  MasssigmaO;
      k = 1;
      f1 = 69;
      f2 = 69;
      h=1;
   }
   else if( processID == "sgluons-gg_OO" ) {
      matrixelementTree = &Process::matrixSgluonTree_gg_OO;
      matrixelementVirt = &Process::f;
      m1 = m2 = MasssigmaO;
      k = 1;
      f1 = 0;
      f2 = 0;
      h=1;
   }
   else {
      std::cout << "Process not implemented.\n";
   }
}

double Process::f(double S, double T, double x, double y, int z) {
   return 0.;
}

double Process::g(double S, double T) {
   return 0.;
}

double Process::matrixMSSMTree_uu_suLsuR( double s ) {
   double MGl2 = pow(MassGlu, 2);
   double a = pdf->alphasQ( mu_r );
   return (-4.*pow(a, 2)*pi*(2.*sqrt(s*(-4*pow(m1,2) + s)) + (2.*pow(m1,2) - 2.*MGl2 - s)*
        log((4.*MGl2 + pow(1. + sqrt(1. - (4.*pow(m1,2))/s),2.)*s)/(4.*MGl2 + pow(-1. + sqrt(1. - (4.*pow(m1,2))/s),2)*s))))/(9.*pow(s,2));
}

double Process::matrixMRSSMTree_uubar_suLsuLdagger(double alphaS, double T, double U, double S) {

   double msquaredtree =  (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/(S*S) + 355.3057584392169*(alphaS*alphaS)*pow(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U) - 236.8705056261446*(alphaS*alphaS)*(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T))*(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T))*(-1.*pow(MassSq,4) + T*U) + 355.3057584392169*(alphaS*alphaS)*pow(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U);

   return msquaredtree;
}

double Process::matrixMSSMTree_ud_suLsdR(double alphaS, double T, double U, double S) {
   return (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(MassGlu*MassGlu - 1.*T,2); /*  ONLY HAVE OF THE VALUE WHICH IS COMMENTED: left-right + right-left or 2 of left-right*/
                   // + 2.*(315.82734083485946*(alphaS*alphaS)*(MassGlu*MassGlu)*S)/pow(MassGlu*MassGlu - 1.*T,2); /* left-left + right-right or 2 of left-left */
}

double Process::matrixMRSSMTree_ud_suLsdR(double alphaS, double T, double U, double S) {
   return (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2);
}

double Process::matrixMRSSMTree_GG_suLsuLdagger(double alphaS, double T, double U, double S) {
   return  -157.91367041742973*(alphaS*alphaS)*((96.*(pow(MassSq,4) + T*U - 1.*(MassSq*MassSq)*(T + U)))/(S*S) + MassSq*MassSq*(37.333333333333336/(-1.*(MassSq*MassSq) + U) - 85.33333333333333*(T/pow(-1.*(MassSq*MassSq) + T,2) + U/pow(-1.*(MassSq*MassSq) + U,2))) + (37.333333333333336*(MassSq*MassSq) + (5.333333333333333*(9.*pow(MassSq,4) - 3.*(MassSq*MassSq)*(2.*(MassSq*MassSq) + S) + (S + T)*(S + U)))/(-1.*(MassSq*MassSq) + U))/(-1.*(MassSq*MassSq) + T) - (48.*((5.*pow(MassSq,4) + MassSq*MassSq*U - 1.*T*(-1.*(MassSq*MassSq) + U) - 1.*(MassSq*MassSq)*(3.*S + 4.*T + 2.*U))/(-1.*(MassSq*MassSq) + U) + (5.*pow(MassSq,4) + MassSq*MassSq*U - 1.*T*(-1.*(MassSq*MassSq) + U) - 1.*(MassSq*MassSq)*(3.*S + 2.*T + 4.*U))/(-1.*(MassSq*MassSq) + T)))/S);
}

inline double Process::matrixSgluonTree_qqbar_OO( double s ) {
   double a = pdf->alphasQ( mu_r );
   return 2. * pi * pow(a, 2)/(9. * s) * pow(1. - 4. * m1*m1/s, 3./2. );
}

inline double Process::matrixSgluonTree_gg_OO( double s ) {
   double b = sqrt( 1. - 4. * m1*m1/s );
   double a = pdf->alphasQ( mu_r );
   return 3. * pi * pow(a, 2)/(32. * s) 
      * (27.*b - 17*pow(b, 3) + 6.*(-3. + 2. * pow(b, 2) + pow(b, 4) ) * atanh(b) ) ;
}
/*
 
 Virtual
 
 */

double Process::matrixMSSMVirt_uu_suLsuR(double S, double T, 
   const double FiniteGs, const double Dminus4, int divergence) {
   
	ltini();
	setmudim(pow(mu_r,2));
	setlambda(divergence);
   
   double alphaS = pdf->alphasQ( mu_r );
   double mu = mu_r;
   double U = pow(m1, 2) + pow(m2, 2) - S - T;
   
	int Divergence;    // UV divergence from gauge coupling
	if(divergence == 0 || divergence == -2) {		    
	    Divergence = 0;     
	}
	else if(divergence == -1) {
		Divergence = 1;
	}
   
    std::complex<double> matrix = (402.1238596594935*pow(alphaS,3)*(-1.*pow(MassSq,4) + T*U)*(0.21875*(C0i(cc0,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + C0i(cc0,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.)) + ((0.5625 + 0.09375*Dminus4)*B0i(bb0,T,0.,MassGlu*MassGlu) - 0.1875*(-1.*(MassGlu*MassGlu) - 1.*(MassSq*MassSq))*C0i(cc0,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + 0.1875*((MassGlu*MassGlu + MassSq*MassSq)*C0i(cc0,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + (-1.*(MassSq*MassSq) + T)*C0i(cc0,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.)) + 0.020833333333333332*(2.*(-1.*(MassSq*MassSq) + T)*(C0i(cc0,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq) + C0i(cc1,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq)) + MassGlu*MassGlu*(C0i(cc0,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq))) - 1.*(0.09375*Dminus4*(MassGlu*MassGlu - 1.*(MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu - 1.*(MassSq*MassSq) - 1.*T))*C0i(cc1,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.010416666666666666*(-2.*B0i(bb0,T,0.,MassSq*MassSq) + 2.*(-3.*(MassSq*MassSq) + T)*C0i(cc2,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu + MassSq*MassSq)*(C0i(cc1,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + C0i(cc2,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq)) + 0.1875*(-1.*(MassSq*MassSq) + T)*C0i(cc2,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.))/(-1.*(MassGlu*MassGlu) + T) - 0.14583333333333334*(2.*(-1.*(MassSq*MassSq) + U)*D0i(dd0,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd00,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-2.*(MassSq*MassSq) + U)*D0i(dd1,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) + 0.07291666666666667*(C0i(cc1,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + C0i(cc1,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + C0i(cc2,T,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) - 2.*(MassSq*MassSq)*D0i(dd11,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + U)*(D0i(dd12,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd2,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd23,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq))) + 0.14583333333333334*((MassSq*MassSq - 1.*U)*D0i(dd2,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + MassSq*MassSq*D0i(dd3,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) - 0.07291666666666667*((2.*(MassSq*MassSq) + 2.*U)*D0i(dd1,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (3.*(MassSq*MassSq) + U)*D0i(dd13,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassSq*MassSq + 3.*U)*D0i(dd3,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassSq*MassSq + U)*D0i(dd33,T,MassSq*MassSq,U,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq))))/(-1.*(MassGlu*MassGlu) + T) + (402.1238596594935*pow(alphaS,3)*(-1.*pow(MassSq,4) + T*U)*(0.21875*C0i(cc0,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq) + 0.14583333333333334*C0i(cc0,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + ((0.5625 + 0.09375*Dminus4)*B0i(bb0,U,0.,MassGlu*MassGlu) - 0.1875*(-1.*(MassGlu*MassGlu) - 1.*(MassSq*MassSq))*C0i(cc0,U,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + 0.1875*((MassGlu*MassGlu + MassSq*MassSq)*C0i(cc0,U,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + (-1.*(MassSq*MassSq) + U)*C0i(cc0,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.)) + 0.020833333333333332*(2.*(-1.*(MassSq*MassSq) + U)*(C0i(cc0,0.,U,MassSq*MassSq,0.,0.,MassSq*MassSq) + C0i(cc1,0.,U,MassSq*MassSq,0.,0.,MassSq*MassSq)) + MassGlu*MassGlu*(C0i(cc0,MassSq*MassSq,U,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,U,0.,MassGlu*MassGlu,0.,MassSq*MassSq))) - 1.*(0.09375*Dminus4*(MassGlu*MassGlu - 1.*(MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu - 1.*(MassSq*MassSq) - 1.*U))*C0i(cc1,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.010416666666666666*(-2.*B0i(bb0,U,0.,MassSq*MassSq) + 2.*(-3.*(MassSq*MassSq) + U)*C0i(cc2,0.,U,MassSq*MassSq,0.,0.,MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu + MassSq*MassSq)*(C0i(cc1,U,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq) + C0i(cc2,U,0.,MassSq*MassSq,0.,MassGlu*MassGlu,MassSq*MassSq)) + 0.1875*(-1.*(MassSq*MassSq) + U)*C0i(cc2,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.))/(-1.*(MassGlu*MassGlu) + U) + (0.08333333333333333 + 0.041666666666666664*Dminus4)*(C0i(cc0,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) - 1.*D0i(dd00,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*U*D0i(dd11,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + (0.041666666666666664 + 0.020833333333333332*Dminus4)*(C0i(cc1,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) + C0i(cc2,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) - 1.*(-1.*(MassSq*MassSq) + U)*D0i(dd12,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*(-1.*(MassSq*MassSq) + U)*D0i(dd13,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) - 0.14583333333333334*((-1.*(MassSq*MassSq) + T)*D0i(dd0,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd00,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + T)*D0i(dd2,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) + 0.041666666666666664*(-1.*S*(D0i(dd0,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) + D0i(dd1,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) - 1.*S*D0i(dd2,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*S*D0i(dd3,U,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + 0.07291666666666667*(C0i(cc2,0.,T,MassSq*MassSq,0.,0.,MassSq*MassSq) + (3.*(MassGlu*MassGlu) - 2.*(MassSq*MassSq) + 2.*S - 1.*U)*D0i(dd0,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassGlu*MassGlu - 3.*(MassSq*MassSq) + 2.*S - 2.*U)*D0i(dd1,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-3.*(MassSq*MassSq) + 2.*S + 5.*U)*D0i(dd1,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + U)*D0i(dd11,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-3.*(MassSq*MassSq) + S - 1.*U)*D0i(dd13,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-2.*(MassSq*MassSq) + S + 2.*U)*D0i(dd13,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassSq*MassSq + 2.*S - 1.*U)*D0i(dd2,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassGlu*MassGlu + 2.*(MassSq*MassSq) - 3.*T - 4.*U)*D0i(dd3,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-2.*(MassSq*MassSq) + S)*D0i(dd33,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) - 0.07291666666666667*(C0i(cc1,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (MassSq*MassSq + U)*D0i(dd11,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + U)*(D0i(dd12,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + D0i(dd23,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) + (-5.*(MassSq*MassSq) + 3.*T)*D0i(dd3,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + (-1.*(MassSq*MassSq) + T)*D0i(dd33,U,MassSq*MassSq,T,MassSq*MassSq,0.,0.,0.,MassGlu*MassGlu,0.,MassSq*MassSq))))/(-1.*(MassGlu*MassGlu) + U) - (25.132741228718345*alphaS*(-1.*pow(MassSq,4) + T*U)*(-16.*(alphaS*alphaS)*(((0.1875*Dminus4*(MassGlu*MassGlu - 1.*T) + 0.375*(3.*(MassGlu*MassGlu) - 1.*T))*B0i(bb0,T,0.,MassGlu*MassGlu) + (MassGlu*MassGlu + T)*(-1.*(0.375 + 0.1875*Dminus4)*B0i(bb1,T,0.,MassGlu*MassGlu) + 0.0625*(10.*B0i(bb1,T,0.,MassSq*MassSq) + 2.*B0i(bb1,T,MassTop*MassTop,MassSq*MassSq))))/pow(-1.*(MassGlu*MassGlu) + T,2) + ((0.375 + 0.09375*Dminus4)*B0i(bb0,T,0.,MassGlu*MassGlu) + 0.1875*(-1.*(MassSq*MassSq) + T)*C0i(cc0,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.020833333333333332*(MassGlu*MassGlu)*(C0i(cc0,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,T,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) - 1.*(0.09375*Dminus4*(MassGlu*MassGlu - 1.*(MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu - 1.*(MassSq*MassSq) - 1.*T))*C0i(cc1,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.1875*(-1.*(MassSq*MassSq) + T)*C0i(cc2,T,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.))/(-1.*(MassGlu*MassGlu) + T) + 0.041666666666666664*(MassSq*MassSq)*D0i(dd11,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (0.08333333333333333 + 0.041666666666666664*Dminus4)*(C0i(cc0,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) - 1.*D0i(dd00,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*T*D0i(dd11,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + (0.041666666666666664 + 0.020833333333333332*Dminus4)*(C0i(cc1,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) + C0i(cc2,MassSq*MassSq,S,MassSq*MassSq,MassGlu*MassGlu,0.,0.) - 1.*(-1.*(MassSq*MassSq) + T)*D0i(dd12,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*(-1.*(MassSq*MassSq) + T)*D0i(dd13,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + 0.041666666666666664*(D0i(dd00,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) - 1.*S*(D0i(dd0,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) + D0i(dd1,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) - 1.*S*D0i(dd2,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.) - 1.*S*D0i(dd3,T,MassSq*MassSq,S,0.,0.,MassSq*MassSq,0.,MassGlu*MassGlu,0.,0.)) + 0.020833333333333332*(C0i(cc0,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + C0i(cc2,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + T + 2.*U)*D0i(dd0,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + S + 2.*T + 3.*U)*D0i(dd1,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (3.*(MassSq*MassSq) + U)*D0i(dd12,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (3.*(MassSq*MassSq) + U)*D0i(dd13,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + MassSq*MassSq + T + 3.*U)*D0i(dd2,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassSq*MassSq + U)*D0i(dd22,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (2.*(MassSq*MassSq) + 2.*U)*D0i(dd23,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + MassSq*MassSq + T + 3.*U)*D0i(dd3,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassSq*MassSq + U)*D0i(dd33,T,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq))) - 12.566370614359172*alphaS*((-1.*(-1.*MassGlu*(-0.954929658551372*alphaS*MassGlu*Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*MassGlu*Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + 0.15915494309189535*alphaS*MassGlu*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))) + MassGlu*(-0.6366197723675814*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 0.954929658551372*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.15915494309189535*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))))) + T*(-0.6366197723675814*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 0.954929658551372*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.15915494309189535*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))))))/pow(-1.*(MassGlu*MassGlu) + T,2) + (0.01688686394038963*(12.566370614359172*alphaS*FiniteGs + 29.608813203268074*(4.*(-0.1193662073189215*alphaS*Divergence - 0.07957747154594767*alphaS*FiniteGs*(log((MassGlu*MassGlu)/(mu*mu)) + log((MassSq*MassSq)/(mu*mu)) + 0.3333333333333333*log((MassTop*MassTop)/(mu*mu)))) - 2.5464790894703255*alphaS*((0.16666666666666666 + 0.08333333333333333*Dminus4)*(Re(B0i(bb0,0.,0.,0.)) + Re(B0i(bb1,0.,0.,0.))) - 0.16666666666666666*Re(B0i(bb1,0.,MassGlu*MassGlu,MassSq*MassSq))) - 1.2732395447351628*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 1.909859317102744*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 2.5464790894703255*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.3183098861837907*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))) + 2.5464790894703255*alphaS*(-0.3333333333333333*Re(B0i(bb0,MassSq*MassSq,0.,MassGlu*MassGlu)) + 0.25*Re(B0i(bb0,MassSq*MassSq,0.,MassSq*MassSq)) + 0.16666666666666666*Re(B0i(bb1,MassSq*MassSq,0.,MassSq*MassSq)) - 0.3333333333333333*Re(B0i(bb1,MassSq*MassSq,MassGlu*MassGlu,0.)) - 0.3333333333333333*(MassSq*MassSq)*Re(B0i(dbb0,MassSq*MassSq,0.,MassGlu*MassGlu)) + 0.3333333333333333*(MassSq*MassSq)*Re(B0i(dbb0,MassSq*MassSq,0.,MassSq*MassSq)) + MassSq*MassSq*(0.16666666666666666*Re(B0i(dbb1,MassSq*MassSq,0.,MassSq*MassSq)) - 0.3333333333333333*Re(B0i(dbb1,MassSq*MassSq,MassGlu*MassGlu,0.)))))))/(-1.*(MassGlu*MassGlu) + T))))/(-1.*(MassGlu*MassGlu) + T) + (25.132741228718345*alphaS*(-1.*pow(MassSq,4) + T*U)*(16.*(alphaS*alphaS)*(((0.1875*Dminus4*(MassGlu*MassGlu - 1.*U) + 0.375*(3.*(MassGlu*MassGlu) - 1.*U))*B0i(bb0,U,0.,MassGlu*MassGlu) + (MassGlu*MassGlu + U)*(-1.*(0.375 + 0.1875*Dminus4)*B0i(bb1,U,0.,MassGlu*MassGlu) + 0.0625*(10.*B0i(bb1,U,0.,MassSq*MassSq) + 2.*B0i(bb1,U,MassTop*MassTop,MassSq*MassSq))))/pow(-1.*(MassGlu*MassGlu) + U,2) + ((0.375 + 0.09375*Dminus4)*B0i(bb0,U,0.,MassGlu*MassGlu) + 0.1875*(-1.*(MassSq*MassSq) + U)*C0i(cc0,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.020833333333333332*(MassGlu*MassGlu)*(C0i(cc0,MassSq*MassSq,U,0.,MassGlu*MassGlu,0.,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,U,0.,MassGlu*MassGlu,0.,MassSq*MassSq)) - 1.*(0.09375*Dminus4*(MassGlu*MassGlu - 1.*(MassSq*MassSq)) + 0.1875*(MassGlu*MassGlu - 1.*(MassSq*MassSq) - 1.*U))*C0i(cc1,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.) + 0.1875*(-1.*(MassSq*MassSq) + U)*C0i(cc2,U,MassSq*MassSq,0.,0.,MassGlu*MassGlu,0.))/(-1.*(MassGlu*MassGlu) + U) + 0.041666666666666664*D0i(dd00,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + 0.041666666666666664*(MassSq*MassSq)*D0i(dd11,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + 0.020833333333333332*(C0i(cc0,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + C0i(cc1,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + C0i(cc2,MassSq*MassSq,S,MassSq*MassSq,0.,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + 2.*T + U)*D0i(dd0,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + S + 3.*T + 2.*U)*D0i(dd1,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (3.*(MassSq*MassSq) + T)*D0i(dd12,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (3.*(MassSq*MassSq) + T)*D0i(dd13,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + MassSq*MassSq + 3.*T + U)*D0i(dd2,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassSq*MassSq + T)*D0i(dd22,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (2.*(MassSq*MassSq) + 2.*T)*D0i(dd23,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassGlu*MassGlu + MassSq*MassSq + 3.*T + U)*D0i(dd3,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq) + (MassSq*MassSq + T)*D0i(dd33,U,0.,S,MassSq*MassSq,MassSq*MassSq,0.,0.,MassGlu*MassGlu,MassSq*MassSq,MassSq*MassSq))) + 12.566370614359172*alphaS*((-1.*(-1.*MassGlu*(-0.954929658551372*alphaS*MassGlu*Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*MassGlu*Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + 0.15915494309189535*alphaS*MassGlu*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))) + MassGlu*(-0.6366197723675814*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 0.954929658551372*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.15915494309189535*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))))) + U*(-0.6366197723675814*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 0.954929658551372*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 1.2732395447351628*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.15915494309189535*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))))))/pow(-1.*(MassGlu*MassGlu) + U,2) + (0.01688686394038963*(12.566370614359172*alphaS*FiniteGs + 29.608813203268074*(4.*(-0.1193662073189215*alphaS*Divergence - 0.07957747154594767*alphaS*FiniteGs*(log((MassGlu*MassGlu)/(mu*mu)) + log((MassSq*MassSq)/(mu*mu)) + 0.3333333333333333*log((MassTop*MassTop)/(mu*mu)))) - 2.5464790894703255*alphaS*((0.16666666666666666 + 0.08333333333333333*Dminus4)*(Re(B0i(bb0,0.,0.,0.)) + Re(B0i(bb1,0.,0.,0.))) - 0.16666666666666666*Re(B0i(bb1,0.,MassGlu*MassGlu,MassSq*MassSq))) - 1.2732395447351628*alphaS*((0.75 + 0.375*Dminus4)*(Re(B0i(bb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) + Re(B0i(bb1,MassGlu*MassGlu,0.,MassGlu*MassGlu))) + 0.125*(-10.*Re(B0i(bb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(bb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq)))) + 1.909859317102744*alphaS*(MassGlu*MassGlu)*Re(B0i(dbb0,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 2.5464790894703255*alphaS*(0.75 + 0.375*Dminus4)*(MassGlu*MassGlu)*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassGlu*MassGlu)) - 0.3183098861837907*alphaS*(MassGlu*MassGlu)*(-10.*Re(B0i(dbb1,MassGlu*MassGlu,0.,MassSq*MassSq)) - 2.*Re(B0i(dbb1,MassGlu*MassGlu,MassTop*MassTop,MassSq*MassSq))) + 2.5464790894703255*alphaS*(-0.3333333333333333*Re(B0i(bb0,MassSq*MassSq,0.,MassGlu*MassGlu)) + 0.25*Re(B0i(bb0,MassSq*MassSq,0.,MassSq*MassSq)) + 0.16666666666666666*Re(B0i(bb1,MassSq*MassSq,0.,MassSq*MassSq)) - 0.3333333333333333*Re(B0i(bb1,MassSq*MassSq,MassGlu*MassGlu,0.)) - 0.3333333333333333*(MassSq*MassSq)*Re(B0i(dbb0,MassSq*MassSq,0.,MassGlu*MassGlu)) + 0.3333333333333333*(MassSq*MassSq)*Re(B0i(dbb0,MassSq*MassSq,0.,MassSq*MassSq)) + MassSq*MassSq*(0.16666666666666666*Re(B0i(dbb1,MassSq*MassSq,0.,MassSq*MassSq)) - 0.3333333333333333*Re(B0i(dbb1,MassSq*MassSq,MassGlu*MassGlu,0.)))))))/(-1.*(MassGlu*MassGlu) + U))))/(-1.*(MassGlu*MassGlu) + U);/*left - right*/
					// /* left-left */
					// /* right-right */
//    ltexi();    // error message of LoopTools 
    double matrixReal = matrix.real();
    return matrixReal;
}

// soft 

double Process::matrixMRSSMSoft_uu_suLsuRg(double s12, double th) {
   
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow(Alfas, 2);
   double beta = sqrt(1. - 4. * pow(m1, 2)/s12);
   double MGl2 = pow(MassGlu, 2);
   
   std::complex<double> temp = ((-42.666666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (42.666666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*MGl2*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
   (2.3703703703703704*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) - 
   (0.59259259259259259*Alfas*Alfas2*beta*(-1. + beta*beta)*(s12*s12)*(-1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
   (1.1851851851851852*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
   (0.59259259259259259*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (1.1851851851851852*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (42.666666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (42.666666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) - 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.14814814814814815*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (1.1851851851851852*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.59259259259259259*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.074074074074074074*Alfas*Alfas2*beta*Power(-1. + beta*beta,4)*Power(s12,3)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)
     - (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)
     - (9.4814814814814815*Alfas*Alfas2*(-1. + beta*beta)*(MGl2*MGl2)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (9.4814814814814815*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (4.7407407407407407*Alfas*Alfas2*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.59259259259259259*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (9.4814814814814815*Alfas*Alfas2*(MGl2*MGl2)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (9.4814814814814815*Alfas*Alfas2*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (4.7407407407407407*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.59259259259259259*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(s12,3)*Power(-1. + beta*Cos(th),3)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),3)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((mu_r*mu_r)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    + (0.018518518518518519*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    + (0.018518518518518519*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.018518518518518519*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.037037037037037037*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.018518518518518519*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    + (0.037037037037037037*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    + (0.037037037037037037*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    - (0.037037037037037037*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    - (1.1851851851851852*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.59259259259259259*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.074074074074074074*Alfas*Alfas2*beta*Power(-1. + beta*beta,4)*Power(s12,3)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*beta*(-1. + beta*beta)*(MGl2*MGl2)*s12*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.59259259259259259*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.074074074074074074*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.14814814814814815*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),3)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((mu_r*mu_r)/s12) + Power(Log((mu_r*mu_r)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.074074074074074074*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*(-1. + beta*beta)*s12*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.018518518518518519*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.037037037037037037*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.037037037037037037*Alfas*Alfas2*(-1. + beta*beta)*s12*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.018518518518518519*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.018518518518518519*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((mu_r*mu_r)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.0023148148148148148*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0046296296296296296*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0023148148148148148*Alfas*Alfas2*beta*(1. + beta*Cos(th))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.0040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) - 
   (0.032407407407407407*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.0040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.0034722222222222222*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.024305555555555556*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.016203703703703704*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.016203703703703704*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.032407407407407407*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),3)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0040509259259259259*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0081018518518518519*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.064814814814814815*Alfas*Alfas2*beta*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.016203703703703704*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.064814814814814815*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.016203703703703704*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.032407407407407407*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0081018518518518519*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0081018518518518519*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0011574074074074074*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.0046296296296296296*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((mu_r*mu_r)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((mu_r*mu_r)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((mu_r*mu_r)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2));
   
   return temp.real();
}
