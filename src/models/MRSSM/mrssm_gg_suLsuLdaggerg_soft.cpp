double MRSSM::matrixMRSSMSoft_gg_suLsuLdaggerg(double Alfas, double s12, double th, double dS, double muR) {
   const double Alfas2 = pow(Alfas, 2);
   const double b = sqrt(1. - 4.*Sqr(MassSq)/s12);
   double lndS = std::log(dS);
   double muR2 = Sqr(muR);
   const double m1 = MassSq;
   std::complex<double> temp = 0.000013563368055555555556*Alfas*Alfas2*pow(s12,-5.)*(2.*s12*PolyLog(2.,2.*b*pow(1. + b,-1.))*(-22. + 9.*pow(b,2.) + 9.*cos(2.*th)*pow(b,2.))*(s12 - 2.*pow(m1,2.))*(-64.*s12*pow(m1,2.) + 32.*s12*pow(b,2.)*pow(m1,2.) + 4.*s12*cos(2.*th)*pow(b,2.)*(s12*(-2. + pow(b,2.)) + 8.*pow(m1,2.)) + 256.*pow(m1,4.) + 8.*pow(s12,2.) - 8.*pow(b,2.)*pow(s12,2.) + 3.*pow(b,4.)*pow(s12,2.) + cos(4.*th)*pow(b,4.)*pow(s12,2.))*pow(-2. + pow(b,2.) + cos(2.*th)*pow(b,2.),2.) + 9.*b*PolyLog(2.,b*(1. + cos(th))*pow(-1. + b,-1.))*(10. - 36.*b*cos(th) + 9.*pow(b,2.) + 9.*cos(2.*th)*pow(b,2.))*pow(s12,2.)*(-64.*s12*pow(m1,2.) + 32.*s12*pow(b,2.)*pow(m1,2.) + 4.*s12*cos(2.*th)*pow(b,2.)*(s12*(-2. + pow(b,2.)) + 8.*pow(m1,2.)) + 256.*pow(m1,4.) + 8.*pow(s12,2.) - 8.*pow(b,2.)*pow(s12,2.) + 3.*pow(b,4.)*pow(s12,2.) + cos(4.*th)*pow(b,4.)*pow(s12,2.))*pow(-2. + pow(b,2.) + cos(2.*th)*pow(b,2.),2.) + 32.*(4352.*b*s12*pow(m1,4.) - 37888.*b*pow(m1,6.) - 2816.*lndS*s12*log((1. + b)*pow(1. - 1.*b,-1.))*pow(m1,6.) + 1408.*s12*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(m1,6.) + 7364.*b*pow(m1,2.)*pow(s12,2.) - 7168.*b*lndS*pow(m1,4.)*pow(s12,2.) + 3584.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(m1,4.)*pow(s12,2.) + 2112.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(m1,4.)*pow(s12,2.) + 3584.*b*log(muR2*pow(s12,-1.))*pow(m1,4.)*pow(s12,2.) - 16128.*b*lndS*log(muR2*pow(s12,-1.))*pow(m1,4.)*pow(s12,2.) - 1056.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(m1,4.)*pow(s12,2.) + 2880.*b*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(m1,4.)*pow(s12,2.) - 1440.*b*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(m1,4.)*pow(s12,2.) + 2880.*b*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(m1,4.)*pow(s12,2.) - 1440.*b*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(m1,4.)*pow(s12,2.) - 2880.*b*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(m1,4.)*pow(s12,2.) + 10368.*lndS*cos(th)*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(m1,4.)*pow(s12,2.) - 5184.*cos(th)*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(m1,4.)*pow(s12,2.) - 10368.*lndS*cos(th)*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(m1,4.)*pow(s12,2.) + 5184.*cos(th)*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(m1,4.)*pow(s12,2.) + 10368.*cos(th)*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,2.)*pow(m1,4.)*pow(s12,2.) + 16128.*b*pow(lndS,2.)*pow(m1,4.)*pow(s12,2.) - 1269.*b*pow(s12,3.) + 1792.*b*lndS*pow(m1,2.)*pow(s12,3.) - 896.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(m1,2.)*pow(s12,3.) - 440.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(m1,2.)*pow(s12,3.) - 896.*b*log(muR2*pow(s12,-1.))*pow(m1,2.)*pow(s12,3.) + 4032.*b*lndS*log(muR2*pow(s12,-1.))*pow(m1,2.)*pow(s12,3.) + 220.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(m1,2.)*pow(s12,3.) - 720.*b*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(m1,2.)*pow(s12,3.) + 360.*b*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(m1,2.)*pow(s12,3.) - 720.*b*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(m1,2.)*pow(s12,3.) + 360.*b*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(m1,2.)*pow(s12,3.) + 720.*b*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(m1,2.)*pow(s12,3.) - 2592.*lndS*cos(th)*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(m1,2.)*pow(s12,3.) + 1296.*cos(th)*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(m1,2.)*pow(s12,3.) + 2592.*lndS*cos(th)*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(m1,2.)*pow(s12,3.) - 1296.*cos(th)*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(m1,2.)*pow(s12,3.) - 2592.*cos(th)*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,2.)*pow(m1,2.)*pow(s12,3.) - 4032.*b*pow(lndS,2.)*pow(m1,2.)*pow(s12,3.) - 224.*b*lndS*pow(s12,4.) + 112.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(s12,4.) + 44.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(s12,4.) + 112.*b*log(muR2*pow(s12,-1.))*pow(s12,4.) - 504.*b*lndS*log(muR2*pow(s12,-1.))*pow(s12,4.) - 22.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(s12,4.) + 90.*b*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(s12,4.) - 45.*b*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(s12,4.) + 90.*b*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(s12,4.) - 45.*b*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(s12,4.) - 90.*b*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(s12,4.) + 324.*lndS*cos(th)*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(s12,4.) - 162.*cos(th)*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(s12,4.) - 324.*lndS*cos(th)*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(s12,4.) + 162.*cos(th)*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,2.)*pow(s12,4.) + 324.*cos(th)*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,2.)*pow(s12,4.) + 504.*b*pow(lndS,2.)*pow(s12,4.) - 512.*s12*pow(b,3.)*pow(m1,4.)*pow(cos(th),2.) + 7936.*lndS*s12*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,2.)*pow(m1,6.)*pow(cos(th),2.) - 3968.*s12*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,2.)*pow(m1,6.)*pow(cos(th),2.) - 186368.*pow(b,3.)*pow(m1,6.)*pow(cos(th),2.) - 22576.*pow(b,3.)*pow(m1,2.)*pow(s12,2.)*pow(cos(th),2.) - 2560.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,2.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) - 6656.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,2.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) + 3328.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,2.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) + 5120.*lndS*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) - 2560.*log(muR2*pow(s12,-1.))*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) + 11520.*lndS*log(muR2*pow(s12,-1.))*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) - 576.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) + 288.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) - 576.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) + 288.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) + 576.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) - 11520.*pow(b,3.)*pow(lndS,2.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.) + 4744.*pow(b,3.)*pow(s12,3.)*pow(cos(th),2.) + 1536.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,2.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) + 1768.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,2.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) - 884.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,2.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) - 3072.*lndS*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) + 1536.*log(muR2*pow(s12,-1.))*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) - 6912.*lndS*log(muR2*pow(s12,-1.))*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) + 864.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) - 432.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) + 864.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) - 432.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) - 864.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) + 6912.*pow(b,3.)*pow(lndS,2.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.) - 304.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,2.)*pow(s12,4.)*pow(cos(th),2.) - 212.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,2.)*pow(s12,4.)*pow(cos(th),2.) + 106.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,2.)*pow(s12,4.)*pow(cos(th),2.) + 608.*lndS*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.) - 304.*log(muR2*pow(s12,-1.))*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.) + 1368.*lndS*log(muR2*pow(s12,-1.))*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.) - 198.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.) + 99.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.) - 198.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.) + 99.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.) + 198.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.) - 1368.*pow(b,3.)*pow(lndS,2.)*pow(s12,4.)*pow(cos(th),2.) - 20736.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),3.) + 10368.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),3.) + 20736.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),3.) - 10368.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),3.) - 20736.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),3.) + 7776.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),3.) - 3888.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),3.) - 7776.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),3.) + 3888.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),3.) + 7776.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),3.) - 1296.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(s12,4.)*pow(cos(th),3.) + 648.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(s12,4.)*pow(cos(th),3.) + 1296.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(s12,4.)*pow(cos(th),3.) - 648.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,4.)*pow(s12,4.)*pow(cos(th),3.) - 1296.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,4.)*pow(s12,4.)*pow(cos(th),3.) - 3840.*s12*pow(b,5.)*pow(m1,4.)*pow(cos(th),4.) - 7424.*lndS*s12*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,4.)*pow(m1,6.)*pow(cos(th),4.) + 3712.*s12*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,4.)*pow(m1,6.)*pow(cos(th),4.) - 37888.*pow(b,5.)*pow(m1,6.)*pow(cos(th),4.) + 28568.*pow(b,5.)*pow(m1,2.)*pow(s12,2.)*pow(cos(th),4.) - 5632.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) + 7552.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) - 3776.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) + 11264.*lndS*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) - 5632.*log(muR2*pow(s12,-1.))*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) + 25344.*lndS*log(muR2*pow(s12,-1.))*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) - 7488.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) + 3744.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) - 7488.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) + 3744.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) + 7488.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) - 25344.*pow(b,5.)*pow(lndS,2.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.) - 5672.*pow(b,5.)*pow(s12,3.)*pow(cos(th),4.) + 768.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) - 2736.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) + 1368.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) - 1536.*lndS*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) + 768.*log(muR2*pow(s12,-1.))*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) - 3456.*lndS*log(muR2*pow(s12,-1.))*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) + 1728.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) - 864.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) + 1728.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) - 864.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) - 1728.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) + 3456.*pow(b,5.)*pow(lndS,2.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.) + 96.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,4.)*pow(s12,4.)*pow(cos(th),4.) + 408.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,4.)*pow(s12,4.)*pow(cos(th),4.) - 204.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,4.)*pow(s12,4.)*pow(cos(th),4.) - 192.*lndS*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.) + 96.*log(muR2*pow(s12,-1.))*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.) - 432.*lndS*log(muR2*pow(s12,-1.))*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.) - 108.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.) + 54.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.) - 108.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.) + 54.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.) + 108.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.) + 432.*pow(b,5.)*pow(lndS,2.)*pow(s12,4.)*pow(cos(th),4.) + 10368.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),5.) - 5184.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),5.) - 10368.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),5.) + 5184.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),5.) + 10368.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),5.) - 7776.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),5.) + 3888.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),5.) + 7776.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),5.) - 3888.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),5.) - 7776.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),5.) + 1944.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(s12,4.)*pow(cos(th),5.) - 972.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(s12,4.)*pow(cos(th),5.) - 1944.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(s12,4.)*pow(cos(th),5.) + 972.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,6.)*pow(s12,4.)*pow(cos(th),5.) + 1944.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,6.)*pow(s12,4.)*pow(cos(th),5.) + 2304.*lndS*s12*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,6.)*pow(m1,6.)*pow(cos(th),6.) - 1152.*s12*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,6.)*pow(m1,6.)*pow(cos(th),6.) - 18864.*pow(b,7.)*pow(m1,2.)*pow(s12,2.)*pow(cos(th),6.) + 4608.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) - 3584.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) + 1792.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) - 9216.*lndS*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) + 4608.*log(muR2*pow(s12,-1.))*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) - 20736.*lndS*log(muR2*pow(s12,-1.))*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) + 5184.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) - 2592.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) + 5184.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) - 2592.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) - 5184.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) + 20736.*pow(b,7.)*pow(lndS,2.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.) + 2026.*pow(b,7.)*pow(s12,3.)*pow(cos(th),6.) - 2560.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) + 2000.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) - 1000.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) + 5120.*lndS*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) - 2560.*log(muR2*pow(s12,-1.))*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) + 11520.*lndS*log(muR2*pow(s12,-1.))*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) - 3168.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) + 1584.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) - 3168.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) + 1584.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) + 3168.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) - 11520.*pow(b,7.)*pow(lndS,2.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.) + 416.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,6.)*pow(s12,4.)*pow(cos(th),6.) - 392.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,6.)*pow(s12,4.)*pow(cos(th),6.) + 196.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,6.)*pow(s12,4.)*pow(cos(th),6.) - 832.*lndS*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.) + 416.*log(muR2*pow(s12,-1.))*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.) - 1872.*lndS*log(muR2*pow(s12,-1.))*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.) + 612.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.) - 306.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.) + 612.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.) - 306.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.) - 612.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.) + 1872.*pow(b,7.)*pow(lndS,2.)*pow(s12,4.)*pow(cos(th),6.) + 2592.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),7.) - 1296.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),7.) - 2592.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),7.) + 1296.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),7.) + 2592.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),7.) - 1296.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,8.)*pow(s12,4.)*pow(cos(th),7.) + 648.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,8.)*pow(s12,4.)*pow(cos(th),7.) + 1296.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,8.)*pow(s12,4.)*pow(cos(th),7.) - 648.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,8.)*pow(s12,4.)*pow(cos(th),7.) - 1296.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,8.)*pow(s12,4.)*pow(cos(th),7.) + 5508.*pow(b,9.)*pow(m1,2.)*pow(s12,2.)*pow(cos(th),8.) + 576.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,8.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),8.) - 288.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,8.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),8.) + 333.*pow(b,9.)*pow(s12,3.)*pow(cos(th),8.) + 1152.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) - 664.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) + 332.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) - 2304.*lndS*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) + 1152.*log(muR2*pow(s12,-1.))*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) - 5184.*lndS*log(muR2*pow(s12,-1.))*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) + 1296.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) - 648.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) + 1296.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) - 648.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) - 1296.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) + 5184.*pow(b,9.)*pow(lndS,2.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.) - 464.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,8.)*pow(s12,4.)*pow(cos(th),8.) + 188.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,8.)*pow(s12,4.)*pow(cos(th),8.) - 94.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,8.)*pow(s12,4.)*pow(cos(th),8.) + 928.*lndS*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.) - 464.*log(muR2*pow(s12,-1.))*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.) + 2088.*lndS*log(muR2*pow(s12,-1.))*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.) - 558.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.) + 279.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.) - 558.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.) + 279.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.) + 558.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.) - 2088.*pow(b,9.)*pow(lndS,2.)*pow(s12,4.)*pow(cos(th),8.) + 324.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,10.)*pow(s12,4.)*pow(cos(th),9.) - 162.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,10.)*pow(s12,4.)*pow(cos(th),9.) - 324.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,10.)*pow(s12,4.)*pow(cos(th),9.) + 162.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,10.)*pow(s12,4.)*pow(cos(th),9.) + 324.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,10.)*pow(s12,4.)*pow(cos(th),9.) - 162.*pow(b,11.)*pow(s12,3.)*pow(cos(th),10.) + 72.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,10.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),10.) - 36.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,10.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),10.) + 144.*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,10.)*pow(s12,4.)*pow(cos(th),10.) - 36.*lndS*log((1. + b)*pow(1. - 1.*b,-1.))*pow(b,10.)*pow(s12,4.)*pow(cos(th),10.) + 18.*log((1. + b)*pow(1. - 1.*b,-1.))*log(muR2*pow(s12,-1.))*pow(b,10.)*pow(s12,4.)*pow(cos(th),10.) - 288.*lndS*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.) + 144.*log(muR2*pow(s12,-1.))*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.) - 648.*lndS*log(muR2*pow(s12,-1.))*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.) + 162.*lndS*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.) - 81.*log(muR2*pow(s12,-1.))*log(-1.*pow(-1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.) + 162.*lndS*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.) - 81.*log(muR2*pow(s12,-1.))*log(-1.*pow(1. + b*cos(th),2.)*pow(-1. + pow(b,2.),-1.))*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.) - 162.*PolyLog(2.,b*(-1. + cos(th))*pow(1. + b*cos(th),-1.))*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.) + 648.*pow(b,11.)*pow(lndS,2.)*pow(s12,4.)*pow(cos(th),10.) + 704.*s12*pow(m1,6.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 528.*pow(m1,4.)*pow(s12,2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 1440.*b*pow(m1,4.)*pow(s12,2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 110.*pow(m1,2.)*pow(s12,3.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 360.*b*pow(m1,2.)*pow(s12,3.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 11.*pow(s12,4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 45.*b*pow(s12,4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 1984.*s12*pow(b,2.)*pow(m1,6.)*pow(cos(th),2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 1664.*pow(b,2.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 288.*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 442.*pow(b,2.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 432.*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 53.*pow(b,2.)*pow(s12,4.)*pow(cos(th),2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 99.*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 1856.*s12*pow(b,4.)*pow(m1,6.)*pow(cos(th),4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 1888.*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 3744.*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 684.*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 864.*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 102.*pow(b,4.)*pow(s12,4.)*pow(cos(th),4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 54.*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 576.*s12*pow(b,6.)*pow(m1,6.)*pow(cos(th),6.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 896.*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 2592.*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 500.*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 1584.*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 98.*pow(b,6.)*pow(s12,4.)*pow(cos(th),6.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 306.*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 144.*pow(b,8.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),8.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 166.*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 648.*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 47.*pow(b,8.)*pow(s12,4.)*pow(cos(th),8.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 279.*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 18.*pow(b,10.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),10.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 9.*pow(b,10.)*pow(s12,4.)*pow(cos(th),10.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) - 81.*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.)*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.) + 4032.*b*pow(m1,4.)*pow(s12,2.)*pow(log(muR2*pow(s12,-1.)),2.) - 1008.*b*pow(m1,2.)*pow(s12,3.)*pow(log(muR2*pow(s12,-1.)),2.) + 126.*b*pow(s12,4.)*pow(log(muR2*pow(s12,-1.)),2.) - 2880.*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.)*pow(log(muR2*pow(s12,-1.)),2.) + 1728.*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.)*pow(log(muR2*pow(s12,-1.)),2.) - 342.*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.)*pow(log(muR2*pow(s12,-1.)),2.) - 6336.*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.)*pow(log(muR2*pow(s12,-1.)),2.) + 864.*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.)*pow(log(muR2*pow(s12,-1.)),2.) + 108.*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.)*pow(log(muR2*pow(s12,-1.)),2.) + 5184.*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.)*pow(log(muR2*pow(s12,-1.)),2.) - 2880.*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.)*pow(log(muR2*pow(s12,-1.)),2.) + 468.*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.)*pow(log(muR2*pow(s12,-1.)),2.) + 1296.*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.)*pow(log(muR2*pow(s12,-1.)),2.) - 522.*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.)*pow(log(muR2*pow(s12,-1.)),2.) + 162.*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.)*pow(log(muR2*pow(s12,-1.)),2.) + 1440.*b*pow(m1,4.)*pow(s12,2.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 5184.*cos(th)*pow(b,2.)*pow(m1,4.)*pow(s12,2.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 360.*b*pow(m1,2.)*pow(s12,3.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 1296.*cos(th)*pow(b,2.)*pow(m1,2.)*pow(s12,3.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 45.*b*pow(s12,4.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 162.*cos(th)*pow(b,2.)*pow(s12,4.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 288.*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 432.*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 99.*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 10368.*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),3.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 3888.*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),3.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 648.*pow(b,4.)*pow(s12,4.)*pow(cos(th),3.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 3744.*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 864.*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 54.*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 5184.*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),5.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 3888.*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),5.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 972.*pow(b,6.)*pow(s12,4.)*pow(cos(th),5.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 2592.*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 1584.*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 306.*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 1296.*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),7.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 648.*pow(b,8.)*pow(s12,4.)*pow(cos(th),7.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 648.*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) - 279.*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 162.*pow(b,10.)*pow(s12,4.)*pow(cos(th),9.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 81.*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.)*pow(log((-1. + b)*pow(-1. + b*cos(th),-1.)),2.) + 1440.*b*pow(m1,4.)*pow(s12,2.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 5184.*cos(th)*pow(b,2.)*pow(m1,4.)*pow(s12,2.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 360.*b*pow(m1,2.)*pow(s12,3.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 1296.*cos(th)*pow(b,2.)*pow(m1,2.)*pow(s12,3.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 45.*b*pow(s12,4.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 162.*cos(th)*pow(b,2.)*pow(s12,4.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 288.*pow(b,3.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),2.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 432.*pow(b,3.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),2.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 99.*pow(b,3.)*pow(s12,4.)*pow(cos(th),2.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 10368.*pow(b,4.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),3.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 3888.*pow(b,4.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),3.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 648.*pow(b,4.)*pow(s12,4.)*pow(cos(th),3.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 3744.*pow(b,5.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),4.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 864.*pow(b,5.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),4.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 54.*pow(b,5.)*pow(s12,4.)*pow(cos(th),4.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 5184.*pow(b,6.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),5.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 3888.*pow(b,6.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),5.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 972.*pow(b,6.)*pow(s12,4.)*pow(cos(th),5.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 2592.*pow(b,7.)*pow(m1,4.)*pow(s12,2.)*pow(cos(th),6.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 1584.*pow(b,7.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),6.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 306.*pow(b,7.)*pow(s12,4.)*pow(cos(th),6.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 1296.*pow(b,8.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),7.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 648.*pow(b,8.)*pow(s12,4.)*pow(cos(th),7.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 648.*pow(b,9.)*pow(m1,2.)*pow(s12,3.)*pow(cos(th),8.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 279.*pow(b,9.)*pow(s12,4.)*pow(cos(th),8.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) - 162.*pow(b,10.)*pow(s12,4.)*pow(cos(th),9.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 81.*pow(b,11.)*pow(s12,4.)*pow(cos(th),10.)*pow(log((1. - 1.*b)*pow(1. + b*cos(th),-1.)),2.) + 0.28125*b*PolyLog(2.,(b - 1.*b*cos(th))*pow(-1. + b,-1.))*(10. + 36.*b*cos(th) + 9.*pow(b,2.) + 9.*cos(2.*th)*pow(b,2.))*pow(s12,2.)*(-64.*s12*pow(m1,2.) + 32.*s12*pow(b,2.)*pow(m1,2.) + 4.*s12*cos(2.*th)*pow(b,2.)*(s12*(-2. + pow(b,2.)) + 8.*pow(m1,2.)) + 256.*pow(m1,4.) + 8.*pow(s12,2.) - 8.*pow(b,2.)*pow(s12,2.) + 3.*pow(b,4.)*pow(s12,2.) + cos(4.*th)*pow(b,4.)*pow(s12,2.))*pow(-2. + pow(b,2.) + cos(2.*th)*pow(b,2.),2.) - 0.28125*b*PolyLog(2.,b*(1. + cos(th))*pow(-1. + b*cos(th),-1.))*(10. + 36.*b*cos(th) + 9.*pow(b,2.) + 9.*cos(2.*th)*pow(b,2.))*pow(s12,2.)*(-64.*s12*pow(m1,2.) + 32.*s12*pow(b,2.)*pow(m1,2.) + 4.*s12*cos(2.*th)*pow(b,2.)*(s12*(-2. + pow(b,2.)) + 8.*pow(m1,2.)) + 256.*pow(m1,4.) + 8.*pow(s12,2.) - 8.*pow(b,2.)*pow(s12,2.) + 3.*pow(b,4.)*pow(s12,2.) + cos(4.*th)*pow(b,4.)*pow(s12,2.))*pow(-2. + pow(b,2.) + cos(2.*th)*pow(b,2.),2.)))*pow(-1. + pow(b,2.)*pow(cos(th),2.),-4.)*Sin(th);
   return temp.real();
}