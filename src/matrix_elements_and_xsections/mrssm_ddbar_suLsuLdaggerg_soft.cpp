double Process::matrixMRSSMSoft_ddbar_suLsuLdaggerg( double s12, double th) {
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow(Alfas, 2);
   double b = sqrt(1. - 4. * pow(m1, 2)/s12);
   double muR = mu_r;
   std::complex<double> temp = (-0.0023148148148148147*Alfas*Alfas2*(8.*(m1*m1) + (-2. + b*b)*s12 + b*b*s12*cos(2.*th))*(16.*s12*log((1. + b)/(1. - 1.*b)) + (-2.*(m1*m1) + s12 - 9.*b*s12)*Power(log((1. + b)/(1. - 1.*b)),2) + 32.*b*s12*Power(log(dS),2) + 4.*log(dS)*(2.*(m1*m1)*log((1. + b)/(1. - 1.*b)) + s12*(log(-1. + 2./(1. + b)) + b*(-8. + 8.*log(s12/(muR*muR)) + 7.*log((-1.*Power(-1. + b*cos(th),2))/(-1. + b*b))))) + 2.*((-2.*(m1*m1) + s12)*log(-1. + 2./(1. + b))*log(s12/(muR*muR)) + b*(-36. + s12*(log(Power(muR,16)/Power(s12,8)) + 4.*Power(log(s12/(muR*muR)),2) + 7.*Power(log((-1. + b)/(-1. + b*cos(th))),2) - 7.*log((muR*muR)/s12)*log((-1.*Power(-1. + b*cos(th),2))/(-1. + b*b)) + 2.*Power(log((1. - 1.*b)/(1. + b*cos(th))),2) + 2.*log((dS*dS*s12)/(muR*muR))*log((-1.*Power(1. + b*cos(th),2))/(-1. + b*b))))) + 4.*((-2.*(m1*m1) + s12)*PolyLog(2.,(2.*b)/(1. + b)) + b*s12*(2.*PolyLog(2.,(b*(1. + cos(th)))/(-1. + b)) + 7.*PolyLog(2.,(b - 1.*b*cos(th))/(-1. + b)) - 7.*PolyLog(2.,(b*(1. + cos(th)))/(-1. + b*cos(th))) - 2.*PolyLog(2.,(b*(-1. + cos(th)))/(1. + b*cos(th))))))*Sin(th))/Power(s12,3);
   //std::cout << std::scientific;
   //std::cout << std::setprecision(18);
   //std::cout << Alfas << ' ' << s12 << ' '<< th << ' ' << temp.real() << ' ' << temp.imag() << ' ' << temp.imag()/temp.real() << std::endl;
   return temp.real();
}