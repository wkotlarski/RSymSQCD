double Process::matrixSimplifiedSoft_uubar_OOg( double s12, double th) {
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow(Alfas, 2);
   double b = sqrt(1. - 4. * pow(m1, 2)/s12);
   double lndS = log(dS);
   double muR = mu_r;
   std::complex<double> temp = (0.0069444444444444444444*Alfas*Alfas2*(b*b)*(-144.*b - 144.*b*lndS*s12 + 64.*b*(lndS*lndS)*s12 + 72.*s12*Log((1. + b)/(1. - 1.*b)) + 36.*lndS*s12*Log((1. + b)/(1. - 1.*b)) + 36.*(b*b)*lndS*s12*Log((1. + b)/(1. - 1.*b)) - 9.*s12*Power(Log((1. + b)/(1. - 1.*b)),2) - 18.*b*s12*Power(Log((1. + b)/(1. - 1.*b)),2) - 9.*(b*b)*s12*Power(Log((1. + b)/(1. - 1.*b)),2) + 72.*b*s12*Log((muR*muR)/s12) - 64.*b*lndS*s12*Log((muR*muR)/s12) - 18.*s12*Log((1. + b)/(1. - 1.*b))*Log((muR*muR)/s12) - 18.*(b*b)*s12*Log((1. + b)/(1. - 1.*b))*Log((muR*muR)/s12) + 16.*b*s12*Power(Log((muR*muR)/s12),2) + 18.*b*s12*Power(Log((-1. + b)/(-1. + b*Cos(th))),2) + 36.*b*lndS*s12*Log((-1.*Power(-1. + b*Cos(th),2))/(-1. + b*b)) - 18.*b*s12*Log((muR*muR)/s12)*Log((-1.*Power(-1. + b*Cos(th),2))/(-1. + b*b)) + 18.*b*s12*Power(Log((1. - 1.*b)/(1. + b*Cos(th))),2) + 36.*b*lndS*s12*Log((-1.*Power(1. + b*Cos(th),2))/(-1. + b*b)) - 18.*b*s12*Log((muR*muR)/s12)*Log((-1.*Power(1. + b*Cos(th),2))/(-1. + b*b)) - 36.*(1. + b*b)*s12*PolyLog(2.,(2.*b)/(1. + b)) + 36.*b*s12*PolyLog(2.,(b*(1. + Cos(th)))/(-1. + b)) + 36.*b*s12*PolyLog(2.,(b - 1.*b*Cos(th))/(-1. + b)) - 36.*b*s12*PolyLog(2.,(b*(1. + Cos(th)))/(-1. + b*Cos(th))) - 36.*b*s12*PolyLog(2.,(b*(-1. + Cos(th)))/(1. + b*Cos(th))))*Power(Sin(th),3))/(s12*s12);
   return temp.real();
}