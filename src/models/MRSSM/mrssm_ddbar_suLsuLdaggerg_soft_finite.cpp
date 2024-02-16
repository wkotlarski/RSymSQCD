double MRSSM::matrixSoft_ddbar_suLsuLdaggerg_finite(const double alphas, const double S, const double th, const double dS, const double mu) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);
   const double lndS = std::log(dS);
   const double cth = std::cos(th);
   const double mu2 = pow<2>(mu);
   const double res = ((pow<2>(b) - pow<2>(b)*cos(2*th))*(4*(1 + pow<2>(b))*polylogarithm::Li2((2*b)/(1 + b)) + 8*b*(2*polylogarithm::Li2((b*(1 + cth))/(-1 + b)) + 7*polylogarithm::Li2((b - b*cth)/(-1 + b)) - 7*polylogarithm::Li2((b*(1 + cth))/(-1 + b*cth)) - 2*polylogarithm::Li2((b*(-1 + cth))/(1 + b*cth))) + 4*(1 + pow<2>(b))*lndS*log(-1 + 2/(1 + b)) + (1 - 18*b + pow<2>(b))*pow<2>(log(-1 + 2/(1 + b))) + 4*atanh(b)*(16 + (1 + pow<2>(b))*log(mu2/S)) + 4*b*(16*pow<2>(lndS) - 2*lndS*(8 + 9*log(1 - pow<2>(b))) + 7*pow<2>(log((-1 + b)/(-1 + b*cth))) + 2*pow<2>(log((1 - b)/(1 + b*cth))) + log(mu2/S)*(8 + 9*log(1 - pow<2>(b)) + 4*log(mu2/(pow<4>(dS)*S))) + 2*(7*log(1 - b*cth) + 2*log(1 + b*cth))*log((pow<2>(dS)*S)/mu2))))/(27.*b*pow<2>(pi));
   return pow<3>(alphas*pi)*res;
}