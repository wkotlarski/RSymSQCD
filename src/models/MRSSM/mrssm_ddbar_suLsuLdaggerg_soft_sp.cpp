double MRSSM::matrixSoft_ddbar_suLsuLdaggerg_sp(const double alphas, const double S, const double th, const double dS, const double mu) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);
   const double lndS = std::log(dS);
   const double cth = std::cos(th);
   const double mu2 = pow<2>(mu);
   const double res = (-2*(-pow<2>(b) + pow<2>(b)*cos(2*th))*((1 + pow<2>(b))*log((1 + b)/(1 - b)) + 18*b*log(1 - pow<2>(b)) - 4*b*(-4 + 8*lndS + 7*log(1 - b*cth) + 2*log(1 + b*cth) - 4*log(mu2) + 4*log(S))))/(27.*b*pow<2>(pi));
   return pow<3>(alphas*pi)*res;
}