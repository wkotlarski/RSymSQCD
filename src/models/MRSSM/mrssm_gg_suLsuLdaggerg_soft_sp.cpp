double MRSSM::matrixSoft_gg_suLsuLdaggerg_sp(const double alphas, const double S, const double th, const double dS, const double mu) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);
   const double lndS = std::log(dS);
   const double cth = std::cos(th);
   const double mu2 = pow<2>(mu);
   const double res = -0.013888888888888888*((1 - 2*pow<2>(b) + 2*pow<4>(b) - 2*pow<4>(b)*pow<2>(cth) + pow<4>(b)*pow<4>(cth))*(22*(1 + pow<2>(b))*atanh(b) + b*(45*log(pow<2>(-1 + pow<2>(b)*pow<2>(cth))/pow<2>(-1 + pow<2>(b))) + 28*(-4 + 18*lndS - 9*log(mu2) + 9*log(S)) + 9*b*cth*(36*log(-1 + 2/(1 + b*cth)) + cth*((1 + pow<2>(b))*log(-1 + 2/(1 + b)) + b*(-16 + 72*lndS + 9*log(pow<2>(-1 + pow<2>(b)*pow<2>(cth))/pow<2>(-1 + pow<2>(b))) + 36*log(S/mu2)))))))/(b*pow<2>(-1 + pow<2>(b)*pow<2>(cth))*pow<2>(pi));
   return pow<3>(alphas*pi)*res;
}