double MRSSM::matrixSoft_uu_suLsuRg_sp(const double alphas, const double S, const double th, const double dS, const double mu) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);
   const double lndS = std::log(dS);
   const double cth = std::cos(th);
   const double MassGlu2 = pow<2>(MassGlu);
   const double mu2 = pow<2>(mu);
   const double res = (-64*pow<2>(S)*(pow<2>(b) - pow<2>(b)*cos(2*th))*(-32*pow<2>(b)*cth*S*(4*MassGlu2 + S + pow<2>(b)*S)*log(-1 + 2/(1 + b*cth)) + (16*pow<2>(MassGlu2) + 8*(1 + pow<2>(b))*MassGlu2*S + (2*pow<2>(b) + pow<2>(1 + pow<2>(b)))*pow<2>(S))*(2*(1 + pow<2>(b))*atanh(b) + 2*b*(-4 + 8*lndS - 3*log(1 - pow<2>(b)) + 3*log(1 - pow<2>(b)*pow<2>(cth)) - 4*log(mu2) + 4*log(S))) + 4*pow<2>(b)*pow<2>(S)*cos(2*th)*((1 + pow<2>(b))*atanh(b) + b*(-4 + 8*lndS - 3*log(1 - pow<2>(b)) + 3*log(1 - pow<2>(b)*pow<2>(cth)) + 4*log(S/mu2)))))/(27.*b*pow<2>(pi)*pow<2>(-4*pow<2>(b)*pow<2>(cth)*pow<2>(S) + pow<2>(4*MassGlu2 + S + pow<2>(b)*S)));
   return pow<3>(alphas*pi)*res;
}