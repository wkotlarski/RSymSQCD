double MRSSM::matrixSoft_uu_suLsuRg_dp(const double alphas, const double S, const double th) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);
   const double cth = std::cos(th);
   const double MassGlu2 = pow<2>(MassGlu);
   const double res = (512*pow<2>(S)*(pow<2>(b) - pow<2>(b)*cos(2*th))*(16*pow<2>(MassGlu2) + 8*(1 + pow<2>(b))*MassGlu2*S + (2*pow<2>(b) + pow<2>(1 + pow<2>(b)))*pow<2>(S) + 2*pow<2>(b)*pow<2>(S)*cos(2*th)))/(27.*pow<2>(pi)*pow<2>(-4*pow<2>(b)*pow<2>(cth)*pow<2>(S) + pow<2>(4*MassGlu2 + S + pow<2>(b)*S)));
   return pow<3>(alphas*pi)*res;
}