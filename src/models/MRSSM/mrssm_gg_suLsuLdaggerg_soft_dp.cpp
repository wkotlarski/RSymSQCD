double MRSSM::matrixSoft_gg_suLsuLdaggerg_dp(const double alphas, const double S, const double th) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);
   const double cth = std::cos(th);
   const double res = ((7 + 9*pow<2>(b)*pow<2>(cth))*(1 - 2*pow<2>(b) + 2*pow<4>(b) - 2*pow<4>(b)*pow<2>(cth) + pow<4>(b)*pow<4>(cth)))/(2.*pow<2>(-1 + pow<2>(b)*pow<2>(cth))*pow<2>(pi));
   return pow<3>(alphas*pi)*res;
}