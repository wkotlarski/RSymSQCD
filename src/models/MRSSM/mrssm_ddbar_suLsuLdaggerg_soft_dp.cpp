double MRSSM::matrixSoft_ddbar_suLsuLdaggerg_dp(const double alphas, const double S, const double th) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);
   const double cth = std::cos(th);
   const double res = (64*(pow<2>(b) - pow<2>(b)*pow<2>(cth)))/(27.*pow<2>(pi));
   return pow<3>(alphas*pi)*res;
}