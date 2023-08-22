double MRSSM::matrixSoft_uubar_suLsuLdaggerg_dp(const double alphas, const double S, const double th) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);
   const double cth = std::cos(th);
   const double MassGlu2 = pow<2>(MassGlu);
   const double res = (32*(pow<2>(b) - pow<2>(b)*cos(2*th))*(48*pow<2>(MassGlu2) + 8*(5 + 3*pow<2>(b))*MassGlu2*S + (31 + 16*pow<2>(b) + 3*pow<4>(b))*pow<2>(S) + 2*b*S*(-2*cth*(12*MassGlu2 + (5 + 3*pow<2>(b))*S) + 3*b*S*cos(2*th))))/(81.*pow<2>(pi)*pow<2>(4*MassGlu2 + S + pow<2>(b)*S - 2*b*cth*S));
   return pow<3>(alphas*pi)*res;
}