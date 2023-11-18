double MRSSM::matrixTree_gg_glglbar(const double alphas, const double S, const double T, const int Dminus4SeriesCoeff) const {
   const double MassGlu2 = pow<2>(MassGlu);
   const double U = 2.*MassGlu2 - T - S;
   double res = 0.;
   if (Dminus4SeriesCoeff == 0) {
      res = 4.*(-9*(alphas*alphas)*(pi*pi)*(3*(S*S) + pow<2>(T - U))*(3*pow<4>(S) + pow<4>(T - U) + 12*pow<3>(S)*(T + U) + 4*S*pow<2>(T - U)*(T + U) + 4*(S*S)*(3*(T*T) + 2*T*U + 3*(U*U))))/(32.*(S*S)*pow<2>(MassGlu2 - T)*pow<2>(MassGlu2 - U));
   }
   else if (Dminus4SeriesCoeff == 1) {
      res = 0;
   }
   else if (Dminus4SeriesCoeff == 2) {
      res = 0;
   }
   if(std::isinf(res)) {
      const double beta = std::abs(1- 4*MassGlu2/S);
      if (beta < 1e-10) {
         spdlog::get("console")->warn("Infinity in the 2->2 |ME|^2 at threshold (β = {}). Returning 0.", beta);
         return 0.;
      }
      else {
         spdlog::get("console")->error("Infinity in the 2->2 |ME|^2 for β = {}", beta);
         return res;
      }
   }
   else {
      return res;
   }
}
