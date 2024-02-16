double MRSSM::matrixTree_uubar_glglbar(const double alphas, const double S, const double T, const int Dminus4SeriesCoeff) const {
   const double MassGlu2 = pow<2>(MassGlu);
   const double MassSq2 = pow<2>(MassSq);
   const double U = 2.*MassGlu2 - T - S;
   double res = 0.;
   if (Dminus4SeriesCoeff == 0) {
      res = 4.*(8*(alphas*alphas)*(pi*pi)*((18*(3*(S*S) + pow<2>(T - U) + 2*S*(T + U)))/(S*S) + (2*pow<2>(S + T - U))/pow<2>(U - MassSq2) + (9*(3*(S*S) + 4*S*T + pow<2>(T - U)))/(2.*S*(U - MassSq2)) + (2*pow<2>(S - T + U))/pow<2>(T - MassSq2) + (9*(3*(S*S) + pow<2>(T - U) + 4*S*U))/(2.*S*(T - MassSq2))))/27.;
   }
   else if (Dminus4SeriesCoeff == 1) {
      res = 4.*(4*(alphas*alphas)*(pi*pi)*(4 + (S + T - U)/(U - MassSq2) + (S - T + U)/(T - MassSq2)))/3.;
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
