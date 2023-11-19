double MRSSM::matrixTree_ug_suLglbar(const double alphas, const double S, const double T, const int Dminus4SeriesCoeff) const {
   const double MassGlu2 = pow<2>(MassGlu);
   const double MassSq2 = pow<2>(MassSq);
   const double U = MassSq2 + MassGlu2 - T - S;
   double res = 0.;
   if (Dminus4SeriesCoeff == 0) {
      res = 2.*(-4*(alphas*alphas)*(pi*pi)*(4*(S*S) + S*(MassGlu2 - U) + 4*pow<2>(-MassGlu2 + U))*(5*pow<3>(MassGlu2) - pow<2>(MassGlu2)*(7*S + 6*T + 9*U) - U*(2*T*(S + T) + (S + 2*T)*U + U*U) + MassGlu2*(2*pow<2>(S + T) + 6*S*U + 8*T*U + 5*(U*U))))/(9.*S*pow<2>(MassGlu2 - T)*pow<2>(-MassGlu2 + U));
   }
   else if (Dminus4SeriesCoeff == 1) {
      res = 0;
   }
   else if (Dminus4SeriesCoeff == 2) {
      res = 0;
   }
   if(std::isinf(res)) {
      const double beta = std::abs(1- pow<2>(MassSq + MassGlu)/S);
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
