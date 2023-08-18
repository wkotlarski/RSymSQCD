double MRSSM::matrixTree_gg_suLsuLdagger(const double alphas, const double S, const double T, const int Dminus4SeriesCoeff) const {
   const double MassSq2 = pow<2>(MassSq);
   const double U = 2.*MassSq2 - T - S;
   double res = 0.;
   if (Dminus4SeriesCoeff == 0) {
      res = (alphas*alphas*(pi*pi)*(35*pow<6>(S) + 9*pow<6>(T - U) + 84*pow<5>(S)*(T + U) + 136*pow<3>(S)*pow<2>(T - U)*(T + U) + 36*S*pow<4>(T - U)*(T + U) + S*S*pow<2>(T - U)*(97*(T*T) + 94*T*U + 97*(U*U)) + pow<4>(S)*(115*(T*T) - 6*T*U + 115*(U*U))))/(6.*(S*S)*pow<2>(S*S - pow<2>(T - U)));
   }
   else if (Dminus4SeriesCoeff == 1) {
      res = 0;
   }
   else if (Dminus4SeriesCoeff == 2) {
      res = 0;
   }
   if(std::isinf(res)) {
      const double beta = std::abs(1- 4*MassSq2/S);
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