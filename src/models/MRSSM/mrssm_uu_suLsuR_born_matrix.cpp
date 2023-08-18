double MRSSM::matrixTree_uu_suLsuR(const double alphas, const double S, const double T, const int Dminus4SeriesCoeff) const {
   const double MassGlu2 = pow<2>(MassGlu);
   const double MassSq2 = pow<2>(MassSq);
   const double U = 2.*MassSq2 - T - S;
   double res = 0.;
   if (Dminus4SeriesCoeff == 0) {
      res = (32*(alphas*alphas)*(pi*pi)*(-pow<2>(MassSq2) + T*U)*(2*pow<2>(MassGlu2) + T*T + U*U - 2*MassGlu2*(T + U)))/(9.*pow<2>(MassGlu2 - T)*pow<2>(MassGlu2 - U));
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
