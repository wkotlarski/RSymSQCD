double Process::matrixMRSSMTree_GG_suLsuLdagger(double alphas, double S, double T) const { // agrees with Philip
   const double U = 2*MassSq*MassSq - S - T;
   static constexpr double k = 2.*2*8*8;
   static constexpr double h = 1.;
   return  h/k*(-157.91367041742973*Sqr(alphas)*((96.*(pow(MassSq,4) + T*U - 1.*(MassSq*MassSq)*(T + U)))/(S*S) + MassSq*MassSq*(37.333333333333336/(-1.*(MassSq*MassSq) + U) - 85.33333333333333*(T/pow(-1.*(MassSq*MassSq) + T,2) + U/pow(-1.*(MassSq*MassSq) + U,2))) + (37.333333333333336*(MassSq*MassSq) + (5.333333333333333*(9.*pow(MassSq,4) - 3.*(MassSq*MassSq)*(2.*(MassSq*MassSq) + S) + (S + T)*(S + U)))/(-1.*(MassSq*MassSq) + U))/(-1.*(MassSq*MassSq) + T) - (48.*((5.*pow(MassSq,4) + MassSq*MassSq*U - 1.*T*(-1.*(MassSq*MassSq) + U) - 1.*(MassSq*MassSq)*(3.*S + 4.*T + 2.*U))/(-1.*(MassSq*MassSq) + U) + (5.*pow(MassSq,4) + MassSq*MassSq*U - 1.*T*(-1.*(MassSq*MassSq) + U) - 1.*(MassSq*MassSq)*(3.*S + 2.*T + 4.*U))/(-1.*(MassSq*MassSq) + T)))/S));
}

double Process::matrixMRSSMTree_ddbar_suLsuLdagger(double alphas, double S, double T ) const { // agrees with Philip
   const double U = 2*MassSq*MassSq - S - T;
   static constexpr double k = 2.*2*3*3;
   static constexpr double h = 2.*2;
   return  h/k*(-631.6546816697189*Sqr(alphas)*(pow(MassSq,4) - T*U))/Sqr(S);
}

double Process::matrixMRSSMTree_uubar_suLsuLdagger(double alphas, double S, double T) const {
   const double U = 2*MassSq*MassSq - S - T;
   static constexpr double k = 2.*2*3*3;
   static constexpr double h = 2.*2;
   return  h/k*((315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/(S*S) + 355.3057584392169*Sqr(alphas)*pow(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U) - 236.8705056261446*Sqr(alphas)*(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T))*(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T))*(-1.*pow(MassSq,4) + T*U) + 355.3057584392169*Sqr(alphas)*pow(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U));
}

double Process::matrixMRSSMTree_udbar_suLsdLdagger(double alphas, double S, double T) const {
   double U = 2*MassSq*MassSq - S - T;
   return  (315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2);
}

double Process::matrixMRSSMTree_uu_suLsuR(double alphas, double S, double T) const {
   static constexpr double h = 2.*2;
   static constexpr double k = 2.*2*3*3;
   const double U = 2*MassSq*MassSq - S - T;
   return h/k*((315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2) +
           (315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + U,2));
}
