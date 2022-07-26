inline double Process::matrixMRSSMTree_uubar_suLsuLdagger(double alphas, double S, double T) const {
   double U = 2*MassSq*MassSq - S - T;
   return  (315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/(S*S) + 355.3057584392169*Sqr(alphas)*pow(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U) - 236.8705056261446*Sqr(alphas)*(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T))*(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T))*(-1.*pow(MassSq,4) + T*U) + 355.3057584392169*Sqr(alphas)*pow(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U);
}
