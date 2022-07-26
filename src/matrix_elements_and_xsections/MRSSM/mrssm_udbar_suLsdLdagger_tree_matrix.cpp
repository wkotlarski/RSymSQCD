inline double Process::matrixMRSSMTree_udbar_suLsdLdagger(double alphas, double S, double T) const {
   double U = 2*MassSq*MassSq - S - T;
   return  (315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2);
}
