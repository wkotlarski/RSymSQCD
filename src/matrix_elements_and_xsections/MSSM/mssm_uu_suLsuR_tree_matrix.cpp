/*
 *  agrees with Philip + checked with MadGraph (same as MRSSM)
 */

inline double Process::matrixMSSMTree_uu_suLsuR(double alphas, double S, double T) const {
   static constexpr double h = 2.*2;
   static constexpr double k = 2.*2*3*3;
   const double U = 2*MassSq*MassSq - S - T;
   return h/k*((315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2) +
           (315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + U,2));
}
