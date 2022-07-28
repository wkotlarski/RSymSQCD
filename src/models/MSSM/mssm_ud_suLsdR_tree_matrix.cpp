/*
 *    agrees with Philip + checked with MadGraph (same as MRSSM) 
 */

inline double Process::matrixMSSMTree_ud_suLsdR(double alphas, double S, double T ) const {
   double U = 2*MassSq*MassSq - S - T;
   return (315.82734083485946*Sqr(alphas)*(-1.*pow(MassSq,4) + T*U))/pow(MassGlu*MassGlu - 1.*T,2); 
}

