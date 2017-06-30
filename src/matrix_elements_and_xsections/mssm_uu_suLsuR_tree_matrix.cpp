/*
 *  agrees with Philip + checked with MadGraph (same as MRSSM)
 */

inline double Process::matrixMSSMTree_uu_suLsuR( double S, double T ) const { 
   double alphaS = pdf->alphasQ( mu_r );
   double U = 2*MassSq*MassSq - S - T;
   return (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2) + 
           (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + U,2);
}
