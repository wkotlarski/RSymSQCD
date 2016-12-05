/*
 *    agrees with Philip + checked with MadGraph (same as MRSSM) 
 */

inline double Process::matrixMSSMTree_ud_suLsdR( double S, double T ) {
	double alphaS = pdf->alphasQ( mu_r );
	double U = 2*MassSq*MassSq - S - T;
   return (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(MassGlu*MassGlu - 1.*T,2); 
} 

