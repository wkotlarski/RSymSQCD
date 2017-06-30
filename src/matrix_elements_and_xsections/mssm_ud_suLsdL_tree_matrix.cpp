inline double Process::matrixMSSMTree_ud_suLsdL( double S, double T ) const { // checked with MadGraph + checked with MadGraph
   double alphaS = pdf->alphasQ( mu_r );
   return (315.82734083485946*(alphaS*alphaS)*(MassGlu*MassGlu)*S)/pow(MassGlu*MassGlu - 1.*T,2);
}

