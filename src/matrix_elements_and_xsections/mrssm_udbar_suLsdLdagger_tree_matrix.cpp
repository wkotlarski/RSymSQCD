double Process::matrixMRSSMTree_udbar_suLsdLdagger( double S, double T ) { 
   double alphaS = pdf->alphasQ( mu_r );
   double U = 2*MassSq*MassSq - S - T;
   //std:cout << S << '\n';
   return  (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2);
} 
