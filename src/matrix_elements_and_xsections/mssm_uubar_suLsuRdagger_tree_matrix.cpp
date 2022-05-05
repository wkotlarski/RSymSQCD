inline double Process::matrixMSSMTree_uubar_suLsuRdagger( double S, double T ) const { // agrees with Philip
   double alphaS = pdf->alphasQ( mu_r );
   double U = 2*MassSq*MassSq - S - T;
   return (-631.6546816697189*(alphaS*alphaS)*(pow(MassSq,4) - 1.*T*U))/(S*S);
} 
