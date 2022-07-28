inline double Process::matrixMSSMTree_uubar_suLsuRdagger(double alphas, double S, double T ) const { // agrees with Philip
   double U = 2*MassSq*MassSq - S - T;
   return (-631.6546816697189*Sqr(alphas)*(pow(MassSq,4) - 1.*T*U))/(S*S);
} 
