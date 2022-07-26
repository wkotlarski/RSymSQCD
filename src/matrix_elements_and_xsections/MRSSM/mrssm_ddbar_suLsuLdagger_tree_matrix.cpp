inline double Process::matrixMRSSMTree_ddbar_suLsuLdagger(double alphas, double S, double T ) const { // agrees with Philip
   const double U = 2*MassSq*MassSq - S - T;
   return  (-631.6546816697189*Sqr(alphas)*(pow(MassSq,4) - T*U))/Sqr(S);
}
