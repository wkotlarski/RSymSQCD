inline double Process::matrixSMTree_eebar_ttbar (double S, double T) const {
   const double MB2 = 5*5;
   const double Alfa2 = pow(137., -2);
   return 4.*(8*Alfa2*(pi*pi)*(S*S + 2*(MB2*MB2 - 2*MB2*T + T*(S + T))))/(3.*(S*S));
}