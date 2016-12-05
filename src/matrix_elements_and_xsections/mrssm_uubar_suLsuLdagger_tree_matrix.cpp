inline double Process::matrixMRSSMTree_uubar_suLsuLdagger(double alphaS, double T, double U, double S) {

      return  (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/(S*S) + 355.3057584392169*(alphaS*alphaS)*pow(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U) - 236.8705056261446*(alphaS*alphaS)*(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T))*(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T))*(-1.*pow(MassSq,4) + T*U) + 355.3057584392169*(alphaS*alphaS)*pow(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U);
}
