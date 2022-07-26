double Process::matrixMRSSMHard_gu_suLsuRubar(double Alfas, std::array<std::array<double, 4>, 5> const& p) const {
   const double Alfas2 = Sqr(Alfas);
   const double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   const double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   const double k45 = p[3][0]*p[4][0]-p[3][1]*p[4][1]-p[3][2]*p[4][2]-p[3][3]*p[4][3];
   const double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   const double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   const double S12 = 2.*k12;
   const double sqrm1 = Sqr(MassSq);
   const double sqrm2 = Sqr(MassSq);
   const double m1_4 = Sqr(sqrm1);
   const double m2_4 = Sqr(sqrm2);
   const double S35 = sqrm1 + 2*k35;
   const double S45 = sqrm2 + 2*k45;
   const double T = sqrm1 - 2.*k13;
   const double T14 = sqrm2 - 2.*k14;
   const double sqrT = Sqr(T);
   const double sqrS12 = Sqr(S12);
   const double cubeS12 = sqrS12*S12;
   const double sqrMassGlu = Sqr(MassGlu);
   const double sqrT14 = Sqr(T14);
   const double sqrS35 = Sqr(S35);
   const double sqrS45 = Sqr(S45);
   return (8*Alfas*Alfas2*Power3(pi)*(-(pow(S12,-1)*pow(S12 - S45 + T,-1)*(S12*(S12 - S45 + T)*(-S45 + sqrMassGlu)*(sqrm1*(-((S12 - S35 + T14)*(17*S12 - 17*S35 - S45 + 17*T14)) + (15*S12 - 15*S35 - 2*T + 15*T14)*sqrm2 + 2*m2_4) + T*(-17*S12*S35 + S35*S45 - 17*S35*T14 + T*T14 + (S12 - 2*S35 - S45 + T)*sqrm2 + 17*sqrS35)) + S12*(-9*(S12 - S45 + T)*((S12 - S35 - S45 + T + T14)*(-2*S12*S35*T - S12*S45*T + S35*S45*T + S12*T*T14 - S35*T*T14 - S45*T*T14 + S12*S35*sqrMassGlu + S12*S45*sqrMassGlu - S12*T14*sqrMassGlu + S45*T14*sqrMassGlu + T*T14*sqrMassGlu + sqrm2*(2*S35*T + (S12 - S45 + T)*sqrMassGlu) - sqrm1*((S12 - S35 + T14)*(S12 - S35 - S45 + T14) + 2*m2_4 + sqrm2*(-3*(S12 - S35 + T14) + 2*sqrMassGlu)) + T*sqrS12 - sqrMassGlu*sqrS12 + T*sqrS35) + (S12 - S35 + T14)*(-2*S12*S35*T - S12*S45*T + S35*S45*T + S12*T*T14 - S35*T*T14 - S45*T*T14 + S12*S35*sqrMassGlu + S12*S45*sqrMassGlu - S12*T*sqrMassGlu - S35*T*sqrMassGlu - S12*T14*sqrMassGlu + S45*T14*sqrMassGlu + sqrm2*(2*(-S12 + S35 + S45 - T)*T + (S12 - S45 + T)*sqrMassGlu) - sqrm1*((-3*S12 + 3*S35 + 2*S45 - 4*T - 3*T14)*sqrm2 + 2*m2_4 + (S12 - S35 + T14)*(S12 - S35 - S45 + T + T14 + sqrMassGlu)) + T*sqrS12 - sqrMassGlu*sqrS12 + T*sqrS35 + S12*sqrT - S35*sqrT - T14*sqrT)) + (-T + sqrMassGlu)*(2*(S12 - S35 - S45)*(S12 - S35 - S45 + T + T14)*(S12*S35 + S12*S45 - S12*T + S35*T - S12*T14 + S45*T14 + sqrm1*(S12 - S35 + T14 - 2*sqrm2) + (S12 - S45 + T)*sqrm2 - sqrS12) - 7*(S12 - S35 + T14)*(2*S12*S35*S45 - 3*S12*S35*T - 2*S12*S45*T + S35*S45*T - S12*S35*T14 - 2*S12*S45*T14 + S35*S45*T14 - (S12 - S45 + T)*(S12 - S35 - S45 + T)*sqrm2 + 2*m1_4*sqrm2 - 2*S35*sqrS12 - 2*S45*sqrS12 + 2*T*sqrS12 + T14*sqrS12 + cubeS12 + sqrm1*(-(S35*S45) + S35*T + S12*(2*S35 + S45 - T - T14) + S35*T14 - S45*T14 - T*T14 + (S12 - 2*S35 - S45 + T)*sqrm2 - sqrS12 - sqrS35) + S12*sqrS35 + T*sqrS35 + S12*sqrS45 + T14*sqrS45 + S12*sqrT - S35*sqrT - T14*sqrT))) + (S12 - S45 + T)*(-T + sqrMassGlu)*(7*(S12 - S35 - S45 + T + T14)*(-(sqrm1*(-T14 + sqrm2)*(-S12 + S35 - T14 + 2*sqrm2)) + sqrm2*(-(S12*(S35 + S45 - T)) + S35*(S45 + T) + sqrS12) - (S12 - S35)*(-(S12*(S35 + S45 - T14)) - (S45 + T)*T14 + sqrS12)) - 2*(S12 - S35 + T14)*(S12*S35*S45 - S12*S35*T + S35*S45*T - S12*S35*T14 - S12*S45*T14 + S35*S45*T14 - S12*T*T14 + S35*T*T14 + S45*T*T14 + 2*m1_4*sqrm2 - 2*S35*sqrS12 - S45*sqrS12 + T*sqrS12 + T14*sqrS12 - sqrm2*(-(S12*(S35 + S45)) + (S45 - T)*T + S35*(S45 + T) + sqrS12) + cubeS12 + S12*sqrS35 + sqrm1*(-(S35*S45) - S35*T14 - S45*T14 + S12*(S45 + T14) + (-2*S12 + S35 + S45 - 3*T - 3*T14)*sqrm2 + 2*m2_4 + sqrT14))))*pow(S12 - S35 + T14,-1)*pow(S12 - S35 - S45 + T + T14,-1)*pow(-S45 + sqrMassGlu,-1)*pow(-T + sqrMassGlu,-2)) + pow(S45 - sqrMassGlu,-1)*(16*((S12 - S35)*S45 + sqrm1*sqrm2)*pow(S12,-1)*pow(S45 - sqrMassGlu,-1) + pow(T - sqrMassGlu,-1)*(pow(-S12 + S45 - T,-1)*(-2*(S12 - S35 - S45)*(-(S12*S35) - S12*S45 + S12*T - S35*T + S12*T14 - S45*T14 - (S12 - S45 + T)*sqrm2 + sqrm1*(-S12 + S35 - T14 + 2*sqrm2) + sqrS12)*pow(S12 - S35 + T14,-1) - 7*(2*S12*S35*S45 - 3*S12*S35*T - 2*S12*S45*T + S35*S45*T - S12*S35*T14 - 2*S12*S45*T14 + S35*S45*T14 - (S12 - S45 + T)*(S12 - S35 - S45 + T)*sqrm2 + 2*m1_4*sqrm2 - 2*S35*sqrS12 - 2*S45*sqrS12 + 2*T*sqrS12 + T14*sqrS12 + cubeS12 + sqrm1*(-(S35*S45) + S35*T + S12*(2*S35 + S45 - T - T14) + S35*T14 - S45*T14 - T*T14 + (S12 - 2*S35 - S45 + T)*sqrm2 - sqrS12 - sqrS35) + S12*sqrS35 + T*sqrS35 + S12*sqrS45 + T14*sqrS45 + S12*sqrT - S35*sqrT - T14*sqrT)*pow(S12 - S35 - S45 + T + T14,-1)) - 9*((-2*S12*S35*T - S12*S45*T + S35*S45*T + S12*T*T14 - S35*T*T14 - S45*T*T14 + S12*S35*sqrMassGlu + S12*S45*sqrMassGlu - S12*T14*sqrMassGlu + S45*T14*sqrMassGlu + T*T14*sqrMassGlu + sqrm2*(2*S35*T + (S12 - S45 + T)*sqrMassGlu) - sqrm1*((S12 - S35 + T14)*(S12 - S35 - S45 + T14) + 2*m2_4 + sqrm2*(-3*(S12 - S35 + T14) + 2*sqrMassGlu)) + T*sqrS12 - sqrMassGlu*sqrS12 + T*sqrS35)*pow(S12 - S35 + T14,-1) - (2*S12*S35*T + S12*S45*T - S35*S45*T - S12*T*T14 + S35*T*T14 + S45*T*T14 - S12*S35*sqrMassGlu - S12*S45*sqrMassGlu + S12*T*sqrMassGlu + S35*T*sqrMassGlu + S12*T14*sqrMassGlu - S45*T14*sqrMassGlu - sqrm2*(2*(-S12 + S35 + S45 - T)*T + (S12 - S45 + T)*sqrMassGlu) + sqrm1*((-3*S12 + 3*S35 + 2*S45 - 4*T - 3*T14)*sqrm2 + 2*m2_4 + (S12 - S35 + T14)*(S12 - S35 - S45 + T + T14 + sqrMassGlu)) - T*sqrS12 + sqrMassGlu*sqrS12 - T*sqrS35 - S12*sqrT + S35*sqrT + T14*sqrT)*pow(S12 - S35 - S45 + T + T14,-1))*pow(T - sqrMassGlu,-1)) + pow(S12,-1)*((7*(sqrm1*(-T14 + sqrm2)*(-S12 + S35 - T14 + 2*sqrm2) - sqrm2*(-(S12*(S35 + S45 - T)) + S35*(S45 + T) + sqrS12) + (S12 - S35)*(-(S12*(S35 + S45 - T14)) - (S45 + T)*T14 + sqrS12))*pow(S12 - S35 + T14,-1) + 2*(S12*S35*S45 - S12*S35*T + S35*S45*T - S12*S35*T14 - S12*S45*T14 + S35*S45*T14 - S12*T*T14 + S35*T*T14 + S45*T*T14 + 2*m1_4*sqrm2 - 2*S35*sqrS12 - S45*sqrS12 + T*sqrS12 + T14*sqrS12 - sqrm2*(-(S12*(S35 + S45)) + (S45 - T)*T + S35*(S45 + T) + sqrS12) + cubeS12 + S12*sqrS35 + sqrm1*(-(S35*S45) - S35*T14 - S45*T14 + S12*(S45 + T14) + (-2*S12 + S35 + S45 - 3*T - 3*T14)*sqrm2 + 2*m2_4 + sqrT14))*pow(S12 - S35 - S45 + T + T14,-1))*pow(T - sqrMassGlu,-1) + 2*pow(S45 - sqrMassGlu,-1)*((-2*m1_4*sqrm2 - T*(S45*(S35 + T14) + (-S45 + T)*sqrm2) + sqrm1*(S45*(-S12 + S35 + T14) + (3*S12 - S45 + 3*T)*sqrm2) + S45*sqrS12 - S12*(S35*S45 - S45*T + T*sqrm2 + sqrS45))*pow(S12 - S45 + T,-1) + 9*(2*m1_4*sqrm2 + T*(S45*(-S12 + S35 + T14) + (S12 - S45 + T)*sqrm2) - sqrm1*(S45*(-S12 + S35 + T14) + (S12 - S45 + 3*T)*sqrm2) + S12*(S12 - S35 + S45)*sqrMassGlu)*pow(T - sqrMassGlu,-1))) - 2*pow(S45 - sqrMassGlu,-1)*(16*sqrm1*(-(S45*T14) - (S12 - S45 + T)*sqrm2 + sqrm1*sqrm2)*pow(S12 - S45 + T,-2) - 18*(T*(S12 - S45 + T)*sqrm2 + T*(-S12 + S35 + T14)*sqrMassGlu + sqrm1*(sqrm2*(S45 - T - 2*sqrMassGlu) + (S12 - S35 + T14)*sqrMassGlu) + S12*Sqr(sqrMassGlu))*pow(-T + sqrMassGlu,-2) + 9*(-(T*(S12 - S45 + T)*sqrm2) + 2*m1_4*sqrm2 + sqrm1*(-2*S45*T14 - (S12 - S35 + T14)*sqrMassGlu + sqrm2*(-S12 + S45 - T + 2*sqrMassGlu)) + sqrMassGlu*(-(S12*(S35 + S45 - T)) - T*(S35 + T14) + sqrS12))*pow(S12 - S45 + T,-1)*pow(-T + sqrMassGlu,-1))) - 32*(sqrm2*(-(S35*T) + sqrm1*(-S12 + S35 - T14 + sqrm2))*pow(S12 - S35 + T14,-2)*pow(-T + sqrMassGlu,-2) + sqrm1*(-(S45*T14) - (S12 - S45 + T)*sqrm2 + sqrm1*sqrm2)*pow(S12 - S45 + T,-2)*pow(-T14 + sqrMassGlu,-2)) + pow(S12 - S35 - S45 + T + T14,-1)*((sqrm1*((S12 - S35 + T14)*(S12 - S35 - S45 + T14) + (S12 - S35 + 2*T + T14)*sqrm2 - 2*m2_4) - T*(-(S12*S35) + S35*S45 - S35*T14 + T*T14 + (S12 - 2*S35 - S45 + T)*sqrm2 + sqrS35))*pow(S12 - S35 + T14,-1)*pow(-T + sqrMassGlu,-2) + 2*((S12 - S45 + T)*(S12 - S35 - S45 + T)*sqrm2 - 2*m1_4*sqrm2 + sqrm1*(-(T14*(S12 - S35 - 2*S45 + T14)) + (S12 - S45 + T + 2*T14)*sqrm2) - T14*(-(S12*S45) + S35*S45 - S45*T + T*T14 + sqrS45))*pow(S12 - S45 + T,-1)*pow(-T14 + sqrMassGlu,-2)) - pow(S35 - sqrMassGlu,-1)*(pow(S12,-1)*((-7*(2*m1_4*sqrm2 + T*(S45*(S35 + T14) + (-S45 + T)*sqrm2) + (-S35 - 2*S45 + T)*sqrS12 - sqrm1*(-(S12*(S35 + S45 - T14)) + S45*(S35 + T14) + (S12 - S45 + 3*T)*sqrm2 + sqrS12) + cubeS12 + S12*(S35*(S45 - T) - S45*T + T*(-T14 + sqrm2) + sqrS45))*pow(S12 - S45 + T,-1) - 2*(S35*S45*T + S35*S45*T14 + S35*T*T14 + S45*T*T14 - S35*S45*sqrm2 - S35*T*sqrm2 - S45*T*sqrm2 + 2*m1_4*sqrm2 + (-S35 - 2*S45 + T + T14)*sqrS12 + cubeS12 + S12*(S35*(S45 - T) - S45*T - S45*T14 - T*T14 + (S35 + T)*sqrm2 + sqrS45) + sqrm2*sqrT + sqrm1*(-(S35*S45) + S12*(S35 + S45) - S35*T14 - S45*T14 + (-2*S12 + S35 + S45 - 3*T - 3*T14)*sqrm2 + 2*m2_4 - sqrS12 + sqrT14))*pow(S12 - S35 - S45 + T + T14,-1))*pow(T14 - sqrMassGlu,-1) - 2*pow(S35 - sqrMassGlu,-1)*((-(sqrm1*(T14*(S12 - S35 + T14) + (-3*S12 + S35 - 3*T14)*sqrm2 + 2*m2_4)) + S35*(-(S12*(S35 + S45 - T14)) - (S45 + T)*T14 + (-S12 + S45 + T)*sqrm2 + sqrS12))*pow(S12 - S35 + T14,-1) + 9*(S35*(S12 - S45 - T)*(-T14 + sqrm2) + sqrm1*(-T14 + sqrm2)*(-S12 + S35 - T14 + 2*sqrm2) + S12*(S12 + S35 - S45)*sqrMassGlu)*pow(T14 - sqrMassGlu,-1))) + 16*(S35*(S12 - S45) + sqrm1*sqrm2)*pow(S12,-1)*pow(-S35 + sqrMassGlu,-1) + 2*pow(S35 - sqrMassGlu,-1)*(16*sqrm2*(-(S35*T) + sqrm1*(-S12 + S35 - T14 + sqrm2))*pow(S12 - S35 + T14,-2) - 18*(sqrm1*(T14*(S12 - S35 + T14) + sqrm2*(S35 - T14 - 2*sqrMassGlu)) + sqrMassGlu*((-S12 + S45 + T)*T14 + (S12 - S45 + T)*sqrm2 + S12*sqrMassGlu))*pow(-T14 + sqrMassGlu,-2) + 9*(-(sqrm2*(2*S35*T + (S12 - S45 + T)*sqrMassGlu)) + sqrm1*(-(T14*(S12 - S35 + T14)) + 2*m2_4 + sqrm2*(-S12 + S35 - T14 + 2*sqrMassGlu)) + sqrMassGlu*(-(S12*(S35 + S45 - T14)) - (S45 + T)*T14 + sqrS12))*pow(S12 - S35 + T14,-1)*pow(-T14 + sqrMassGlu,-1)) + pow(T14 - sqrMassGlu,-1)*(pow(-S12 + S45 - T,-1)*(2*(S12 - S35 - S45)*(-(S12*S35) - S12*S45 + S12*T - S35*T + S12*T14 - S45*T14 - (S12 - S45 + T)*sqrm2 + sqrm1*(-S12 + S35 - T14 + 2*sqrm2) + sqrS12)*pow(S12 - S35 + T14,-1) + 9*(-(S12*S35*T14) - 2*S12*S45*T14 + S35*S45*T14 + S12*T*T14 - S35*T*T14 - S45*T*T14 - (S12 - S45 + T)*(S12 - S35 - S45 + T)*sqrm2 - 2*m1_4*sqrm2 + S12*S35*sqrMassGlu + S12*S45*sqrMassGlu - S12*T*sqrMassGlu + S35*T*sqrMassGlu + T*T14*sqrMassGlu + sqrm1*(2*S45*T14 + sqrm2*(3*(S12 - S45 + T) - 2*sqrMassGlu) + (S12 - S35 + T14)*sqrMassGlu) + T14*sqrS12 - sqrMassGlu*sqrS12 + T14*sqrS45)*pow(-T14 + sqrMassGlu,-1)) - pow(S12 - S35 - S45 + T + T14,-1)*(7*(2*S12*S35*S45 - 2*S12*S35*T - S12*S45*T + S35*S45*T - 2*S12*S35*T14 - 3*S12*S45*T14 + S35*S45*T14 + sqrm1*(S12 - S35 - S45 + T14 + sqrm2)*(-S12 + S35 - T14 + 2*sqrm2) - 2*S35*sqrS12 - 2*S45*sqrS12 + T*sqrS12 + 2*T14*sqrS12 - sqrm2*(S35*(S45 + T) + (S45 - T)*(S45 - T14) + S12*(-S35 - 2*S45 + T + T14) + sqrS12) + cubeS12 + S12*sqrS35 + T*sqrS35 + S12*sqrS45 + T14*sqrS45 + S12*sqrT14 - S45*sqrT14 - T*sqrT14)*pow(S12 - S35 + T14,-1) - 9*(S12*S35*T14 + 2*S12*S45*T14 - S35*S45*T14 - S12*T*T14 + S35*T*T14 + S45*T*T14 + 2*m1_4*sqrm2 - S12*S35*sqrMassGlu - S12*S45*sqrMassGlu + S12*T*sqrMassGlu - S35*T*sqrMassGlu + S12*T14*sqrMassGlu + S45*T14*sqrMassGlu + (S12 - S45 + T)*sqrm2*(S12 - S35 - S45 + T + T14 + sqrMassGlu) - sqrm1*(2*(-S12 + S35 + S45 - T14)*T14 + (3*S12 - 2*S35 - 3*S45 + 3*T + 4*T14)*sqrm2 + (S12 - S35 + T14)*sqrMassGlu) - T14*sqrS12 + sqrMassGlu*sqrS12 - T14*sqrS45 - S12*sqrT14 + S45*sqrT14 + T*sqrT14)*pow(-T14 + sqrMassGlu,-1)))) + pow(T14 - sqrMassGlu,-1)*(pow(S12,-1)*(7*(2*m1_4*sqrm2 + T*(S45*(S35 + T14) + (-S45 + T)*sqrm2) + (-S35 - 2*S45 + T)*sqrS12 - sqrm1*(-(S12*(S35 + S45 - T14)) + S45*(S35 + T14) + (S12 - S45 + 3*T)*sqrm2 + sqrS12) + cubeS12 + S12*(S35*(S45 - T) - S45*T + T*(-T14 + sqrm2) + sqrS45))*pow(S12 - S45 + T,-1) + 2*(S35*S45*T + S35*S45*T14 + S35*T*T14 + S45*T*T14 - S35*S45*sqrm2 - S35*T*sqrm2 - S45*T*sqrm2 + 2*m1_4*sqrm2 + (-S35 - 2*S45 + T + T14)*sqrS12 + cubeS12 + S12*(S35*(S45 - T) - S45*T - S45*T14 - T*T14 + (S35 + T)*sqrm2 + sqrS45) + sqrm2*sqrT + sqrm1*(-(S35*S45) + S12*(S35 + S45) - S35*T14 - S45*T14 + (-2*S12 + S35 + S45 - 3*T - 3*T14)*sqrm2 + 2*m2_4 - sqrS12 + sqrT14))*pow(S12 - S35 - S45 + T + T14,-1))*pow(S35 - sqrMassGlu,-1) + 16*(S45*T14 + (S12 - S45 + T)*sqrm2 - sqrm1*sqrm2)*pow(S12 - S35 - S45 + T + T14,-1)*pow(T14 - sqrMassGlu,-1) + pow(S35 - sqrMassGlu,-1)*(pow(-S12 + S45 - T,-1)*(-2*(S12 - S35 - S45)*(-(S12*S35) - S12*S45 + S12*T - S35*T + S12*T14 - S45*T14 - (S12 - S45 + T)*sqrm2 + sqrm1*(-S12 + S35 - T14 + 2*sqrm2) + sqrS12)*pow(S12 - S35 + T14,-1) - 9*(-(S12*S35*T14) - 2*S12*S45*T14 + S35*S45*T14 + S12*T*T14 - S35*T*T14 - S45*T*T14 - (S12 - S45 + T)*(S12 - S35 - S45 + T)*sqrm2 - 2*m1_4*sqrm2 + S12*S35*sqrMassGlu + S12*S45*sqrMassGlu - S12*T*sqrMassGlu + S35*T*sqrMassGlu + T*T14*sqrMassGlu + sqrm1*(2*S45*T14 + sqrm2*(3*(S12 - S45 + T) - 2*sqrMassGlu) + (S12 - S35 + T14)*sqrMassGlu) + T14*sqrS12 - sqrMassGlu*sqrS12 + T14*sqrS45)*pow(-T14 + sqrMassGlu,-1)) + pow(S12 - S35 - S45 + T + T14,-1)*(7*(2*S12*S35*S45 - 2*S12*S35*T - S12*S45*T + S35*S45*T - 2*S12*S35*T14 - 3*S12*S45*T14 + S35*S45*T14 + sqrm1*(S12 - S35 - S45 + T14 + sqrm2)*(-S12 + S35 - T14 + 2*sqrm2) - 2*S35*sqrS12 - 2*S45*sqrS12 + T*sqrS12 + 2*T14*sqrS12 - sqrm2*(S35*(S45 + T) + (S45 - T)*(S45 - T14) + S12*(-S35 - 2*S45 + T + T14) + sqrS12) + cubeS12 + S12*sqrS35 + T*sqrS35 + S12*sqrS45 + T14*sqrS45 + S12*sqrT14 - S45*sqrT14 - T*sqrT14)*pow(S12 - S35 + T14,-1) - 9*(S12*S35*T14 + 2*S12*S45*T14 - S35*S45*T14 - S12*T*T14 + S35*T*T14 + S45*T*T14 + 2*m1_4*sqrm2 - S12*S35*sqrMassGlu - S12*S45*sqrMassGlu + S12*T*sqrMassGlu - S35*T*sqrMassGlu + S12*T14*sqrMassGlu + S45*T14*sqrMassGlu + (S12 - S45 + T)*sqrm2*(S12 - S35 - S45 + T + T14 + sqrMassGlu) - sqrm1*(2*(-S12 + S35 + S45 - T14)*T14 + (3*S12 - 2*S35 - 3*S45 + 3*T + 4*T14)*sqrm2 + (S12 - S35 + T14)*sqrMassGlu) - T14*sqrS12 + sqrMassGlu*sqrS12 - T14*sqrS45 - S12*sqrT14 + S45*sqrT14 + T*sqrT14)*pow(-T14 + sqrMassGlu,-1))))))/9.;
}
