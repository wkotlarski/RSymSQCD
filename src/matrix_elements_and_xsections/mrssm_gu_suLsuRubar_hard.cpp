double Process::matrixMRSSMHard_gu_suLsuRubar( std::vector< double* >& p ) {
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow(Alfas, 2);
   double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   double k45 = p[3][0]*p[4][0]-p[3][1]*p[4][1]-p[3][2]*p[4][2]-p[3][3]*p[4][3];
   double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   double S12 = 2.*k12;
   double S35 = m1*m1 + 2*k35;
   double S45 = m2*m2 + 2*k45;
   double T = m1*m1 - 2.*k13;
   double T14 = m2*m2 - 2.*k14;
   return (8*Alfas*Alfas2*pow(pi,3)*(-(pow(S12,-1)*pow(S12 - S45 + T,-1)*(S12*(S12 - S45 + T)*(-S45 + pow(MassGlu,2))*(pow(m1,2)*(-((S12 - S35 + T14)*(17*S12 - 17*S35 - S45 + 17*T14)) + (15*S12 - 15*S35 - 2*T + 15*T14)*pow(m2,2) + 2*pow(m2,4)) + T*(-17*S12*S35 + S35*S45 - 17*S35*T14 + T*T14 + (S12 - 2*S35 - S45 + T)*pow(m2,2) + 17*pow(S35,2))) + S12*(-9*(S12 - S45 + T)*((S12 - S35 - S45 + T + T14)*(-2*S12*S35*T - S12*S45*T + S35*S45*T + S12*T*T14 - S35*T*T14 - S45*T*T14 + S12*S35*pow(MassGlu,2) + S12*S45*pow(MassGlu,2) - S12*T14*pow(MassGlu,2) + S45*T14*pow(MassGlu,2) + T*T14*pow(MassGlu,2) + pow(m2,2)*(2*S35*T + (S12 - S45 + T)*pow(MassGlu,2)) - pow(m1,2)*((S12 - S35 + T14)*(S12 - S35 - S45 + T14) + 2*pow(m2,4) + pow(m2,2)*(-3*(S12 - S35 + T14) + 2*pow(MassGlu,2))) + T*pow(S12,2) - pow(MassGlu,2)*pow(S12,2) + T*pow(S35,2)) + (S12 - S35 + T14)*(-2*S12*S35*T - S12*S45*T + S35*S45*T + S12*T*T14 - S35*T*T14 - S45*T*T14 + S12*S35*pow(MassGlu,2) + S12*S45*pow(MassGlu,2) - S12*T*pow(MassGlu,2) - S35*T*pow(MassGlu,2) - S12*T14*pow(MassGlu,2) + S45*T14*pow(MassGlu,2) + pow(m2,2)*(2*(-S12 + S35 + S45 - T)*T + (S12 - S45 + T)*pow(MassGlu,2)) - pow(m1,2)*((-3*S12 + 3*S35 + 2*S45 - 4*T - 3*T14)*pow(m2,2) + 2*pow(m2,4) + (S12 - S35 + T14)*(S12 - S35 - S45 + T + T14 + pow(MassGlu,2))) + T*pow(S12,2) - pow(MassGlu,2)*pow(S12,2) + T*pow(S35,2) + S12*pow(T,2) - S35*pow(T,2) - T14*pow(T,2))) + (-T + pow(MassGlu,2))*(2*(S12 - S35 - S45)*(S12 - S35 - S45 + T + T14)*(S12*S35 + S12*S45 - S12*T + S35*T - S12*T14 + S45*T14 + pow(m1,2)*(S12 - S35 + T14 - 2*pow(m2,2)) + (S12 - S45 + T)*pow(m2,2) - pow(S12,2)) - 7*(S12 - S35 + T14)*(2*S12*S35*S45 - 3*S12*S35*T - 2*S12*S45*T + S35*S45*T - S12*S35*T14 - 2*S12*S45*T14 + S35*S45*T14 - (S12 - S45 + T)*(S12 - S35 - S45 + T)*pow(m2,2) + 2*pow(m1,4)*pow(m2,2) - 2*S35*pow(S12,2) - 2*S45*pow(S12,2) + 2*T*pow(S12,2) + T14*pow(S12,2) + pow(S12,3) + pow(m1,2)*(-(S35*S45) + S35*T + S12*(2*S35 + S45 - T - T14) + S35*T14 - S45*T14 - T*T14 + (S12 - 2*S35 - S45 + T)*pow(m2,2) - pow(S12,2) - pow(S35,2)) + S12*pow(S35,2) + T*pow(S35,2) + S12*pow(S45,2) + T14*pow(S45,2) + S12*pow(T,2) - S35*pow(T,2) - T14*pow(T,2)))) + (S12 - S45 + T)*(-T + pow(MassGlu,2))*(7*(S12 - S35 - S45 + T + T14)*(-(pow(m1,2)*(-T14 + pow(m2,2))*(-S12 + S35 - T14 + 2*pow(m2,2))) + pow(m2,2)*(-(S12*(S35 + S45 - T)) + S35*(S45 + T) + pow(S12,2)) - (S12 - S35)*(-(S12*(S35 + S45 - T14)) - (S45 + T)*T14 + pow(S12,2))) - 2*(S12 - S35 + T14)*(S12*S35*S45 - S12*S35*T + S35*S45*T - S12*S35*T14 - S12*S45*T14 + S35*S45*T14 - S12*T*T14 + S35*T*T14 + S45*T*T14 + 2*pow(m1,4)*pow(m2,2) - 2*S35*pow(S12,2) - S45*pow(S12,2) + T*pow(S12,2) + T14*pow(S12,2) - pow(m2,2)*(-(S12*(S35 + S45)) + (S45 - T)*T + S35*(S45 + T) + pow(S12,2)) + pow(S12,3) + S12*pow(S35,2) + pow(m1,2)*(-(S35*S45) - S35*T14 - S45*T14 + S12*(S45 + T14) + (-2*S12 + S35 + S45 - 3*T - 3*T14)*pow(m2,2) + 2*pow(m2,4) + pow(T14,2)))))*pow(S12 - S35 + T14,-1)*pow(S12 - S35 - S45 + T + T14,-1)*pow(-S45 + pow(MassGlu,2),-1)*pow(-T + pow(MassGlu,2),-2)) + pow(S45 - pow(MassGlu,2),-1)*(16*((S12 - S35)*S45 + pow(m1,2)*pow(m2,2))*pow(S12,-1)*pow(S45 - pow(MassGlu,2),-1) + pow(T - pow(MassGlu,2),-1)*(pow(-S12 + S45 - T,-1)*(-2*(S12 - S35 - S45)*(-(S12*S35) - S12*S45 + S12*T - S35*T + S12*T14 - S45*T14 - (S12 - S45 + T)*pow(m2,2) + pow(m1,2)*(-S12 + S35 - T14 + 2*pow(m2,2)) + pow(S12,2))*pow(S12 - S35 + T14,-1) - 7*(2*S12*S35*S45 - 3*S12*S35*T - 2*S12*S45*T + S35*S45*T - S12*S35*T14 - 2*S12*S45*T14 + S35*S45*T14 - (S12 - S45 + T)*(S12 - S35 - S45 + T)*pow(m2,2) + 2*pow(m1,4)*pow(m2,2) - 2*S35*pow(S12,2) - 2*S45*pow(S12,2) + 2*T*pow(S12,2) + T14*pow(S12,2) + pow(S12,3) + pow(m1,2)*(-(S35*S45) + S35*T + S12*(2*S35 + S45 - T - T14) + S35*T14 - S45*T14 - T*T14 + (S12 - 2*S35 - S45 + T)*pow(m2,2) - pow(S12,2) - pow(S35,2)) + S12*pow(S35,2) + T*pow(S35,2) + S12*pow(S45,2) + T14*pow(S45,2) + S12*pow(T,2) - S35*pow(T,2) - T14*pow(T,2))*pow(S12 - S35 - S45 + T + T14,-1)) - 9*((-2*S12*S35*T - S12*S45*T + S35*S45*T + S12*T*T14 - S35*T*T14 - S45*T*T14 + S12*S35*pow(MassGlu,2) + S12*S45*pow(MassGlu,2) - S12*T14*pow(MassGlu,2) + S45*T14*pow(MassGlu,2) + T*T14*pow(MassGlu,2) + pow(m2,2)*(2*S35*T + (S12 - S45 + T)*pow(MassGlu,2)) - pow(m1,2)*((S12 - S35 + T14)*(S12 - S35 - S45 + T14) + 2*pow(m2,4) + pow(m2,2)*(-3*(S12 - S35 + T14) + 2*pow(MassGlu,2))) + T*pow(S12,2) - pow(MassGlu,2)*pow(S12,2) + T*pow(S35,2))*pow(S12 - S35 + T14,-1) - (2*S12*S35*T + S12*S45*T - S35*S45*T - S12*T*T14 + S35*T*T14 + S45*T*T14 - S12*S35*pow(MassGlu,2) - S12*S45*pow(MassGlu,2) + S12*T*pow(MassGlu,2) + S35*T*pow(MassGlu,2) + S12*T14*pow(MassGlu,2) - S45*T14*pow(MassGlu,2) - pow(m2,2)*(2*(-S12 + S35 + S45 - T)*T + (S12 - S45 + T)*pow(MassGlu,2)) + pow(m1,2)*((-3*S12 + 3*S35 + 2*S45 - 4*T - 3*T14)*pow(m2,2) + 2*pow(m2,4) + (S12 - S35 + T14)*(S12 - S35 - S45 + T + T14 + pow(MassGlu,2))) - T*pow(S12,2) + pow(MassGlu,2)*pow(S12,2) - T*pow(S35,2) - S12*pow(T,2) + S35*pow(T,2) + T14*pow(T,2))*pow(S12 - S35 - S45 + T + T14,-1))*pow(T - pow(MassGlu,2),-1)) + pow(S12,-1)*((7*(pow(m1,2)*(-T14 + pow(m2,2))*(-S12 + S35 - T14 + 2*pow(m2,2)) - pow(m2,2)*(-(S12*(S35 + S45 - T)) + S35*(S45 + T) + pow(S12,2)) + (S12 - S35)*(-(S12*(S35 + S45 - T14)) - (S45 + T)*T14 + pow(S12,2)))*pow(S12 - S35 + T14,-1) + 2*(S12*S35*S45 - S12*S35*T + S35*S45*T - S12*S35*T14 - S12*S45*T14 + S35*S45*T14 - S12*T*T14 + S35*T*T14 + S45*T*T14 + 2*pow(m1,4)*pow(m2,2) - 2*S35*pow(S12,2) - S45*pow(S12,2) + T*pow(S12,2) + T14*pow(S12,2) - pow(m2,2)*(-(S12*(S35 + S45)) + (S45 - T)*T + S35*(S45 + T) + pow(S12,2)) + pow(S12,3) + S12*pow(S35,2) + pow(m1,2)*(-(S35*S45) - S35*T14 - S45*T14 + S12*(S45 + T14) + (-2*S12 + S35 + S45 - 3*T - 3*T14)*pow(m2,2) + 2*pow(m2,4) + pow(T14,2)))*pow(S12 - S35 - S45 + T + T14,-1))*pow(T - pow(MassGlu,2),-1) + 2*pow(S45 - pow(MassGlu,2),-1)*((-2*pow(m1,4)*pow(m2,2) - T*(S45*(S35 + T14) + (-S45 + T)*pow(m2,2)) + pow(m1,2)*(S45*(-S12 + S35 + T14) + (3*S12 - S45 + 3*T)*pow(m2,2)) + S45*pow(S12,2) - S12*(S35*S45 - S45*T + T*pow(m2,2) + pow(S45,2)))*pow(S12 - S45 + T,-1) + 9*(2*pow(m1,4)*pow(m2,2) + T*(S45*(-S12 + S35 + T14) + (S12 - S45 + T)*pow(m2,2)) - pow(m1,2)*(S45*(-S12 + S35 + T14) + (S12 - S45 + 3*T)*pow(m2,2)) + S12*(S12 - S35 + S45)*pow(MassGlu,2))*pow(T - pow(MassGlu,2),-1))) - 2*pow(S45 - pow(MassGlu,2),-1)*(16*pow(m1,2)*(-(S45*T14) - (S12 - S45 + T)*pow(m2,2) + pow(m1,2)*pow(m2,2))*pow(S12 - S45 + T,-2) - 18*(T*(S12 - S45 + T)*pow(m2,2) + T*(-S12 + S35 + T14)*pow(MassGlu,2) + pow(m1,2)*(pow(m2,2)*(S45 - T - 2*pow(MassGlu,2)) + (S12 - S35 + T14)*pow(MassGlu,2)) + S12*pow(MassGlu,4))*pow(-T + pow(MassGlu,2),-2) + 9*(-(T*(S12 - S45 + T)*pow(m2,2)) + 2*pow(m1,4)*pow(m2,2) + pow(m1,2)*(-2*S45*T14 - (S12 - S35 + T14)*pow(MassGlu,2) + pow(m2,2)*(-S12 + S45 - T + 2*pow(MassGlu,2))) + pow(MassGlu,2)*(-(S12*(S35 + S45 - T)) - T*(S35 + T14) + pow(S12,2)))*pow(S12 - S45 + T,-1)*pow(-T + pow(MassGlu,2),-1))) - 32*(pow(m2,2)*(-(S35*T) + pow(m1,2)*(-S12 + S35 - T14 + pow(m2,2)))*pow(S12 - S35 + T14,-2)*pow(-T + pow(MassGlu,2),-2) + pow(m1,2)*(-(S45*T14) - (S12 - S45 + T)*pow(m2,2) + pow(m1,2)*pow(m2,2))*pow(S12 - S45 + T,-2)*pow(-T14 + pow(MassGlu,2),-2)) + pow(S12 - S35 - S45 + T + T14,-1)*((pow(m1,2)*((S12 - S35 + T14)*(S12 - S35 - S45 + T14) + (S12 - S35 + 2*T + T14)*pow(m2,2) - 2*pow(m2,4)) - T*(-(S12*S35) + S35*S45 - S35*T14 + T*T14 + (S12 - 2*S35 - S45 + T)*pow(m2,2) + pow(S35,2)))*pow(S12 - S35 + T14,-1)*pow(-T + pow(MassGlu,2),-2) + 2*((S12 - S45 + T)*(S12 - S35 - S45 + T)*pow(m2,2) - 2*pow(m1,4)*pow(m2,2) + pow(m1,2)*(-(T14*(S12 - S35 - 2*S45 + T14)) + (S12 - S45 + T + 2*T14)*pow(m2,2)) - T14*(-(S12*S45) + S35*S45 - S45*T + T*T14 + pow(S45,2)))*pow(S12 - S45 + T,-1)*pow(-T14 + pow(MassGlu,2),-2)) - pow(S35 - pow(MassGlu,2),-1)*(pow(S12,-1)*((-7*(2*pow(m1,4)*pow(m2,2) + T*(S45*(S35 + T14) + (-S45 + T)*pow(m2,2)) + (-S35 - 2*S45 + T)*pow(S12,2) - pow(m1,2)*(-(S12*(S35 + S45 - T14)) + S45*(S35 + T14) + (S12 - S45 + 3*T)*pow(m2,2) + pow(S12,2)) + pow(S12,3) + S12*(S35*(S45 - T) - S45*T + T*(-T14 + pow(m2,2)) + pow(S45,2)))*pow(S12 - S45 + T,-1) - 2*(S35*S45*T + S35*S45*T14 + S35*T*T14 + S45*T*T14 - S35*S45*pow(m2,2) - S35*T*pow(m2,2) - S45*T*pow(m2,2) + 2*pow(m1,4)*pow(m2,2) + (-S35 - 2*S45 + T + T14)*pow(S12,2) + pow(S12,3) + S12*(S35*(S45 - T) - S45*T - S45*T14 - T*T14 + (S35 + T)*pow(m2,2) + pow(S45,2)) + pow(m2,2)*pow(T,2) + pow(m1,2)*(-(S35*S45) + S12*(S35 + S45) - S35*T14 - S45*T14 + (-2*S12 + S35 + S45 - 3*T - 3*T14)*pow(m2,2) + 2*pow(m2,4) - pow(S12,2) + pow(T14,2)))*pow(S12 - S35 - S45 + T + T14,-1))*pow(T14 - pow(MassGlu,2),-1) - 2*pow(S35 - pow(MassGlu,2),-1)*((-(pow(m1,2)*(T14*(S12 - S35 + T14) + (-3*S12 + S35 - 3*T14)*pow(m2,2) + 2*pow(m2,4))) + S35*(-(S12*(S35 + S45 - T14)) - (S45 + T)*T14 + (-S12 + S45 + T)*pow(m2,2) + pow(S12,2)))*pow(S12 - S35 + T14,-1) + 9*(S35*(S12 - S45 - T)*(-T14 + pow(m2,2)) + pow(m1,2)*(-T14 + pow(m2,2))*(-S12 + S35 - T14 + 2*pow(m2,2)) + S12*(S12 + S35 - S45)*pow(MassGlu,2))*pow(T14 - pow(MassGlu,2),-1))) + 16*(S35*(S12 - S45) + pow(m1,2)*pow(m2,2))*pow(S12,-1)*pow(-S35 + pow(MassGlu,2),-1) + 2*pow(S35 - pow(MassGlu,2),-1)*(16*pow(m2,2)*(-(S35*T) + pow(m1,2)*(-S12 + S35 - T14 + pow(m2,2)))*pow(S12 - S35 + T14,-2) - 18*(pow(m1,2)*(T14*(S12 - S35 + T14) + pow(m2,2)*(S35 - T14 - 2*pow(MassGlu,2))) + pow(MassGlu,2)*((-S12 + S45 + T)*T14 + (S12 - S45 + T)*pow(m2,2) + S12*pow(MassGlu,2)))*pow(-T14 + pow(MassGlu,2),-2) + 9*(-(pow(m2,2)*(2*S35*T + (S12 - S45 + T)*pow(MassGlu,2))) + pow(m1,2)*(-(T14*(S12 - S35 + T14)) + 2*pow(m2,4) + pow(m2,2)*(-S12 + S35 - T14 + 2*pow(MassGlu,2))) + pow(MassGlu,2)*(-(S12*(S35 + S45 - T14)) - (S45 + T)*T14 + pow(S12,2)))*pow(S12 - S35 + T14,-1)*pow(-T14 + pow(MassGlu,2),-1)) + pow(T14 - pow(MassGlu,2),-1)*(pow(-S12 + S45 - T,-1)*(2*(S12 - S35 - S45)*(-(S12*S35) - S12*S45 + S12*T - S35*T + S12*T14 - S45*T14 - (S12 - S45 + T)*pow(m2,2) + pow(m1,2)*(-S12 + S35 - T14 + 2*pow(m2,2)) + pow(S12,2))*pow(S12 - S35 + T14,-1) + 9*(-(S12*S35*T14) - 2*S12*S45*T14 + S35*S45*T14 + S12*T*T14 - S35*T*T14 - S45*T*T14 - (S12 - S45 + T)*(S12 - S35 - S45 + T)*pow(m2,2) - 2*pow(m1,4)*pow(m2,2) + S12*S35*pow(MassGlu,2) + S12*S45*pow(MassGlu,2) - S12*T*pow(MassGlu,2) + S35*T*pow(MassGlu,2) + T*T14*pow(MassGlu,2) + pow(m1,2)*(2*S45*T14 + pow(m2,2)*(3*(S12 - S45 + T) - 2*pow(MassGlu,2)) + (S12 - S35 + T14)*pow(MassGlu,2)) + T14*pow(S12,2) - pow(MassGlu,2)*pow(S12,2) + T14*pow(S45,2))*pow(-T14 + pow(MassGlu,2),-1)) - pow(S12 - S35 - S45 + T + T14,-1)*(7*(2*S12*S35*S45 - 2*S12*S35*T - S12*S45*T + S35*S45*T - 2*S12*S35*T14 - 3*S12*S45*T14 + S35*S45*T14 + pow(m1,2)*(S12 - S35 - S45 + T14 + pow(m2,2))*(-S12 + S35 - T14 + 2*pow(m2,2)) - 2*S35*pow(S12,2) - 2*S45*pow(S12,2) + T*pow(S12,2) + 2*T14*pow(S12,2) - pow(m2,2)*(S35*(S45 + T) + (S45 - T)*(S45 - T14) + S12*(-S35 - 2*S45 + T + T14) + pow(S12,2)) + pow(S12,3) + S12*pow(S35,2) + T*pow(S35,2) + S12*pow(S45,2) + T14*pow(S45,2) + S12*pow(T14,2) - S45*pow(T14,2) - T*pow(T14,2))*pow(S12 - S35 + T14,-1) - 9*(S12*S35*T14 + 2*S12*S45*T14 - S35*S45*T14 - S12*T*T14 + S35*T*T14 + S45*T*T14 + 2*pow(m1,4)*pow(m2,2) - S12*S35*pow(MassGlu,2) - S12*S45*pow(MassGlu,2) + S12*T*pow(MassGlu,2) - S35*T*pow(MassGlu,2) + S12*T14*pow(MassGlu,2) + S45*T14*pow(MassGlu,2) + (S12 - S45 + T)*pow(m2,2)*(S12 - S35 - S45 + T + T14 + pow(MassGlu,2)) - pow(m1,2)*(2*(-S12 + S35 + S45 - T14)*T14 + (3*S12 - 2*S35 - 3*S45 + 3*T + 4*T14)*pow(m2,2) + (S12 - S35 + T14)*pow(MassGlu,2)) - T14*pow(S12,2) + pow(MassGlu,2)*pow(S12,2) - T14*pow(S45,2) - S12*pow(T14,2) + S45*pow(T14,2) + T*pow(T14,2))*pow(-T14 + pow(MassGlu,2),-1)))) + pow(T14 - pow(MassGlu,2),-1)*(pow(S12,-1)*(7*(2*pow(m1,4)*pow(m2,2) + T*(S45*(S35 + T14) + (-S45 + T)*pow(m2,2)) + (-S35 - 2*S45 + T)*pow(S12,2) - pow(m1,2)*(-(S12*(S35 + S45 - T14)) + S45*(S35 + T14) + (S12 - S45 + 3*T)*pow(m2,2) + pow(S12,2)) + pow(S12,3) + S12*(S35*(S45 - T) - S45*T + T*(-T14 + pow(m2,2)) + pow(S45,2)))*pow(S12 - S45 + T,-1) + 2*(S35*S45*T + S35*S45*T14 + S35*T*T14 + S45*T*T14 - S35*S45*pow(m2,2) - S35*T*pow(m2,2) - S45*T*pow(m2,2) + 2*pow(m1,4)*pow(m2,2) + (-S35 - 2*S45 + T + T14)*pow(S12,2) + pow(S12,3) + S12*(S35*(S45 - T) - S45*T - S45*T14 - T*T14 + (S35 + T)*pow(m2,2) + pow(S45,2)) + pow(m2,2)*pow(T,2) + pow(m1,2)*(-(S35*S45) + S12*(S35 + S45) - S35*T14 - S45*T14 + (-2*S12 + S35 + S45 - 3*T - 3*T14)*pow(m2,2) + 2*pow(m2,4) - pow(S12,2) + pow(T14,2)))*pow(S12 - S35 - S45 + T + T14,-1))*pow(S35 - pow(MassGlu,2),-1) + 16*(S45*T14 + (S12 - S45 + T)*pow(m2,2) - pow(m1,2)*pow(m2,2))*pow(S12 - S35 - S45 + T + T14,-1)*pow(T14 - pow(MassGlu,2),-1) + pow(S35 - pow(MassGlu,2),-1)*(pow(-S12 + S45 - T,-1)*(-2*(S12 - S35 - S45)*(-(S12*S35) - S12*S45 + S12*T - S35*T + S12*T14 - S45*T14 - (S12 - S45 + T)*pow(m2,2) + pow(m1,2)*(-S12 + S35 - T14 + 2*pow(m2,2)) + pow(S12,2))*pow(S12 - S35 + T14,-1) - 9*(-(S12*S35*T14) - 2*S12*S45*T14 + S35*S45*T14 + S12*T*T14 - S35*T*T14 - S45*T*T14 - (S12 - S45 + T)*(S12 - S35 - S45 + T)*pow(m2,2) - 2*pow(m1,4)*pow(m2,2) + S12*S35*pow(MassGlu,2) + S12*S45*pow(MassGlu,2) - S12*T*pow(MassGlu,2) + S35*T*pow(MassGlu,2) + T*T14*pow(MassGlu,2) + pow(m1,2)*(2*S45*T14 + pow(m2,2)*(3*(S12 - S45 + T) - 2*pow(MassGlu,2)) + (S12 - S35 + T14)*pow(MassGlu,2)) + T14*pow(S12,2) - pow(MassGlu,2)*pow(S12,2) + T14*pow(S45,2))*pow(-T14 + pow(MassGlu,2),-1)) + pow(S12 - S35 - S45 + T + T14,-1)*(7*(2*S12*S35*S45 - 2*S12*S35*T - S12*S45*T + S35*S45*T - 2*S12*S35*T14 - 3*S12*S45*T14 + S35*S45*T14 + pow(m1,2)*(S12 - S35 - S45 + T14 + pow(m2,2))*(-S12 + S35 - T14 + 2*pow(m2,2)) - 2*S35*pow(S12,2) - 2*S45*pow(S12,2) + T*pow(S12,2) + 2*T14*pow(S12,2) - pow(m2,2)*(S35*(S45 + T) + (S45 - T)*(S45 - T14) + S12*(-S35 - 2*S45 + T + T14) + pow(S12,2)) + pow(S12,3) + S12*pow(S35,2) + T*pow(S35,2) + S12*pow(S45,2) + T14*pow(S45,2) + S12*pow(T14,2) - S45*pow(T14,2) - T*pow(T14,2))*pow(S12 - S35 + T14,-1) - 9*(S12*S35*T14 + 2*S12*S45*T14 - S35*S45*T14 - S12*T*T14 + S35*T*T14 + S45*T*T14 + 2*pow(m1,4)*pow(m2,2) - S12*S35*pow(MassGlu,2) - S12*S45*pow(MassGlu,2) + S12*T*pow(MassGlu,2) - S35*T*pow(MassGlu,2) + S12*T14*pow(MassGlu,2) + S45*T14*pow(MassGlu,2) + (S12 - S45 + T)*pow(m2,2)*(S12 - S35 - S45 + T + T14 + pow(MassGlu,2)) - pow(m1,2)*(2*(-S12 + S35 + S45 - T14)*T14 + (3*S12 - 2*S35 - 3*S45 + 3*T + 4*T14)*pow(m2,2) + (S12 - S35 + T14)*pow(MassGlu,2)) - T14*pow(S12,2) + pow(MassGlu,2)*pow(S12,2) - T14*pow(S45,2) - S12*pow(T14,2) + S45*pow(T14,2) + T*pow(T14,2))*pow(-T14 + pow(MassGlu,2),-1))))))/9.;
}