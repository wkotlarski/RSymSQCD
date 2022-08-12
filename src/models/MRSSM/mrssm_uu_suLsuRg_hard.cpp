double MRSSM::matrixHard_uu_suLsuRg(
   double Alfas, std::array<std::array<double, 4>, 5> const &p) const {
   const double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   const double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   const double k45 = p[3][0]*p[4][0]-p[3][1]*p[4][1]-p[3][2]*p[4][2]-p[3][3]*p[4][3];
   const double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   const double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   const double S12 = 2*k12;
   const double m1 = MassSq;
   const double m2 = MassSq;
   const double S35 = m1*m1 + 2*k35;
   const double S45 = m2*m2 + 2*k45;
   const double T = m1*m1 - 2*k13;
   const double T14 = m2*m2 - 2*k14;
   return (64*pow<3>(Alfas)*pi_cubed*(pow<-1>(T14 - pow<2>(MassGlu))*(7*(-(S12*S35*S45) - S12*S45*T + S12*S35*T14 - S35*S45*T14 + S12*T*T14 + S35*T*T14 - S45*T*T14 + S12*pow<4>(m1) - (S35 + T)*(S12 - S45 + T)*pow<2>(m2) + S35*pow<2>(S12) + T*pow<2>(S12) + S12*pow<2>(T) - S35*pow<2>(T) + 2*T14*pow<2>(T) + pow<2>(m1)*(S35*T - 2*S35*T14 + S45*T14 - 2*T*T14 + S12*(-S35 + S45 - 2*T + T14) + (S12 + 2*S35 - S45 + T - 2*T14)*pow<2>(m2) - pow<2>(S12) + 2*pow<2>(T14)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1))*pow<-1>(-S35 + pow<2>(m1))*pow<-1>(-S12 + S45 - T + pow<2>(m1) - pow<2>(MassGlu)) + 4*pow<-1>(T14 - pow<2>(MassGlu))*(4*(T*T14 - pow<2>(m1)*pow<2>(m2))*pow<-1>(S12 - S35 - S45 + T + T14) - (9*(S12*S35*T14 - S35*S45*T14 + S35*T*T14 + S12*S45*pow<2>(MassGlu) - S12*T*pow<2>(MassGlu) + S35*T*pow<2>(MassGlu) - S12*T14*pow<2>(MassGlu) + S45*T14*pow<2>(MassGlu) - 2*T*T14*pow<2>(MassGlu) + (S12 - S45 + T)*pow<2>(m2)*(-S35 + pow<2>(MassGlu)) + pow<2>(m1)*(T14*(S12 - 2*S35 + S45 - T + 2*T14) + (S12 + 2*S35 - S45 + T - 2*T14)*pow<2>(m2) + (2*S12 - S35 + T14)*pow<2>(MassGlu)) - pow<2>(MassGlu)*pow<2>(S12))*pow<-1>(-S35 + pow<2>(m1))*pow<-1>(-S12 + S45 - T + pow<2>(m1) - pow<2>(MassGlu)))/4.)) + pow<-1>(S12 - S35 - S45 + T + T14)*pow<-1>(T - pow<2>(MassGlu))*(-((-2*(-(S12*S35*S45) + S12*S35*T + 3*S12*S45*T + 2*S12*S45*T14 - 4*S12*T*T14 + S35*T*T14 + 3*S45*T*T14 + S35*pow<2>(S12) + 2*S45*pow<2>(S12) - 2*T*pow<2>(S12) - T14*pow<2>(S12) - pow<3>(S12) - S12*pow<2>(S45) - T*pow<2>(S45) - T14*pow<2>(S45) + pow<2>(m2)*(-3*S45*T + S12*(-2*S45 + 3*T) + T*(-S35 + 2*T + T14) + pow<2>(S12) + pow<2>(S45)) - 2*S12*pow<2>(T) + 2*S45*pow<2>(T) - 2*T14*pow<2>(T) - pow<3>(T) - T*pow<2>(T14) + pow<2>(m1)*(-2*S45*T - S12*(S35 + S45 - T - 2*T14) - S35*T14 - S45*T14 + T*T14 + (S35 + S45 - T - T14)*pow<2>(m2) + pow<2>(S12) + pow<2>(S45) + pow<2>(T) + pow<2>(T14)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)) + 9*(T*(2*S35*S45 - 2*S35*T + 2*S45*T - 2*S35*T14 - S45*T14 + T*T14 + S12*(-3*S35 + T14) + (S12 - S35 - S45 + T + T14)*pow<2>(MassGlu) + pow<2>(S12) + 2*pow<2>(S35) - pow<2>(S45) - pow<2>(T)) + pow<2>(m1)*(4*S12*S35 + 2*S12*S45 - 2*S35*S45 - 2*S12*T + 2*S35*T - 2*S45*T - 4*S12*T14 + 4*S35*T14 + 2*S45*T14 - 2*T*T14 + (2*S12 - 2*S35 - S45 + T + 2*T14)*pow<2>(m2) - (S12 - S35 - S45 + T + T14)*pow<2>(MassGlu) - 2*pow<2>(S12) - 2*pow<2>(S35) + pow<2>(S45) + pow<2>(T) - 2*pow<2>(T14)))*pow<-1>(T - pow<2>(MassGlu)))*pow<-1>(-S12 + S35 - T14 + pow<2>(m2) - pow<2>(MassGlu))) + 16*(T*(S12 - S35 - S45 + T) + pow<2>(m1)*(-S12 + S35 + S45 - T - T14 + pow<2>(m2)))*pow<-1>(-T + pow<2>(MassGlu))) + pow<-1>(S45 - pow<2>(m2))*(pow<-1>(-S12 + S45 - T + pow<2>(m1) - pow<2>(MassGlu))*((7*(T*(-(S35*T) + S45*T - S45*T14 + 2*T*T14 + S12*(-S35 + T + T14) + pow<2>(S12) - pow<2>(S45)) + pow<2>(m1)*(-(S45*T) - T*(S12 - S35 + T14) + (S45 - T)*pow<2>(m2) + pow<2>(S45)))*pow<-1>(S12 - S35 - S45 + T + T14) - 2*(-S12 + S45 + pow<2>(m1))*(S12*S35 + S12*S45 - S12*T + S35*T - S12*T14 + S45*T14 - 2*T*T14 + (S12 - S35 + T14)*pow<2>(m1) + (S12 - S45 + T)*pow<2>(m2) - pow<2>(S12))*pow<-1>(-S35 + pow<2>(m1)))*pow<-1>(T14 - pow<2>(MassGlu)) + (-(((S12 - S35 + 2*S45 - 2*T + T14)*pow<4>(m1) - (S12 - S45 + T)*(S35*T + S45*T + S12*(S35 + S45 + T - T14) + S45*T14 + (S12 - S45)*pow<2>(m2) - pow<2>(S12) + 2*pow<2>(T)) + pow<2>(m1)*(-(S35*S45) + 2*S35*T - 3*S45*T + 2*S12*(S35 + T - T14) - T*T14 + S12*pow<2>(m2) - 2*pow<2>(S12) - pow<2>(S45) + 4*pow<2>(T)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1))) + 9*((S12 - S35 + T14)*pow<4>(m1) + (S12 - S45 + T)*(-(S35*T) + S45*T - S45*T14 + S12*(-S35 - S45 + T + T14) + (-S12 + S45)*pow<2>(m2) + 2*T*pow<2>(MassGlu) + pow<2>(S12)) + pow<2>(m1)*(-(S35*S45) + 2*S35*T - S45*T + 2*S12*(S35 + S45 - T - T14) - T*T14 + S12*pow<2>(m2) + 2*S45*pow<2>(MassGlu) - 2*T*pow<2>(MassGlu) - 2*pow<2>(S12) + pow<2>(S45)))*pow<-1>(T14 - pow<2>(MassGlu)))*pow<-1>(-S12 + S45 - T + pow<2>(m1) - pow<2>(MassGlu))) + pow<-1>(T - pow<2>(MassGlu))*((7*((S12 - S35 + T14)*pow<4>(m1) - (S12 - S45 + T)*(S35*T + S45*T + S12*(S35 + S45 - T - T14) + S45*T14 - 2*T*T14 + (S12 - S45 + 2*T)*pow<2>(m2) - pow<2>(S12)) + pow<2>(m1)*(-(S35*S45) + 2*S35*T + S45*T + 2*S12*(S35 + S45 - T - T14) + 2*S45*T14 - 3*T*T14 + (S12 - 2*S45 + 2*T)*pow<2>(m2) - 2*pow<2>(S12) - pow<2>(S45)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)) - 2*(-S12 + S45 + pow<2>(m1))*(S12*S35 + S12*S45 - S12*T + S35*T - S12*T14 + S45*T14 - 2*T*T14 + (S12 - S35 + T14)*pow<2>(m1) + (S12 - S45 + T)*pow<2>(m2) - pow<2>(S12))*pow<-1>(-S35 + pow<2>(m1)) + 9*(T*(2*S35*S45 - S35*T + S45*T - S45*T14 + S12*(-S35 + T + T14) + 2*(S12 - S45 + T)*pow<2>(MassGlu) + pow<2>(S12) - pow<2>(S45)) + pow<2>(m1)*(2*S12*S45 - 2*S35*S45 - S12*T + S35*T - S45*T + 2*S45*T14 - T*T14 + (-S45 + T)*pow<2>(m2) + 2*(S45 - T)*pow<2>(MassGlu) + pow<2>(S45)))*pow<-1>(T - pow<2>(MassGlu)))*pow<-1>(-S12 + S35 - T14 + pow<2>(m2) - pow<2>(MassGlu)) + (T*(-2*S35*S45 + S35*T + 5*S45*T + S12*(S35 + 4*S45 - 3*T - T14) + S45*T14 - pow<2>(S12) - 3*pow<2>(S45) - 2*pow<2>(T)) + pow<2>(m1)*(-2*S12*S45 + 2*S35*S45 + S12*T - S35*T - 5*S45*T - 2*S45*T14 + T*T14 + (S45 - T)*pow<2>(m2) + 3*pow<2>(S45) + 2*pow<2>(T)))*pow<-1>(S12 - S35 - S45 + T + T14)*pow<-1>(-T + pow<2>(MassGlu)))) + (16*((S35 + T)*(S12 - S45 + T) + (-S35 + S45 - T + T14)*pow<2>(m1))*(S12 + T + T14 - pow<2>(m1) - pow<2>(m2))*pow<-2>(S12 + T + T14 - 2*pow<2>(m1)) + 9*((-S35 - S45 + T + T14)*pow<4>(m1) + (S12 - S45 + T)*(-(S12*S35) + S12*T - S35*T - 2*S35*T14 - T*T14 + pow<2>(m2)*(2*S35 + T - pow<2>(MassGlu)) + (S12 + T + T14)*pow<2>(MassGlu) + pow<2>(T)) - pow<2>(m1)*(-2*S12*S35 - S12*S45 + S35*S45 + 2*S12*T - 2*S35*T - 2*S45*T + S12*T14 - 2*S35*T14 + S45*T14 + (2*S35 - S45 + T - 2*T14)*pow<2>(m2) + (S12 - S45 + T)*pow<2>(MassGlu) + 2*pow<2>(T) + 2*pow<2>(T14)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1))*pow<-1>(T14 - pow<2>(MassGlu)) + 18*pow<-1>(T14 - pow<2>(MassGlu))*(-(((S35 + S45 - T - T14)*pow<4>(m1) - (S12 - S45 + T)*(-(S12*S35) + S12*T - S35*T - 2*S35*T14 - T*T14 + pow<2>(m2)*(2*S35 + T - pow<2>(MassGlu)) + (S12 + T + T14)*pow<2>(MassGlu) + pow<2>(T)) + pow<2>(m1)*(-2*S12*S35 - S12*S45 + S35*S45 + 2*S12*T - 2*S35*T - 2*S45*T + S12*T14 - 2*S35*T14 + S45*T14 + (2*S35 - S45 + T - 2*T14)*pow<2>(m2) + (S12 - S45 + T)*pow<2>(MassGlu) + 2*pow<2>(T) + 2*pow<2>(T14)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)))/2. + 2*(pow<2>(m1)*((S35 - T14)*(-T14 + pow<2>(m2)) + (S12 + S45 - T)*pow<2>(MassGlu)) + (S12 - S45 + T)*(S35*T14 - (S12 - T + T14)*pow<2>(MassGlu) + pow<2>(m2)*(-S35 + pow<2>(MassGlu))))*pow<-1>(T14 - pow<2>(MassGlu))))*pow<-2>(S12 - S45 + T - pow<2>(m1) + pow<2>(MassGlu)) + (-16*(-S12 - T - T14 + pow<2>(m1) + pow<2>(m2))*((S12 - S35 + T14)*pow<2>(m1) + (S12 - S45 + T)*(-S12 + S35 - T14 + pow<2>(m2)))*pow<-2>(S12 + T + T14 - 2*pow<2>(m1)) + 18*pow<-1>(T - pow<2>(MassGlu))*((S12*S35*T + S12*S45*T - S35*S45*T - S35*T*T14 + S45*T*T14 - (S12 - S35 + T14)*pow<4>(m1) - S12*S45*pow<2>(MassGlu) + 2*S12*T*pow<2>(MassGlu) - S45*T*pow<2>(MassGlu) + S12*T14*pow<2>(MassGlu) - S45*T14*pow<2>(MassGlu) + T*T14*pow<2>(MassGlu) + pow<2>(m2)*(T*(S12 + S35 - S45 + T - T14) - (S12 - S45 + T)*pow<2>(MassGlu)) - pow<2>(m1)*(-((S12 - S45 + 2*T - T14)*(S12 - S35 + T14)) + (S12 + S35 - S45 + T - T14)*pow<2>(m2) + (S12 - S45 + T)*pow<2>(MassGlu)) - T*pow<2>(S12) + pow<2>(MassGlu)*pow<2>(S12) - S12*pow<2>(T) + S35*pow<2>(T) - T14*pow<2>(T) + pow<2>(MassGlu)*pow<2>(T) + T*pow<2>(T14))*pow<-1>(-2*(S12 + T + T14) + 4*pow<2>(m1)) + 2*(T*(S35*(S12 - S35 + T14) + (S12 + S35 - S45 + T - T14)*pow<2>(MassGlu)) + pow<2>(m1)*(-((S12 - S35 + T14)*pow<2>(m2)) + (S12 - S35 + S45 - T + T14)*pow<2>(MassGlu) + pow<2>(S12 - S35 + T14)))*pow<-1>(T - pow<2>(MassGlu))) - 9*(-(S12*S35*T) - S12*S45*T + S35*S45*T + S35*T*T14 - S45*T*T14 + (S12 - S35 + T14)*pow<4>(m1) + S12*S45*pow<2>(MassGlu) - 2*S12*T*pow<2>(MassGlu) + S45*T*pow<2>(MassGlu) - S12*T14*pow<2>(MassGlu) + S45*T14*pow<2>(MassGlu) - T*T14*pow<2>(MassGlu) + pow<2>(m2)*(-(T*(S12 + S35 - S45 + T - T14)) + (S12 - S45 + T)*pow<2>(MassGlu)) + pow<2>(m1)*(-((S12 - S45 + 2*T - T14)*(S12 - S35 + T14)) + (S12 + S35 - S45 + T - T14)*pow<2>(m2) + (S12 - S45 + T)*pow<2>(MassGlu)) + T*pow<2>(S12) - pow<2>(MassGlu)*pow<2>(S12) + S12*pow<2>(T) - S35*pow<2>(T) + T14*pow<2>(T) - pow<2>(MassGlu)*pow<2>(T) - T*pow<2>(T14))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1))*pow<-1>(-T + pow<2>(MassGlu)))*pow<-2>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu)) + pow<-1>(S12 - S35 - S45 + T + T14)*(-2*pow<-1>(T14 - pow<2>(MassGlu))*pow<-1>(-S12 + S45 - T + pow<2>(m1) - pow<2>(MassGlu))*(-2*(-2*S12*S35*T - 2*S12*S45*T + S35*S45*T - S12*S35*T14 - S12*S45*T14 + S35*S45*T14 + 2*S12*T*T14 - S35*T*T14 - S45*T*T14 + (S12 - S35 - S45 + T + T14)*pow<4>(m1) + (S12 - S45 + T)*(S35 - T14)*pow<2>(m2) + 2*T*pow<2>(S12) + T14*pow<2>(S12) + pow<2>(m1)*(-(S35*S45) + 3*S35*T + 2*S45*T + S12*(S35 + 2*S45 - 4*T - T14) + S45*T14 - 3*T*T14 + S12*pow<2>(m2) - pow<2>(S12) - 2*pow<2>(T)) + 3*S12*pow<2>(T) - 2*S35*pow<2>(T) - S45*pow<2>(T) + 2*T14*pow<2>(T) + pow<3>(T) + S12*pow<2>(T14) - S45*pow<2>(T14) + T*pow<2>(T14))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)) - 9*(-(S12*S35*T14) - S12*S45*T14 + S35*S45*T14 - S35*T*T14 + S45*T*T14 + (S12 - S45 + T)*(S35 - T14)*pow<2>(m2) - S12*T*pow<2>(MassGlu) + S35*T*pow<2>(MassGlu) + S45*T*pow<2>(MassGlu) - T*T14*pow<2>(MassGlu) + pow<2>(m1)*((S12 - S45 + T)*pow<2>(m2) + (S12 - S35 - S45 + T + T14)*pow<2>(MassGlu)) + T14*pow<2>(S12) - T14*pow<2>(T) - pow<2>(MassGlu)*pow<2>(T) + S12*pow<2>(T14) - S45*pow<2>(T14) + T*pow<2>(T14))*pow<-1>(-T14 + pow<2>(MassGlu))) + pow<-1>(-T + pow<2>(MassGlu))*(2*(-(S12*S35*S45) + S12*S35*T + 3*S12*S45*T + 2*S12*S45*T14 - 4*S12*T*T14 + S35*T*T14 + 3*S45*T*T14 + S35*pow<2>(S12) + 2*S45*pow<2>(S12) - 2*T*pow<2>(S12) - T14*pow<2>(S12) - pow<3>(S12) - S12*pow<2>(S45) - T*pow<2>(S45) - T14*pow<2>(S45) + pow<2>(m2)*(-3*S45*T + S12*(-2*S45 + 3*T) + T*(-S35 + 2*T + T14) + pow<2>(S12) + pow<2>(S45)) - 2*S12*pow<2>(T) + 2*S45*pow<2>(T) - 2*T14*pow<2>(T) - pow<3>(T) - T*pow<2>(T14) + pow<2>(m1)*(-2*S45*T - S12*(S35 + S45 - T - 2*T14) - S35*T14 - S45*T14 + T*T14 + (S35 + S45 - T - T14)*pow<2>(m2) + pow<2>(S12) + pow<2>(S45) + pow<2>(T) + pow<2>(T14)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)) + 9*(T*(2*S35*S45 - 2*S35*T + 2*S45*T - 2*S35*T14 - S45*T14 + T*T14 + S12*(-3*S35 + T14) + (S12 - S35 - S45 + T + T14)*pow<2>(MassGlu) + pow<2>(S12) + 2*pow<2>(S35) - pow<2>(S45) - pow<2>(T)) + pow<2>(m1)*(4*S12*S35 + 2*S12*S45 - 2*S35*S45 - 2*S12*T + 2*S35*T - 2*S45*T - 4*S12*T14 + 4*S35*T14 + 2*S45*T14 - 2*T*T14 + (2*S12 - 2*S35 - S45 + T + 2*T14)*pow<2>(m2) - (S12 - S35 - S45 + T + T14)*pow<2>(MassGlu) - 2*pow<2>(S12) - 2*pow<2>(S35) + pow<2>(S45) + pow<2>(T) - 2*pow<2>(T14)))*pow<-1>(-T + pow<2>(MassGlu)))*pow<-1>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu))) + pow<-1>(S45 - pow<2>(m2))*(pow<-1>(S12 - S35 - S45 + T + T14)*(7*(T*(-(S35*T) + S45*T - S45*T14 + 2*T*T14 + S12*(-S35 + T + T14) + pow<2>(S12) - pow<2>(S45)) + pow<2>(m1)*(-(S45*T) - T*(S12 - S35 + T14) + (S45 - T)*pow<2>(m2) + pow<2>(S45)))*pow<-1>(T14 - pow<2>(MassGlu))*pow<-1>(-S12 + S45 - T + pow<2>(m1) - pow<2>(MassGlu)) + (-(pow<2>(m1)*(-2*S12*S45 + 2*S35*S45 + S12*T - S35*T - 5*S45*T - 2*S45*T14 + T*T14 + (S45 - T)*pow<2>(m2) + 3*pow<2>(S45) + 2*pow<2>(T))) + T*(2*S35*S45 - S35*T - 5*S45*T - S45*T14 + S12*(-S35 - 4*S45 + 3*T + T14) + pow<2>(S12) + 3*pow<2>(S45) + 2*pow<2>(T)))*pow<-2>(-T + pow<2>(MassGlu))) + (-(((S12 - S35 + 2*S45 - 2*T + T14)*pow<4>(m1) - (S12 - S45 + T)*(S35*T + S45*T + S12*(S35 + S45 + T - T14) + S45*T14 + (S12 - S45)*pow<2>(m2) - pow<2>(S12) + 2*pow<2>(T)) + pow<2>(m1)*(-(S35*S45) + 2*S35*T - 3*S45*T + 2*S12*(S35 + T - T14) - T*T14 + S12*pow<2>(m2) - 2*pow<2>(S12) - pow<2>(S45) + 4*pow<2>(T)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1))) + 9*((S12 - S35 + T14)*pow<4>(m1) + (S12 - S45 + T)*(-(S35*T) + S45*T - S45*T14 + S12*(-S35 - S45 + T + T14) + (-S12 + S45)*pow<2>(m2) + 2*T*pow<2>(MassGlu) + pow<2>(S12)) + pow<2>(m1)*(-(S35*S45) + 2*S35*T - S45*T + 2*S12*(S35 + S45 - T - T14) - T*T14 + S12*pow<2>(m2) + 2*S45*pow<2>(MassGlu) - 2*T*pow<2>(MassGlu) - 2*pow<2>(S12) + pow<2>(S45)))*pow<-1>(T14 - pow<2>(MassGlu)))*pow<-2>(S12 - S45 + T - pow<2>(m1) + pow<2>(MassGlu)) + 32*S45*(T*(S12 - S45 + T) + (S45 - T)*pow<2>(m1))*pow<-1>(S45 - pow<2>(m2))*(pow<-2>(-T + pow<2>(MassGlu)) + pow<-2>(S12 - S45 + T - pow<2>(m1) + pow<2>(MassGlu))) - pow<-1>(-T + pow<2>(MassGlu))*(-7*((S12 - S35 + T14)*pow<4>(m1) - (S12 - S45 + T)*(S35*T + S45*T + S12*(S35 + S45 - T - T14) + S45*T14 - 2*T*T14 + (S12 - S45 + 2*T)*pow<2>(m2) - pow<2>(S12)) + pow<2>(m1)*(-(S35*S45) + 2*S35*T + S45*T + 2*S12*(S35 + S45 - T - T14) + 2*S45*T14 - 3*T*T14 + (S12 - 2*S45 + 2*T)*pow<2>(m2) - 2*pow<2>(S12) - pow<2>(S45)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)) + 9*(T*(2*S35*S45 - S35*T + S45*T - S45*T14 + S12*(-S35 + T + T14) + 2*(S12 - S45 + T)*pow<2>(MassGlu) + pow<2>(S12) - pow<2>(S45)) + pow<2>(m1)*(2*S12*S45 - 2*S35*S45 - S12*T + S35*T - S45*T + 2*S45*T14 - T*T14 + (-S45 + T)*pow<2>(m2) + 2*(S45 - T)*pow<2>(MassGlu) + pow<2>(S45)))*pow<-1>(-T + pow<2>(MassGlu)))*pow<-1>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu)) - 2*(-S12 + S45 + pow<2>(m1))*(S12*T - S35*T + S12*T14 - S45*T14 + 2*T*T14 - 2*S12*pow<2>(MassGlu) + S35*pow<2>(MassGlu) + S45*pow<2>(MassGlu) + pow<2>(m2)*(-T + pow<2>(MassGlu)) + pow<2>(m1)*(-T14 + pow<2>(MassGlu)) - 2*pow<4>(MassGlu))*(S12*S35 + S12*S45 - S12*T + S35*T - S12*T14 + S45*T14 - 2*T*T14 + (S12 - S35 + T14)*pow<2>(m1) + (S12 - S45 + T)*pow<2>(m2) - pow<2>(S12))*pow<-1>(S35 - pow<2>(m1))*pow<-1>(-T + pow<2>(MassGlu))*pow<-1>(-T14 + pow<2>(MassGlu))*pow<-1>(S12 - S45 + T - pow<2>(m1) + pow<2>(MassGlu))*pow<-1>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu))) + pow<-1>(S35 - pow<2>(m1))*(-((((S12 - S35 + T14)*pow<4>(m1) + pow<2>(m1)*(-((S12 - S35 + T14)*(2*S12 + T + 3*T14)) + (2*S12 - 3*S35 + 3*T14)*pow<2>(m2)) + (S12 - S45 + T)*pow<4>(m2) + pow<2>(m2)*(-(S35*S45) + 2*S35*T + S12*(S35 + 2*S45 - 2*T - 2*T14) + 2*S45*T14 - 3*T*T14 - 2*pow<2>(S12)) + (S12 - S35 + T14)*(-((S45 - 2*T)*T14) + S12*(-S45 + T + T14) + pow<2>(S12)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)) + 9*(S12*S35*T + S35*T*T14 + (S12 - S35 + T14)*pow<4>(m1) + S12*S45*pow<2>(MassGlu) - S12*T*pow<2>(MassGlu) + S35*T*pow<2>(MassGlu) - S12*T14*pow<2>(MassGlu) + S45*T14*pow<2>(MassGlu) - 2*T*T14*pow<2>(MassGlu) + (S12 - S45 + T)*pow<2>(m2)*pow<2>(MassGlu) + pow<2>(m1)*(-((S35 + T - 2*T14)*(S12 - S35 + T14)) + 2*(S35 - T14)*pow<2>(m2) + (2*S12 - S35 + T14)*pow<2>(MassGlu)) - pow<2>(MassGlu)*pow<2>(S12) - T*pow<2>(S35))*pow<-1>(-T + pow<2>(MassGlu)))*pow<-2>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu))) + 4*pow<-1>(S12 - S35 - S45 + T + T14)*(-(((S12 - S45 + T)*(S35 - T14)*pow<2>(m2) + pow<2>(m1)*(T14*(S12 - S35 - S45 + T + T14) - (S12 - 2*S35 - S45 + T + 2*T14)*pow<2>(m2)) + T14*((S45 - 2*T)*(S35 - T14) + S12*(-S35 - S45 + T + T14) + pow<2>(S12)))*pow<-2>(-T14 + pow<2>(MassGlu)))/4. - (7*(S12*S35*S45 - 2*S12*S35*T - 2*S12*S45*T + S35*S45*T - 2*S12*S45*T14 + 3*S12*T*T14 - S35*T*T14 - 3*S45*T*T14 + pow<2>(m1)*(-((S12 - S35 - S45 + T + T14)*(2*S12 - S35 + 2*T14)) + (S12 - 2*S35 - S45 + T + 2*T14)*pow<2>(m2)) - S35*pow<2>(S12) - 2*S45*pow<2>(S12) + 2*T*pow<2>(S12) + T14*pow<2>(S12) + pow<3>(S12) + T*pow<2>(S35) + S12*pow<2>(S45) + T14*pow<2>(S45) + S12*pow<2>(T) - S35*pow<2>(T) + 2*T14*pow<2>(T) - pow<2>(m2)*pow<2>(S12 - S45 + T))*pow<-1>(-T + pow<2>(MassGlu))*pow<-1>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu)))/4.)) + pow<-1>(S35 - pow<2>(m1))*(pow<-1>(T14 - pow<2>(MassGlu))*pow<-1>(-S12 + S45 - T + pow<2>(m1) - pow<2>(MassGlu))*(-7*(-(S12*S35*S45) - S12*S45*T + S12*S35*T14 - S35*S45*T14 + S12*T*T14 + S35*T*T14 - S45*T*T14 + S12*pow<4>(m1) - (S35 + T)*(S12 - S45 + T)*pow<2>(m2) + S35*pow<2>(S12) + T*pow<2>(S12) + S12*pow<2>(T) - S35*pow<2>(T) + 2*T14*pow<2>(T) + pow<2>(m1)*(S35*T - 2*S35*T14 + S45*T14 - 2*T*T14 + S12*(-S35 + S45 - 2*T + T14) + (S12 + 2*S35 - S45 + T - 2*T14)*pow<2>(m2) - pow<2>(S12) + 2*pow<2>(T14)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)) - 9*(S12*S35*T14 - S35*S45*T14 + S35*T*T14 + S12*S45*pow<2>(MassGlu) - S12*T*pow<2>(MassGlu) + S35*T*pow<2>(MassGlu) - S12*T14*pow<2>(MassGlu) + S45*T14*pow<2>(MassGlu) - 2*T*T14*pow<2>(MassGlu) + (S12 - S45 + T)*pow<2>(m2)*(-S35 + pow<2>(MassGlu)) + pow<2>(m1)*(T14*(S12 - 2*S35 + S45 - T + 2*T14) + (S12 + 2*S35 - S45 + T - 2*T14)*pow<2>(m2) + (2*S12 - S35 + T14)*pow<2>(MassGlu)) - pow<2>(MassGlu)*pow<2>(S12))*pow<-1>(-T14 + pow<2>(MassGlu))) - (((S12 - S35 + T14)*pow<4>(m1) + pow<2>(m1)*(-((S12 - S35 + T14)*(2*S12 + T + 3*T14)) + (2*S12 - 3*S35 + 3*T14)*pow<2>(m2)) + (S12 - S45 + T)*pow<4>(m2) + pow<2>(m2)*(-(S35*S45) + 2*S35*T + S12*(S35 + 2*S45 - 2*T - 2*T14) + 2*S45*T14 - 3*T*T14 - 2*pow<2>(S12)) + (S12 - S35 + T14)*(-((S45 - 2*T)*T14) + S12*(-S45 + T + T14) + pow<2>(S12)))*pow<-1>(S12 + T + T14 - 2*pow<2>(m1)) + 9*(S12*S35*T + S35*T*T14 + (S12 - S35 + T14)*pow<4>(m1) + S12*S45*pow<2>(MassGlu) - S12*T*pow<2>(MassGlu) + S35*T*pow<2>(MassGlu) - S12*T14*pow<2>(MassGlu) + S45*T14*pow<2>(MassGlu) - 2*T*T14*pow<2>(MassGlu) + (S12 - S45 + T)*pow<2>(m2)*pow<2>(MassGlu) + pow<2>(m1)*(-((S35 + T - 2*T14)*(S12 - S35 + T14)) + 2*(S35 - T14)*pow<2>(m2) + (2*S12 - S35 + T14)*pow<2>(MassGlu)) - pow<2>(MassGlu)*pow<2>(S12) - T*pow<2>(S35))*pow<-1>(-T + pow<2>(MassGlu)))*pow<-2>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu)) + 32*pow<2>(m1)*(T14*(S12 - S35 + T14) + (S35 - T14)*pow<2>(m2))*pow<-1>(S35 - pow<2>(m1))*(pow<-2>(-T14 + pow<2>(MassGlu)) + pow<-2>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu))) - pow<-1>(S12 - S35 - S45 + T + T14)*(((S12 - S45 + T)*(S35 - T14)*pow<2>(m2) + pow<2>(m1)*(T14*(S12 - S35 - S45 + T + T14) - (S12 - 2*S35 - S45 + T + 2*T14)*pow<2>(m2)) + T14*((S45 - 2*T)*(S35 - T14) + S12*(-S35 - S45 + T + T14) + pow<2>(S12)))*pow<-2>(-T14 + pow<2>(MassGlu)) + 7*(S12*S35*S45 - 2*S12*S35*T - 2*S12*S45*T + S35*S45*T - 2*S12*S45*T14 + 3*S12*T*T14 - S35*T*T14 - 3*S45*T*T14 + pow<2>(m1)*(-((S12 - S35 - S45 + T + T14)*(2*S12 - S35 + 2*T14)) + (S12 - 2*S35 - S45 + T + 2*T14)*pow<2>(m2)) - S35*pow<2>(S12) - 2*S45*pow<2>(S12) + 2*T*pow<2>(S12) + T14*pow<2>(S12) + pow<3>(S12) + T*pow<2>(S35) + S12*pow<2>(S45) + T14*pow<2>(S45) + S12*pow<2>(T) - S35*pow<2>(T) + 2*T14*pow<2>(T) - pow<2>(m2)*pow<2>(S12 - S45 + T))*pow<-1>(-T + pow<2>(MassGlu))*pow<-1>(S12 - S35 + T14 - pow<2>(m2) + pow<2>(MassGlu))))))/27.;
}
