double Sgluons::sgluons_qqbar_OOg_hard(double Alfas, std::array<std::array<double, 4>, 5> const& p) const {
   const double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   const double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   const double k45 = p[3][0]*p[4][0]-p[3][1]*p[4][1]-p[3][2]*p[4][2]-p[3][3]*p[4][3];
   const double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   const double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   const double S12 = 2.*k12;
   const double m1 = mO;
   const double S35 = m1*m1 + 2*k35;
   const double S45 = m1*m1 + 2*k45;
   const double T = m1*m1 - 2*k13;
   const double T14 = m1*m1 - 2*k14;
   return (64*pow<3>(Alfas)*pow<3>(pi)*pow(S12,-2)*(-9*S12*(-2*pow(S12,-1)*(2*(S35 + S45)*(S35*T + (S45 - 2*T)*T14) - 2*(11*S12 - 8*S35 - 4*S45 + 4*T + 8*T14)*pow<4>(m1) + (-3*S45 + T + T14)*pow<2>(S12) + pow<3>(S12) + S12*(-3*S35*T - 2*S45*T - 2*S35*T14 - 3*S45*T14 + 2*T*T14 + 2*pow<2>(S45)) + pow<2>(m1)*(S12*(5*S35 + S45 + 7*T + 15*T14) + 3*pow<2>(S12) - 2*(-(S45*(T - 3*T14)) - 4*T14*(2*T + T14) + S35*(2*S45 + 3*(T + T14)) + pow<2>(S35) + pow<2>(S45)))) + (S35*S45*T - S35*S45*T14 + S35*T*T14 - 3*S45*T*T14 + 2*(S12 - 2*S35 - S45 + T + 2*T14)*pow<4>(m1) - (S35 + 2*S45 - 2*T)*pow<2>(S12) + pow<3>(S12) + T*pow<2>(S35) + T14*pow<2>(S45) - S35*pow<2>(T) + 2*T14*pow<2>(T) + S12*(2*T*T14 + S35*(S45 - 2*T + T14) - S45*(2*T + T14) + pow<2>(S45) + pow<2>(T) - pow<2>(T14)) + S45*pow<2>(T14) - 2*T*pow<2>(T14) - pow<2>(m1)*(-2*S45*T - 2*S12*(S35 + 2*S45 - 2*T - 2*T14) - 4*S35*T14 - 2*S45*T14 + 2*T*T14 + 3*pow<2>(S12) + pow<2>(S35) + pow<2>(S45) + pow<2>(T) + 3*pow<2>(T14)))*pow(S12 - S35 - S45 + T + T14,-1) + (-2*S35*S45*T14 + 3*S35*T*T14 - S45*T*T14 - 2*(S12 - 3*S35 + 3*T14)*pow<4>(m1) + (2*S35 + S45 - 2*T14)*pow<2>(S12) - pow<3>(S12) - S35*pow<2>(T) + 2*T14*pow<2>(T) + S12*(-(S45*T) + 2*S45*T14 - 2*T*T14 + S35*(-2*S45 + T + 2*T14) + pow<2>(T) - pow<2>(T14)) + S45*pow<2>(T14) - 2*T*pow<2>(T14) + pow<2>(m1)*(2*S35*S45 - 3*S35*T + S45*T - S12*(5*S35 + S45 - 8*T14) - 5*S35*T14 - S45*T14 + 2*T*T14 + 3*pow<2>(S12) - pow<2>(T) + 5*pow<2>(T14)))*pow(S12 + T + T14 - 2*pow<2>(m1),-1))*pow(S35 - pow<2>(m1),-1)*pow(S12 - S35 - S45 + 2*pow<2>(m1),-1) + S12*pow(S12 - S35 - S45 + 2*pow<2>(m1),-1)*(-9*(((-S45 + T)*pow<4>(m1) + pow<2>(m1)*(S12*S45 - S35*S45 + 3*S45*T + S45*T14 - 2*pow<2>(S45) - pow<2>(T)) + T*(S35*S45 - 3*S45*T + S12*(-2*S45 + T) - T*T14 + 2*pow<2>(S45) + pow<2>(T)))*pow(S12 - S35 - S45 + T + T14,-1) - ((2*S12 - S35 + T14)*pow<4>(m1) + (S12 - S45 + T)*(-(S35*T) - S45*T - S12*(S35 + S45 - T14) - S45*T14 + T*T14 + pow<2>(S12) - pow<2>(T)) + pow<2>(m1)*(-(S35*S45) + 2*S35*T + S45*T + S45*T14 - 2*T*T14 + S12*(2*S35 + 3*S45 - 2*(T + T14)) - 3*pow<2>(S12) - 2*pow<2>(S45) + pow<2>(T)))*pow(S12 + T + T14 - 2*pow<2>(m1),-1))*pow(S45 - pow<2>(m1),-1) + 2*((68*S12 - 72*(S35 + S45 - T - T14))*pow<6>(m1) + (-40*S35 - 31*S45 + 2*T + 38*T14)*pow<3>(S12) + 11*pow<4>(S12) - 2*pow<4>(m1)*(2*S12*(5*S35 - 22*S45 + 26*T + 8*T14) + 26*pow<2>(S12) - 9*(9*S45*T + S35*(2*S45 - T - T14) + S45*T14 - 2*T*T14 + 2*pow<2>(S35) - 4*pow<2>(S45) - 5*pow<2>(T) - pow<2>(T14))) + pow<2>(S12)*(7*S45*T + S35*(78*S45 - 49*T - 83*T14) - 67*S45*T14 + 58*T*T14 + 29*pow<2>(S35) + 29*pow<2>(S45) - 63*pow<2>(T) + 45*pow<2>(T14)) + pow<2>(m1)*((94*S35 + 22*S45 + 30*T - 78*T14)*pow<2>(S12) - 13*pow<3>(S12) - S12*(147*S45*T + S35*(144*S45 - 83*T - 159*T14) - 65*S45*T14 + 54*T*T14 + 65*pow<2>(S35) + 11*pow<2>(S45) - 131*pow<2>(T) + 85*pow<2>(T14)) + 9*((6*S45 - 3*(T + T14))*pow<2>(S35) + (9*T + T14)*pow<2>(S45) + 2*pow<3>(S45) + S35*(14*T*T14 - 2*S45*(3*T + 7*T14) + 8*pow<2>(S45) + pow<2>(T) + 5*pow<2>(T14)) + S45*(-6*T*T14 - 21*pow<2>(T) + 7*pow<2>(T14)) + 2*(T14*pow<2>(T) + 5*pow<3>(T) - 5*T*pow<2>(T14) - pow<3>(T14)))) + S12*((-29*S45 + 29*T + 27*T14)*pow<2>(S35) + 38*T14*pow<2>(S45) - 9*pow<3>(S45) + S35*(83*S45*T + 119*S45*T14 - 130*T*T14 - 38*pow<2>(S45) + 9*pow<2>(T) - 45*pow<2>(T14)) + 2*S45*(-20*T*T14 + 54*pow<2>(T) - 27*pow<2>(T14)) - 18*(T14*pow<2>(T) + 5*pow<3>(T) - 5*T*pow<2>(T14) - pow<3>(T14))) - 9*((-2*T*T14 + 3*S45*(T + T14))*pow<2>(S35) + (T + T14)*pow<3>(S45) + pow<2>(S45)*(3*pow<2>(T) - pow<2>(T14)) + S35*(4*(T + T14)*pow<2>(S45) - 2*T*(-2*T*T14 + pow<2>(T) - 3*pow<2>(T14)) - S45*(10*T*T14 + pow<2>(T) + 5*pow<2>(T14))) + S45*(-6*T14*pow<2>(T) - 8*pow<3>(T) + 4*T*pow<2>(T14) + 2*pow<3>(T14)) + 4*T*(T - T14)*pow<2>(T + T14)))*pow(S12 - S35 - S45 + T + T14,-1)*pow(S12 + T + T14 - 2*pow<2>(m1),-1)*pow(S12 - S35 - S45 + 2*pow<2>(m1),-1)) + 9*S12*pow(S45 - pow<2>(m1),-1)*(16*pow(S12,-1)*(S45*(T*(S12 - S45 + T) + (S45 - T)*pow<2>(m1))*pow(S45 - pow<2>(m1),-1) + ((-S12 + S45 + pow<2>(m1))*(-(S35*T) - S45*T14 + 2*T*T14 + S12*(-S35 - S45 + T + T14) + (-2*S12 + S35 + S45 - T - T14)*pow<2>(m1) + pow<2>(S12))*pow(-S35 + pow<2>(m1),-1))/8.) - (((-S45 + T)*pow<4>(m1) + pow<2>(m1)*(S12*S45 - S35*S45 + 3*S45*T + S45*T14 - 2*pow<2>(S45) - pow<2>(T)) + T*(S35*S45 - 3*S45*T + S12*(-2*S45 + T) - T*T14 + 2*pow<2>(S45) + pow<2>(T)))*pow(S12 - S35 - S45 + T + T14,-1) - ((2*S12 - S35 + T14)*pow<4>(m1) + (S12 - S45 + T)*(-(S35*T) - S45*T - S12*(S35 + S45 - T14) - S45*T14 + T*T14 + pow<2>(S12) - pow<2>(T)) + pow<2>(m1)*(-(S35*S45) + 2*S35*T + S45*T + S45*T14 - 2*T*T14 + S12*(2*S35 + 3*S45 - 2*(T + T14)) - 3*pow<2>(S12) - 2*pow<2>(S45) + pow<2>(T)))*pow(S12 + T + T14 - 2*pow<2>(m1),-1))*pow(S12 - S35 - S45 + 2*pow<2>(m1),-1)) + 36*((-((S35 + S45)*(S35*T - 2*T*(T + T14) + S45*(2*T + T14))) - 2*(2*S12 - 3*S35 - S45 + T + 3*T14)*pow<4>(m1) + (3*S35 + 3*S45 - 2*(T + T14))*pow<2>(S12) - 2*pow<3>(S12) + S12*(-4*T*T14 + 3*S45*(T + T14) + S35*(-3*S45 + 5*T + T14) - pow<2>(S45)) + pow<2>(m1)*(-3*S45*T - 3*S45*T14 + 4*T*T14 + S12*(-7*S35 - 5*S45 + 4*T + 8*T14) + S35*(4*S45 - 5*(T + T14)) + 6*pow<2>(S12) + pow<2>(S35) + 3*pow<2>(S45) + 4*pow<2>(T14)))*pow(-S12 + S35 + S45 - 2*pow<2>(m1),-2) + 4*pow(S35 - pow<2>(m1),-1)*(pow<2>(m1)*(T14*(S12 - S35 + T14) + (S35 - T14)*pow<2>(m1))*pow(S35 - pow<2>(m1),-1) + ((-S12 + S45 + pow<2>(m1))*(-(S35*T) - S45*T14 + 2*T*T14 + S12*(-S35 - S45 + T + T14) + (-2*S12 + S35 + S45 - T - T14)*pow<2>(m1) + pow<2>(S12))*pow(-S45 + pow<2>(m1),-1))/8.) + (5*S35*S45*T - 2*S45*T*T14 + (2*S12 - S35 + 5*S45 - 5*T + T14)*pow<4>(m1) + pow<2>(m1)*(S12*(2*S35 + 7*S45 + T - 2*T14) - (S45 - T)*(5*S35 + 7*S45 + 6*T - 2*T14) - 3*pow<2>(S12)) + (-S35 - 2*S45 + 4*T + T14)*pow<2>(S12) + pow<3>(S12) + 6*T*pow<2>(S45) + T14*pow<2>(S45) + S12*(S35*(S45 - 5*T) - 2*S45*(5*T + T14) + T*(3*T + 2*T14) + pow<2>(S45)) - 4*S35*pow<2>(T) - 6*S45*pow<2>(T))*pow(-S45 + pow<2>(m1),-1)*pow(S12 - S35 - S45 + 2*pow<2>(m1),-1))))/9.;
}
