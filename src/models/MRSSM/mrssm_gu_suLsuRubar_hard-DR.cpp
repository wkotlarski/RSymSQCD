double MRSSM::matrixMRSSMHard_gu_suLsuRubar_DR(double Alfas, std::array<std::array<double, 4>, 5> const& p) const {
   double Alfas2 = Sqr(Alfas);
   double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   double k45 = p[3][0]*p[4][0]-p[3][1]*p[4][1]-p[3][2]*p[4][2]-p[3][3]*p[4][3];
   double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   double S12 = 2*k12;
   const double m1 = MassSq;
   const double m2 = MassSq;
   double S35 = m1*m1 + 2*k35;
   double S45 = m2*m2 + 2*k45;
   double T = m1*m1 - 2*k13;
   double T14 = m2*m2 - 2*k14;
   return (-16*pow<3>(Alfas)*pow<3>(pi)*pow(S12,-1)*(-((8*(8*S12*S35*S45*T + 7*S12*S35*S45*T14 - 8*S12*S35*T*T14 - 4*S12*S45*T*T14 + 5*S35*S45*T*T14 - 6*S12*S35*S45*pow<2>(m1) + 8*S12*S35*T*pow<2>(m1) + S12*S45*T*pow<2>(m1) - 4*S35*S45*T*pow<2>(m1) + 9*S12*S35*T14*pow<2>(m1) + 3*S12*S45*T14*pow<2>(m1) - 3*S35*S45*T14*pow<2>(m1) - 5*S12*T*T14*pow<2>(m1) + 7*S35*T*T14*pow<2>(m1) + 2*S45*T*T14*pow<2>(m1) - 3*S12*S35*pow<4>(m1) + S35*S45*pow<4>(m1) - S35*T*pow<4>(m1) + 3*S12*T14*pow<4>(m1) - 4*S35*T14*pow<4>(m1) - S45*T14*pow<4>(m1) + T*T14*pow<4>(m1) - 5*S12*S35*S45*pow<2>(m2) + 4*S12*S35*T*pow<2>(m2) + 5*S12*S45*T*pow<2>(m2) - 3*S35*S45*T*pow<2>(m2) + 2*S12*S35*T14*pow<2>(m2) + 4*S12*S45*T14*pow<2>(m2) - 2*S35*S45*T14*pow<2>(m2) - 3*S12*T*T14*pow<2>(m2) + 2*S35*T*T14*pow<2>(m2) + 3*S45*T*T14*pow<2>(m2) - 3*S12*S35*pow<2>(m1)*pow<2>(m2) - 2*S12*S45*pow<2>(m1)*pow<2>(m2) + S35*S45*pow<2>(m1)*pow<2>(m2) + 2*S12*T*pow<2>(m1)*pow<2>(m2) - S35*T*pow<2>(m1)*pow<2>(m2) - 2*S45*T*pow<2>(m1)*pow<2>(m2) + 3*S12*T14*pow<2>(m1)*pow<2>(m2) - 4*S35*T14*pow<2>(m1)*pow<2>(m2) - S45*T14*pow<2>(m1)*pow<2>(m2) + T*T14*pow<2>(m1)*pow<2>(m2) - 2*S12*S45*pow<4>(m2) + 2*S12*T*pow<4>(m2) - 2*S45*T*pow<4>(m2) + 5*S35*S45*pow<2>(S12) - 7*S35*T*pow<2>(S12) - 3*S45*T*pow<2>(S12) - 5*S35*T14*pow<2>(S12) - 4*S45*T14*pow<2>(S12) + 3*T*T14*pow<2>(S12) + 6*S35*pow<2>(m1)*pow<2>(S12) + 2*S45*pow<2>(m1)*pow<2>(S12) - 2*T*pow<2>(m1)*pow<2>(S12) - 5*T14*pow<2>(m1)*pow<2>(S12) + pow<4>(m1)*pow<2>(S12) + 3*S35*pow<2>(m2)*pow<2>(S12) + 4*S45*pow<2>(m2)*pow<2>(S12) - 4*T*pow<2>(m2)*pow<2>(S12) - 2*T14*pow<2>(m2)*pow<2>(S12) + 2*pow<2>(m1)*pow<2>(m2)*pow<2>(S12) + pow<4>(m2)*pow<2>(S12) - 3*S35*pow<3>(S12) - 2*S45*pow<3>(S12) + 2*T*pow<3>(S12) + 2*T14*pow<3>(S12) - 2*pow<2>(m1)*pow<3>(S12) - 2*pow<2>(m2)*pow<3>(S12) + pow<4>(S12) - 2*S12*S45*pow<2>(S35) + 4*S12*T*pow<2>(S35) - 2*S45*T*pow<2>(S35) + 2*S12*T14*pow<2>(S35) - 2*S45*T14*pow<2>(S35) + 2*T*T14*pow<2>(S35) - 4*S12*pow<2>(m1)*pow<2>(S35) + 2*S45*pow<2>(m1)*pow<2>(S35) - 4*T*pow<2>(m1)*pow<2>(S35) - 2*T14*pow<2>(m1)*pow<2>(S35) + 2*pow<4>(m1)*pow<2>(S35) - 2*S12*pow<2>(m2)*pow<2>(S35) + 2*S45*pow<2>(m2)*pow<2>(S35) - 2*T*pow<2>(m2)*pow<2>(S35) + 2*pow<2>(m1)*pow<2>(m2)*pow<2>(S35) + 2*pow<2>(S12)*pow<2>(S35) - 2*S12*S35*pow<2>(S45) + S12*T*pow<2>(S45) - 2*S35*T*pow<2>(S45) + 2*S12*T14*pow<2>(S45) - 2*S35*T14*pow<2>(S45) + T*T14*pow<2>(S45) - S12*pow<2>(m1)*pow<2>(S45) + 2*S35*pow<2>(m1)*pow<2>(S45) - T14*pow<2>(m1)*pow<2>(S45) - 2*S12*pow<2>(m2)*pow<2>(S45) + 2*S35*pow<2>(m2)*pow<2>(S45) - T*pow<2>(m2)*pow<2>(S45) - 2*T14*pow<2>(m2)*pow<2>(S45) + pow<2>(m1)*pow<2>(m2)*pow<2>(S45) + pow<4>(m2)*pow<2>(S45) + pow<2>(S12)*pow<2>(S45) - 5*S12*S35*pow<2>(T) - S12*S45*pow<2>(T) + 3*S35*S45*pow<2>(T) + S12*T14*pow<2>(T) - 3*S35*T14*pow<2>(T) - S45*T14*pow<2>(T) + 2*S35*pow<2>(m1)*pow<2>(T) - T14*pow<2>(m1)*pow<2>(T) - 3*S12*pow<2>(m2)*pow<2>(T) + S35*pow<2>(m2)*pow<2>(T) + 2*S45*pow<2>(m2)*pow<2>(T) - T14*pow<2>(m2)*pow<2>(T) + pow<2>(m1)*pow<2>(m2)*pow<2>(T) + pow<4>(m2)*pow<2>(T) + pow<2>(S12)*pow<2>(T) + 2*pow<2>(S35)*pow<2>(T) - S35*pow<3>(T) - pow<2>(m2)*pow<3>(T) - 2*S12*S35*pow<2>(T14) - 2*S12*S45*pow<2>(T14) + 2*S35*S45*pow<2>(T14) + S12*T*pow<2>(T14) - 2*S35*T*pow<2>(T14) - S45*T*pow<2>(T14) - 5*S12*pow<2>(m1)*pow<2>(T14) + 4*S35*pow<2>(m1)*pow<2>(T14) + S45*pow<2>(m1)*pow<2>(T14) - 3*T*pow<2>(m1)*pow<2>(T14) + 2*pow<4>(m1)*pow<2>(T14) + 2*pow<2>(m1)*pow<2>(m2)*pow<2>(T14) + pow<2>(S12)*pow<2>(T14) + pow<2>(S45)*pow<2>(T14) - 2*pow<2>(m1)*pow<3>(T14))*pow(-S12 - T - T14 + 2*pow<2>(m1),-2) + 16*(-(S35*(S12 - S45 + T)) + (S35 - T14)*pow<2>(m1))*(T14*(S12 - S35 + T14) + (S35 - T14)*pow<2>(m2))*pow(T14 - pow<2>(m2),-2) - (pow<2>(m1)*((S35 - T14)*(T14*(-S45 + 3*T + 4*T14) + 2*S35*(S45 - 2*(T + T14))) + (6*S35 + 2*S45 - 2*T - 5*T14)*pow<2>(S12) + pow<2>(m2)*((4*S35 + S45 - T - 4*T14)*(S35 - T14) + S12*(-3*S35 - 2*S45 + 2*T + 3*T14) + 2*pow<2>(S12)) - 2*pow<3>(S12) + S12*((3*S45 - 4*T - 7*T14)*T14 + S35*(-5*S45 + 6*T + 11*T14) - 4*pow<2>(S35))) + pow<4>(m1)*(-3*S12*(S35 - T14) + pow<2>(S12) + 2*pow<2>(S35 - T14)) + (S12 - S45 + T)*(2*S35*S45*T14 - 3*S35*T*T14 + (S12 - S45 + T)*pow<4>(m2) + pow<2>(m2)*(S12*(3*S35 + 2*S45 - 2*T - 2*T14) - (4*S35 + 2*S45 - T)*(S35 - T14) - 2*pow<2>(S12)) + (-3*S35 - S45 + T + 2*T14)*pow<2>(S12) + pow<3>(S12) + 2*T*pow<2>(S35) + 4*T14*pow<2>(S35) + S12*(S35*(2*S45 - 3*T - 7*T14) + T14*(-2*S45 + T + T14) + 2*pow<2>(S35)) - 4*S35*pow<2>(T14) - S45*pow<2>(T14)))*pow(-S12 - T - T14 + 2*pow<2>(m1),-1)*pow(T14 - pow<2>(m2),-1))*pow(-S12 + S45 - T + pow<2>(m1) - pow<2>(MassGlu),-2)) + (16*(T*(S12 - S45 + T) + (S45 - T)*pow<2>(m1))*(S45*(S12 - S35 + T14) + (-S45 + T)*pow<2>(m2))*pow(T - pow<2>(m1),-2) - 8*(7*S12*S35*S45*T + 8*S12*S35*S45*T14 - 4*S12*S35*T*T14 - 8*S12*S45*T*T14 + 5*S35*S45*T*T14 - 5*S12*S35*S45*pow<2>(m1) + 4*S12*S35*T*pow<2>(m1) + 2*S12*S45*T*pow<2>(m1) - 2*S35*S45*T*pow<2>(m1) + 5*S12*S35*T14*pow<2>(m1) + 4*S12*S45*T14*pow<2>(m1) - 3*S35*S45*T14*pow<2>(m1) - 3*S12*T*T14*pow<2>(m1) + 3*S35*T*T14*pow<2>(m1) + 2*S45*T*T14*pow<2>(m1) - 2*S12*S35*pow<4>(m1) + 2*S12*T14*pow<4>(m1) - 2*S35*T14*pow<4>(m1) - 6*S12*S35*S45*pow<2>(m2) + 3*S12*S35*T*pow<2>(m2) + 9*S12*S45*T*pow<2>(m2) - 3*S35*S45*T*pow<2>(m2) + S12*S35*T14*pow<2>(m2) + 8*S12*S45*T14*pow<2>(m2) - 4*S35*S45*T14*pow<2>(m2) - 5*S12*T*T14*pow<2>(m2) + 2*S35*T*T14*pow<2>(m2) + 7*S45*T*T14*pow<2>(m2) - 2*S12*S35*pow<2>(m1)*pow<2>(m2) - 3*S12*S45*pow<2>(m1)*pow<2>(m2) + S35*S45*pow<2>(m1)*pow<2>(m2) + 3*S12*T*pow<2>(m1)*pow<2>(m2) - S35*T*pow<2>(m1)*pow<2>(m2) - 4*S45*T*pow<2>(m1)*pow<2>(m2) + 2*S12*T14*pow<2>(m1)*pow<2>(m2) - 2*S35*T14*pow<2>(m1)*pow<2>(m2) - S45*T14*pow<2>(m1)*pow<2>(m2) + T*T14*pow<2>(m1)*pow<2>(m2) - 3*S12*S45*pow<4>(m2) + S35*S45*pow<4>(m2) + 3*S12*T*pow<4>(m2) - S35*T*pow<4>(m2) - 4*S45*T*pow<4>(m2) - S45*T14*pow<4>(m2) + T*T14*pow<4>(m2) + 5*S35*S45*pow<2>(S12) - 4*S35*T*pow<2>(S12) - 5*S45*T*pow<2>(S12) - 3*S35*T14*pow<2>(S12) - 7*S45*T14*pow<2>(S12) + 3*T*T14*pow<2>(S12) + 4*S35*pow<2>(m1)*pow<2>(S12) + 3*S45*pow<2>(m1)*pow<2>(S12) - 2*T*pow<2>(m1)*pow<2>(S12) - 4*T14*pow<2>(m1)*pow<2>(S12) + pow<4>(m1)*pow<2>(S12) + 2*S35*pow<2>(m2)*pow<2>(S12) + 6*S45*pow<2>(m2)*pow<2>(S12) - 5*T*pow<2>(m2)*pow<2>(S12) - 2*T14*pow<2>(m2)*pow<2>(S12) + 2*pow<2>(m1)*pow<2>(m2)*pow<2>(S12) + pow<4>(m2)*pow<2>(S12) - 2*S35*pow<3>(S12) - 3*S45*pow<3>(S12) + 2*T*pow<3>(S12) + 2*T14*pow<3>(S12) - 2*pow<2>(m1)*pow<3>(S12) - 2*pow<2>(m2)*pow<3>(S12) + pow<4>(S12) - 2*S12*S45*pow<2>(S35) + 2*S12*T*pow<2>(S35) - 2*S45*T*pow<2>(S35) + S12*T14*pow<2>(S35) - 2*S45*T14*pow<2>(S35) + T*T14*pow<2>(S35) - 2*S12*pow<2>(m1)*pow<2>(S35) + 2*S45*pow<2>(m1)*pow<2>(S35) - 2*T*pow<2>(m1)*pow<2>(S35) - T14*pow<2>(m1)*pow<2>(S35) + pow<4>(m1)*pow<2>(S35) - S12*pow<2>(m2)*pow<2>(S35) + 2*S45*pow<2>(m2)*pow<2>(S35) - T*pow<2>(m2)*pow<2>(S35) + pow<2>(m1)*pow<2>(m2)*pow<2>(S35) + pow<2>(S12)*pow<2>(S35) - 2*S12*S35*pow<2>(S45) + 2*S12*T*pow<2>(S45) - 2*S35*T*pow<2>(S45) + 4*S12*T14*pow<2>(S45) - 2*S35*T14*pow<2>(S45) + 2*T*T14*pow<2>(S45) - 2*S12*pow<2>(m1)*pow<2>(S45) + 2*S35*pow<2>(m1)*pow<2>(S45) - 2*T14*pow<2>(m1)*pow<2>(S45) - 4*S12*pow<2>(m2)*pow<2>(S45) + 2*S35*pow<2>(m2)*pow<2>(S45) - 2*T*pow<2>(m2)*pow<2>(S45) - 4*T14*pow<2>(m2)*pow<2>(S45) + 2*pow<2>(m1)*pow<2>(m2)*pow<2>(S45) + 2*pow<4>(m2)*pow<2>(S45) + 2*pow<2>(S12)*pow<2>(S45) - 2*S12*S35*pow<2>(T) - 2*S12*S45*pow<2>(T) + 2*S35*S45*pow<2>(T) + S12*T14*pow<2>(T) - S35*T14*pow<2>(T) - 2*S45*T14*pow<2>(T) - 5*S12*pow<2>(m2)*pow<2>(T) + S35*pow<2>(m2)*pow<2>(T) + 4*S45*pow<2>(m2)*pow<2>(T) - 3*T14*pow<2>(m2)*pow<2>(T) + 2*pow<2>(m1)*pow<2>(m2)*pow<2>(T) + 2*pow<4>(m2)*pow<2>(T) + pow<2>(S12)*pow<2>(T) + pow<2>(S35)*pow<2>(T) - 2*pow<2>(m2)*pow<3>(T) - S12*S35*pow<2>(T14) - 5*S12*S45*pow<2>(T14) + 3*S35*S45*pow<2>(T14) + S12*T*pow<2>(T14) - S35*T*pow<2>(T14) - 3*S45*T*pow<2>(T14) - 3*S12*pow<2>(m1)*pow<2>(T14) + 2*S35*pow<2>(m1)*pow<2>(T14) + S45*pow<2>(m1)*pow<2>(T14) - T*pow<2>(m1)*pow<2>(T14) + pow<4>(m1)*pow<2>(T14) + 2*S45*pow<2>(m2)*pow<2>(T14) - T*pow<2>(m2)*pow<2>(T14) + pow<2>(m1)*pow<2>(m2)*pow<2>(T14) + pow<2>(S12)*pow<2>(T14) + 2*pow<2>(S45)*pow<2>(T14) - S45*pow<3>(T14) - pow<2>(m1)*pow<3>(T14))*pow(-S12 - T - T14 + 2*pow<2>(m1),-2) + (9*S12*S35*S45*T + 5*S12*S35*S45*T14 - 3*S12*S35*T*T14 - 10*S12*S45*T*T14 + 5*S35*S45*T*T14 + 5*S35*S45*pow<2>(S12) - 4*S35*T*pow<2>(S12) - 7*S45*T*pow<2>(S12) - 2*S35*T14*pow<2>(S12) - 6*S45*T14*pow<2>(S12) + 3*T*T14*pow<2>(S12) + pow<2>(m1)*(pow<2>(m2)*((S45 - T)*(S35 + 4*S45 - 4*T - T14) + S12*(-2*S35 - 3*S45 + 3*T + 2*T14) + 2*pow<2>(S12)) - (S12 - S35 + T14)*((S45 - T)*(2*S35 + 4*S45 - T14) + S12*(-2*S35 - 3*S45 + 2*(T + T14)) + 2*pow<2>(S12))) - 2*S35*pow<3>(S12) - 3*S45*pow<3>(S12) + 2*T*pow<3>(S12) + 2*T14*pow<3>(S12) + pow<4>(S12) - 2*S12*S45*pow<2>(S35) + 2*S12*T*pow<2>(S35) - 2*S45*T*pow<2>(S35) + pow<2>(S12)*pow<2>(S35) - 2*S12*S35*pow<2>(S45) + 4*S12*T*pow<2>(S45) - 4*S35*T*pow<2>(S45) + 4*S12*T14*pow<2>(S45) - 2*S35*T14*pow<2>(S45) + 4*T*T14*pow<2>(S45) + 2*pow<2>(S12)*pow<2>(S45) + pow<4>(m2)*(-3*S12*(S45 - T) + pow<2>(S12) + 2*pow<2>(S45 - T)) - 2*S12*S35*pow<2>(T) - 4*S12*S45*pow<2>(T) + 4*S35*S45*pow<2>(T) + S12*T14*pow<2>(T) - S35*T14*pow<2>(T) - 4*S45*T14*pow<2>(T) + pow<2>(S12)*pow<2>(T) + pow<2>(S35)*pow<2>(T) - pow<2>(m2)*((S45 - T)*(S35*(-2*S45 + T) + 4*S45*(T + T14) - T*(4*T + 3*T14)) + (-2*S35 - 6*S45 + 5*T + 2*T14)*pow<2>(S12) + 2*pow<3>(S12) + S12*(5*S35*S45 - 3*S35*T - 11*S45*T - 6*S45*T14 + 4*T*T14 + 4*pow<2>(S45) + 7*pow<2>(T))) - 3*S12*S45*pow<2>(T14) + S12*T*pow<2>(T14) - 3*S45*T*pow<2>(T14) + pow<2>(S12)*pow<2>(T14) + 2*pow<2>(S45)*pow<2>(T14) + pow<4>(m1)*pow<2>(S12 - S35 + T14))*pow(T - pow<2>(m1),-1)*pow(-S12 - T - T14 + 2*pow<2>(m1),-1))*pow(-S12 + S35 - T14 + pow<2>(m2) - pow<2>(MassGlu),-2)))/9.;
}
