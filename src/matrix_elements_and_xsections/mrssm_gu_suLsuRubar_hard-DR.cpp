double Process::matrixMRSSMHard_gu_suLsuRubar_DR(std::vector<double*> const& p) const {
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
   return (-16*Alfas*Alfas2*pow(pi,3)*pow(S12,-1)*(-((8*(8*S12*S35*S45*T + 7*S12*S35*S45*T14 - 8*S12*S35*T*T14 - 4*S12*S45*T*T14 + 5*S35*S45*T*T14 - 6*S12*S35*S45*pow(m1,2) + 8*S12*S35*T*pow(m1,2) + S12*S45*T*pow(m1,2) - 4*S35*S45*T*pow(m1,2) + 9*S12*S35*T14*pow(m1,2) + 3*S12*S45*T14*pow(m1,2) - 3*S35*S45*T14*pow(m1,2) - 5*S12*T*T14*pow(m1,2) + 7*S35*T*T14*pow(m1,2) + 2*S45*T*T14*pow(m1,2) - 3*S12*S35*pow(m1,4) + S35*S45*pow(m1,4) - S35*T*pow(m1,4) + 3*S12*T14*pow(m1,4) - 4*S35*T14*pow(m1,4) - S45*T14*pow(m1,4) + T*T14*pow(m1,4) - 5*S12*S35*S45*pow(m2,2) + 4*S12*S35*T*pow(m2,2) + 5*S12*S45*T*pow(m2,2) - 3*S35*S45*T*pow(m2,2) + 2*S12*S35*T14*pow(m2,2) + 4*S12*S45*T14*pow(m2,2) - 2*S35*S45*T14*pow(m2,2) - 3*S12*T*T14*pow(m2,2) + 2*S35*T*T14*pow(m2,2) + 3*S45*T*T14*pow(m2,2) - 3*S12*S35*pow(m1,2)*pow(m2,2) - 2*S12*S45*pow(m1,2)*pow(m2,2) + S35*S45*pow(m1,2)*pow(m2,2) + 2*S12*T*pow(m1,2)*pow(m2,2) - S35*T*pow(m1,2)*pow(m2,2) - 2*S45*T*pow(m1,2)*pow(m2,2) + 3*S12*T14*pow(m1,2)*pow(m2,2) - 4*S35*T14*pow(m1,2)*pow(m2,2) - S45*T14*pow(m1,2)*pow(m2,2) + T*T14*pow(m1,2)*pow(m2,2) - 2*S12*S45*pow(m2,4) + 2*S12*T*pow(m2,4) - 2*S45*T*pow(m2,4) + 5*S35*S45*pow(S12,2) - 7*S35*T*pow(S12,2) - 3*S45*T*pow(S12,2) - 5*S35*T14*pow(S12,2) - 4*S45*T14*pow(S12,2) + 3*T*T14*pow(S12,2) + 6*S35*pow(m1,2)*pow(S12,2) + 2*S45*pow(m1,2)*pow(S12,2) - 2*T*pow(m1,2)*pow(S12,2) - 5*T14*pow(m1,2)*pow(S12,2) + pow(m1,4)*pow(S12,2) + 3*S35*pow(m2,2)*pow(S12,2) + 4*S45*pow(m2,2)*pow(S12,2) - 4*T*pow(m2,2)*pow(S12,2) - 2*T14*pow(m2,2)*pow(S12,2) + 2*pow(m1,2)*pow(m2,2)*pow(S12,2) + pow(m2,4)*pow(S12,2) - 3*S35*pow(S12,3) - 2*S45*pow(S12,3) + 2*T*pow(S12,3) + 2*T14*pow(S12,3) - 2*pow(m1,2)*pow(S12,3) - 2*pow(m2,2)*pow(S12,3) + pow(S12,4) - 2*S12*S45*pow(S35,2) + 4*S12*T*pow(S35,2) - 2*S45*T*pow(S35,2) + 2*S12*T14*pow(S35,2) - 2*S45*T14*pow(S35,2) + 2*T*T14*pow(S35,2) - 4*S12*pow(m1,2)*pow(S35,2) + 2*S45*pow(m1,2)*pow(S35,2) - 4*T*pow(m1,2)*pow(S35,2) - 2*T14*pow(m1,2)*pow(S35,2) + 2*pow(m1,4)*pow(S35,2) - 2*S12*pow(m2,2)*pow(S35,2) + 2*S45*pow(m2,2)*pow(S35,2) - 2*T*pow(m2,2)*pow(S35,2) + 2*pow(m1,2)*pow(m2,2)*pow(S35,2) + 2*pow(S12,2)*pow(S35,2) - 2*S12*S35*pow(S45,2) + S12*T*pow(S45,2) - 2*S35*T*pow(S45,2) + 2*S12*T14*pow(S45,2) - 2*S35*T14*pow(S45,2) + T*T14*pow(S45,2) - S12*pow(m1,2)*pow(S45,2) + 2*S35*pow(m1,2)*pow(S45,2) - T14*pow(m1,2)*pow(S45,2) - 2*S12*pow(m2,2)*pow(S45,2) + 2*S35*pow(m2,2)*pow(S45,2) - T*pow(m2,2)*pow(S45,2) - 2*T14*pow(m2,2)*pow(S45,2) + pow(m1,2)*pow(m2,2)*pow(S45,2) + pow(m2,4)*pow(S45,2) + pow(S12,2)*pow(S45,2) - 5*S12*S35*pow(T,2) - S12*S45*pow(T,2) + 3*S35*S45*pow(T,2) + S12*T14*pow(T,2) - 3*S35*T14*pow(T,2) - S45*T14*pow(T,2) + 2*S35*pow(m1,2)*pow(T,2) - T14*pow(m1,2)*pow(T,2) - 3*S12*pow(m2,2)*pow(T,2) + S35*pow(m2,2)*pow(T,2) + 2*S45*pow(m2,2)*pow(T,2) - T14*pow(m2,2)*pow(T,2) + pow(m1,2)*pow(m2,2)*pow(T,2) + pow(m2,4)*pow(T,2) + pow(S12,2)*pow(T,2) + 2*pow(S35,2)*pow(T,2) - S35*pow(T,3) - pow(m2,2)*pow(T,3) - 2*S12*S35*pow(T14,2) - 2*S12*S45*pow(T14,2) + 2*S35*S45*pow(T14,2) + S12*T*pow(T14,2) - 2*S35*T*pow(T14,2) - S45*T*pow(T14,2) - 5*S12*pow(m1,2)*pow(T14,2) + 4*S35*pow(m1,2)*pow(T14,2) + S45*pow(m1,2)*pow(T14,2) - 3*T*pow(m1,2)*pow(T14,2) + 2*pow(m1,4)*pow(T14,2) + 2*pow(m1,2)*pow(m2,2)*pow(T14,2) + pow(S12,2)*pow(T14,2) + pow(S45,2)*pow(T14,2) - 2*pow(m1,2)*pow(T14,3))*pow(-S12 - T - T14 + 2*pow(m1,2),-2) + 16*(-(S35*(S12 - S45 + T)) + (S35 - T14)*pow(m1,2))*(T14*(S12 - S35 + T14) + (S35 - T14)*pow(m2,2))*pow(T14 - pow(m2,2),-2) - (pow(m1,2)*((S35 - T14)*(T14*(-S45 + 3*T + 4*T14) + 2*S35*(S45 - 2*(T + T14))) + (6*S35 + 2*S45 - 2*T - 5*T14)*pow(S12,2) + pow(m2,2)*((4*S35 + S45 - T - 4*T14)*(S35 - T14) + S12*(-3*S35 - 2*S45 + 2*T + 3*T14) + 2*pow(S12,2)) - 2*pow(S12,3) + S12*((3*S45 - 4*T - 7*T14)*T14 + S35*(-5*S45 + 6*T + 11*T14) - 4*pow(S35,2))) + pow(m1,4)*(-3*S12*(S35 - T14) + pow(S12,2) + 2*pow(S35 - T14,2)) + (S12 - S45 + T)*(2*S35*S45*T14 - 3*S35*T*T14 + (S12 - S45 + T)*pow(m2,4) + pow(m2,2)*(S12*(3*S35 + 2*S45 - 2*T - 2*T14) - (4*S35 + 2*S45 - T)*(S35 - T14) - 2*pow(S12,2)) + (-3*S35 - S45 + T + 2*T14)*pow(S12,2) + pow(S12,3) + 2*T*pow(S35,2) + 4*T14*pow(S35,2) + S12*(S35*(2*S45 - 3*T - 7*T14) + T14*(-2*S45 + T + T14) + 2*pow(S35,2)) - 4*S35*pow(T14,2) - S45*pow(T14,2)))*pow(-S12 - T - T14 + 2*pow(m1,2),-1)*pow(T14 - pow(m2,2),-1))*pow(-S12 + S45 - T + pow(m1,2) - pow(MassGlu,2),-2)) + (16*(T*(S12 - S45 + T) + (S45 - T)*pow(m1,2))*(S45*(S12 - S35 + T14) + (-S45 + T)*pow(m2,2))*pow(T - pow(m1,2),-2) - 8*(7*S12*S35*S45*T + 8*S12*S35*S45*T14 - 4*S12*S35*T*T14 - 8*S12*S45*T*T14 + 5*S35*S45*T*T14 - 5*S12*S35*S45*pow(m1,2) + 4*S12*S35*T*pow(m1,2) + 2*S12*S45*T*pow(m1,2) - 2*S35*S45*T*pow(m1,2) + 5*S12*S35*T14*pow(m1,2) + 4*S12*S45*T14*pow(m1,2) - 3*S35*S45*T14*pow(m1,2) - 3*S12*T*T14*pow(m1,2) + 3*S35*T*T14*pow(m1,2) + 2*S45*T*T14*pow(m1,2) - 2*S12*S35*pow(m1,4) + 2*S12*T14*pow(m1,4) - 2*S35*T14*pow(m1,4) - 6*S12*S35*S45*pow(m2,2) + 3*S12*S35*T*pow(m2,2) + 9*S12*S45*T*pow(m2,2) - 3*S35*S45*T*pow(m2,2) + S12*S35*T14*pow(m2,2) + 8*S12*S45*T14*pow(m2,2) - 4*S35*S45*T14*pow(m2,2) - 5*S12*T*T14*pow(m2,2) + 2*S35*T*T14*pow(m2,2) + 7*S45*T*T14*pow(m2,2) - 2*S12*S35*pow(m1,2)*pow(m2,2) - 3*S12*S45*pow(m1,2)*pow(m2,2) + S35*S45*pow(m1,2)*pow(m2,2) + 3*S12*T*pow(m1,2)*pow(m2,2) - S35*T*pow(m1,2)*pow(m2,2) - 4*S45*T*pow(m1,2)*pow(m2,2) + 2*S12*T14*pow(m1,2)*pow(m2,2) - 2*S35*T14*pow(m1,2)*pow(m2,2) - S45*T14*pow(m1,2)*pow(m2,2) + T*T14*pow(m1,2)*pow(m2,2) - 3*S12*S45*pow(m2,4) + S35*S45*pow(m2,4) + 3*S12*T*pow(m2,4) - S35*T*pow(m2,4) - 4*S45*T*pow(m2,4) - S45*T14*pow(m2,4) + T*T14*pow(m2,4) + 5*S35*S45*pow(S12,2) - 4*S35*T*pow(S12,2) - 5*S45*T*pow(S12,2) - 3*S35*T14*pow(S12,2) - 7*S45*T14*pow(S12,2) + 3*T*T14*pow(S12,2) + 4*S35*pow(m1,2)*pow(S12,2) + 3*S45*pow(m1,2)*pow(S12,2) - 2*T*pow(m1,2)*pow(S12,2) - 4*T14*pow(m1,2)*pow(S12,2) + pow(m1,4)*pow(S12,2) + 2*S35*pow(m2,2)*pow(S12,2) + 6*S45*pow(m2,2)*pow(S12,2) - 5*T*pow(m2,2)*pow(S12,2) - 2*T14*pow(m2,2)*pow(S12,2) + 2*pow(m1,2)*pow(m2,2)*pow(S12,2) + pow(m2,4)*pow(S12,2) - 2*S35*pow(S12,3) - 3*S45*pow(S12,3) + 2*T*pow(S12,3) + 2*T14*pow(S12,3) - 2*pow(m1,2)*pow(S12,3) - 2*pow(m2,2)*pow(S12,3) + pow(S12,4) - 2*S12*S45*pow(S35,2) + 2*S12*T*pow(S35,2) - 2*S45*T*pow(S35,2) + S12*T14*pow(S35,2) - 2*S45*T14*pow(S35,2) + T*T14*pow(S35,2) - 2*S12*pow(m1,2)*pow(S35,2) + 2*S45*pow(m1,2)*pow(S35,2) - 2*T*pow(m1,2)*pow(S35,2) - T14*pow(m1,2)*pow(S35,2) + pow(m1,4)*pow(S35,2) - S12*pow(m2,2)*pow(S35,2) + 2*S45*pow(m2,2)*pow(S35,2) - T*pow(m2,2)*pow(S35,2) + pow(m1,2)*pow(m2,2)*pow(S35,2) + pow(S12,2)*pow(S35,2) - 2*S12*S35*pow(S45,2) + 2*S12*T*pow(S45,2) - 2*S35*T*pow(S45,2) + 4*S12*T14*pow(S45,2) - 2*S35*T14*pow(S45,2) + 2*T*T14*pow(S45,2) - 2*S12*pow(m1,2)*pow(S45,2) + 2*S35*pow(m1,2)*pow(S45,2) - 2*T14*pow(m1,2)*pow(S45,2) - 4*S12*pow(m2,2)*pow(S45,2) + 2*S35*pow(m2,2)*pow(S45,2) - 2*T*pow(m2,2)*pow(S45,2) - 4*T14*pow(m2,2)*pow(S45,2) + 2*pow(m1,2)*pow(m2,2)*pow(S45,2) + 2*pow(m2,4)*pow(S45,2) + 2*pow(S12,2)*pow(S45,2) - 2*S12*S35*pow(T,2) - 2*S12*S45*pow(T,2) + 2*S35*S45*pow(T,2) + S12*T14*pow(T,2) - S35*T14*pow(T,2) - 2*S45*T14*pow(T,2) - 5*S12*pow(m2,2)*pow(T,2) + S35*pow(m2,2)*pow(T,2) + 4*S45*pow(m2,2)*pow(T,2) - 3*T14*pow(m2,2)*pow(T,2) + 2*pow(m1,2)*pow(m2,2)*pow(T,2) + 2*pow(m2,4)*pow(T,2) + pow(S12,2)*pow(T,2) + pow(S35,2)*pow(T,2) - 2*pow(m2,2)*pow(T,3) - S12*S35*pow(T14,2) - 5*S12*S45*pow(T14,2) + 3*S35*S45*pow(T14,2) + S12*T*pow(T14,2) - S35*T*pow(T14,2) - 3*S45*T*pow(T14,2) - 3*S12*pow(m1,2)*pow(T14,2) + 2*S35*pow(m1,2)*pow(T14,2) + S45*pow(m1,2)*pow(T14,2) - T*pow(m1,2)*pow(T14,2) + pow(m1,4)*pow(T14,2) + 2*S45*pow(m2,2)*pow(T14,2) - T*pow(m2,2)*pow(T14,2) + pow(m1,2)*pow(m2,2)*pow(T14,2) + pow(S12,2)*pow(T14,2) + 2*pow(S45,2)*pow(T14,2) - S45*pow(T14,3) - pow(m1,2)*pow(T14,3))*pow(-S12 - T - T14 + 2*pow(m1,2),-2) + (9*S12*S35*S45*T + 5*S12*S35*S45*T14 - 3*S12*S35*T*T14 - 10*S12*S45*T*T14 + 5*S35*S45*T*T14 + 5*S35*S45*pow(S12,2) - 4*S35*T*pow(S12,2) - 7*S45*T*pow(S12,2) - 2*S35*T14*pow(S12,2) - 6*S45*T14*pow(S12,2) + 3*T*T14*pow(S12,2) + pow(m1,2)*(pow(m2,2)*((S45 - T)*(S35 + 4*S45 - 4*T - T14) + S12*(-2*S35 - 3*S45 + 3*T + 2*T14) + 2*pow(S12,2)) - (S12 - S35 + T14)*((S45 - T)*(2*S35 + 4*S45 - T14) + S12*(-2*S35 - 3*S45 + 2*(T + T14)) + 2*pow(S12,2))) - 2*S35*pow(S12,3) - 3*S45*pow(S12,3) + 2*T*pow(S12,3) + 2*T14*pow(S12,3) + pow(S12,4) - 2*S12*S45*pow(S35,2) + 2*S12*T*pow(S35,2) - 2*S45*T*pow(S35,2) + pow(S12,2)*pow(S35,2) - 2*S12*S35*pow(S45,2) + 4*S12*T*pow(S45,2) - 4*S35*T*pow(S45,2) + 4*S12*T14*pow(S45,2) - 2*S35*T14*pow(S45,2) + 4*T*T14*pow(S45,2) + 2*pow(S12,2)*pow(S45,2) + pow(m2,4)*(-3*S12*(S45 - T) + pow(S12,2) + 2*pow(S45 - T,2)) - 2*S12*S35*pow(T,2) - 4*S12*S45*pow(T,2) + 4*S35*S45*pow(T,2) + S12*T14*pow(T,2) - S35*T14*pow(T,2) - 4*S45*T14*pow(T,2) + pow(S12,2)*pow(T,2) + pow(S35,2)*pow(T,2) - pow(m2,2)*((S45 - T)*(S35*(-2*S45 + T) + 4*S45*(T + T14) - T*(4*T + 3*T14)) + (-2*S35 - 6*S45 + 5*T + 2*T14)*pow(S12,2) + 2*pow(S12,3) + S12*(5*S35*S45 - 3*S35*T - 11*S45*T - 6*S45*T14 + 4*T*T14 + 4*pow(S45,2) + 7*pow(T,2))) - 3*S12*S45*pow(T14,2) + S12*T*pow(T14,2) - 3*S45*T*pow(T14,2) + pow(S12,2)*pow(T14,2) + 2*pow(S45,2)*pow(T14,2) + pow(m1,4)*pow(S12 - S35 + T14,2))*pow(T - pow(m1,2),-1)*pow(-S12 - T - T14 + 2*pow(m1,2),-1))*pow(-S12 + S35 - T14 + pow(m2,2) - pow(MassGlu,2),-2)))/9.;
}
