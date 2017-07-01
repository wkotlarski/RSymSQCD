inline double Process::matrixMRSSMHard_gubar_suLsuLdaggerubar_DS_CSub1( const std::vector< double* >& p ) const noexcept {
   const double Alfas = pdf->alphasQ( mu_r );
   const double Alfas2 = pow(Alfas, 2);
   const double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   const double S12 = 2.*k12;
   const double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   const double k45 = p[3][0]*p[4][0]-p[3][1]*p[4][1]-p[3][2]*p[4][2]-p[3][3]*p[4][3];
   const double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   const double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   const double S35 = m1*m1 + 2*k35;
   const double S45 = m2*m2 + 2*k45;
   const double T13 = m1*m1 - 2.*k13;
   const double T14 = m2*m2 - 2.*k14;
   return (32*Alfas*Alfas2*pow(pi,3)*pow(S12,-1)*(-(S12*T14) + (S12 - 8*T14)*pow(m1,2) + 4*pow(m1,4) + 4*pow(S12,2) + 4*pow(T14,2))*(pow(m1,8) - pow(m1,6)*(S12 + T14 + pow(MassGlu,2)) + T14*pow(MassGlu,2)*(-((S12 - S45 + T13)*(S12 + T14)) + (S12 - 2*S45 + T13)*pow(MassGlu,2)) + pow(m1,4)*((2*S12 - S45 - T14)*pow(MassGlu,2) + pow(MassGlu,4) + pow(T14,2)) - pow(m1,2)*((S12 - 2*S45 + T13 + T14)*pow(MassGlu,4) + pow(MassGlu,2)*(S12*(S12 - S45 + T13) - (S12 + T13)*T14 - 2*pow(T14,2)) + (S12 + T14)*pow(T14,2)))*pow(S12 + T14 - pow(m1,2),-2)*pow(-T14 + pow(m1,2),-2)*pow(pow(MassGlu,2)*pow(WidthGlu,2) + pow(-S35 + pow(MassGlu,2),2),-1))/9.;
}
