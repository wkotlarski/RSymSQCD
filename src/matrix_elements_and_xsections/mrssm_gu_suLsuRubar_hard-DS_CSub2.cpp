double Process::matrixMRSSMHard_gu_suLsuRubar_DS_CSub2( const std::vector< double* >& p ) const noexcept {
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
   double T13 = m1*m1 - 2.*k13;
   double T14 = m2*m2 - 2.*k14;
   if (S12 > pow (MassGlu + m1,2))
      return (-32*Alfas*Alfas2*pow(pi,3)*pow(S12,-1)*(-(S12*T13) + (S12 - 8*T13)*pow(m1,2) + 4*pow(m1,4) + 4*pow(S12,2) + 4*pow(T13,2))*(-(pow(m1,6)*pow(m2,2)) + pow(m1,4)*((-S12 + S35)*pow(MassGlu,2) + pow(m2,2)*(S12 + T13 + pow(MassGlu,2))) + T13*((S12 + T13)*(S12 - S35 + T14)*pow(MassGlu,2) + pow(m2,2)*(-T13 + pow(MassGlu,2))*(-S12 - T13 + pow(MassGlu,2)) - (S12 - 2*S35 + T14)*pow(MassGlu,4)) + pow(m1,2)*(pow(MassGlu,2)*(-(S12*S35) + S12*T14 - T13*T14 + (S12 - 2*S35 + T14)*pow(MassGlu,2) + pow(S12,2)) - pow(m2,2)*((S12 - T13)*pow(MassGlu,2) + pow(MassGlu,4) + pow(T13,2))))*pow(S12 + T13 - pow(m1,2),-2)*pow(-T13 + pow(m1,2),-2)*pow(pow(MassGlu,2)*pow(WidthGlu,2) + pow(-S45 + pow(MassGlu,2),2),-1))/9.;
   else return 0.;
}
