double MRSSM::matrixHard_gu_suLsuRubar_DS(double Alfas, std::array<std::array<double, 4>, 5> const& p) const {
   double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   double k45 = p[3][0]*p[4][0]-p[3][1]*p[4][1]-p[3][2]*p[4][2]-p[3][3]*p[4][3];
   double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   double S12 = 2*k12;
   const double m1 = MassSq;
   const double m2 = MassSq;
   const int Theta = (MassGlu > m2 ? 1 : 0) * (S12 > pow<2>(MassGlu + m1) ? 1 : 0) ;
   double S35 = m1*m1 + 2*k35;
   double S45 = m2*m2 + 2*k45;
   double T13 = m1*m1 - 2*k13;
   double T14 = m2*m2 - 2*k14;
   return 0.;
}
