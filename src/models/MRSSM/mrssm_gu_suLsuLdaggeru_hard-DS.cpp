// unsyplified ME and fullsimplified subtraction term
double MRSSM::matrixHard_gu_suLsuLdaggeru_DS(double Alfas, std::array<std::array<double, 4>, 5> const& p) const {
   const double m1 = MassSq;
   const double m2 = MassSq;
   const double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   const double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   const double k45 = p[3][0]*p[4][0]-p[3][1]*p[4][1]-p[3][2]*p[4][2]-p[3][3]*p[4][3];
   const double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   const double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   const double S12 = 2.*k12;
   const int Theta = ((S12 > pow<2>(MassGlu + m1)) ? 1 : 0) * (MassGlu > m2 ? 1 : 0);
   const double S35 = m1*m1 + 2*k35;
   const double S45 = m2*m2 + 2*k45;
   const double T13 = m1*m1 - 2*k13;
   const double T14 = m2*m2 - 2*k14;
   return 0.;
}
