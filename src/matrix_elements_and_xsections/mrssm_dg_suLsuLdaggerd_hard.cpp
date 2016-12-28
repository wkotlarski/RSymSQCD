double Process::matrixMRSSMHard_dg_suLsuLdaggerd( std::vector< double* >& p ) {
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow(Alfas, 2);
   double k12 = p[0][0]*p[1][0]-p[0][1]*p[1][1]-p[0][2]*p[1][2]-p[0][3]*p[1][3];
   double k13 = p[0][0]*p[2][0]-p[0][1]*p[2][1]-p[0][2]*p[2][2]-p[0][3]*p[2][3];
   double k14 = p[0][0]*p[3][0]-p[0][1]*p[3][1]-p[0][2]*p[3][2]-p[0][3]*p[3][3];
   double k23 = p[1][0]*p[2][0]-p[1][1]*p[2][1]-p[1][2]*p[2][2]-p[1][3]*p[2][3];
   double k24 = p[1][0]*p[3][0]-p[1][1]*p[3][1]-p[1][2]*p[3][2]-p[1][3]*p[3][3];
   double k34 = p[2][0]*p[3][0]-p[2][1]*p[3][1]-p[2][2]*p[3][2]-p[2][3]*p[3][3];
   double S12 = 2.*k12;
   double T = m1*m1 - 2.*k13;
   double U = m1*m1 - 2.*k23;
   double S34 = m1*m1 + m2*m2 + 2.*k34;
   double T14 = m2*m2 - 2.*k14;
   double T24 = m2*m2 - 2.*k24;
   return (-8*Alfas*Alfas2*pow(pi,3)*(-8*((S12 + 2*(T + U))*(S12 + 2*(S34 + T + U)) - 4*(3*S12 + 2*S34 + 5*T + T14 + T24 + 5*U)*pow(m1,2) + 32*pow(m1,4))*pow(S12,-1)*pow(S34,-2) + 8*pow(S34,-2)*((S34 + T - T14)*(S12*(S34 + 3*T + T14) + 2*T*(2*(S34 + T + T14) + T24)) - 2*(2*S12*T + T14*(S34 + 3*T + T14) + 2*T*T24)*U - 4*(5*S12 + 6*S34 + 27*T + 7*T14 - 5*T24 + 7*U)*pow(m1,4) + 32*pow(m1,6) - 4*T*pow(U,2) + 2*pow(m1,2)*(20*T*T14 + 11*T*T24 - T14*T24 + 23*T*U + 7*T14*U - 2*T24*U + 2*S12*(S34 + 7*T + 2*T14 - T24 + 2*U) + S34*(7*T + 5*T14 - T24 + 5*U) + 2*pow(S12,2) + 2*pow(S34,2) + 5*pow(T,2) + 3*pow(T14,2) - 4*pow(T24,2) + 2*pow(U,2)))*pow(S12 + T24 + U - 2*pow(m1,2),-2) + pow(S34,-1)*pow(S12 + T + T14 - 2*pow(m1,2),-1)*pow(S12 + T24 + U - 2*pow(m1,2),-1)*(5*(S12*(T + 2*T14) + T*(3*(S34 + T + T14) + T24) + (3*T + 2*T14)*U - 2*(2*S12 + 2*S34 + 7*T + 3*T14 + T24 + 3*U)*pow(m1,2) + 16*pow(m1,4)) + 2*(9*pow(S34,-1)*(-((T - T14 - T24 + U)*(-(T14*T24) + T*U)) - S34*(T*(4*T14 + 2*T24 + U) + T14*(T24 + 2*U)) + (3*S12 - 11*S34 + 23*T - 17*T14 + 31*T24 - 17*U)*pow(m1,4) + 28*pow(m1,6) - 2*(T + T14 + T24 + U)*pow(S34,2) - 2*pow(S34,3) - pow(m1,2)*(8*T*T14 + 18*T*T24 + 4*T14*T24 + S12*(7*S34 + 2*(6*T + T14 + 5*T24 - U)) + 2*(T - 2*T14 + 5*T24)*U + S34*(4*T - 6*(T14 - T24 + U)) - 9*pow(S34,2) + 7*pow(T,2) - 3*pow(T14,2) + 9*pow(T24,2) - 3*pow(U,2))) - 7*(-(S12*S34*(S34 + 3*T + T14)) + (S34*(-4*T + T14) + T14*(-T + T14))*U - (8*S12 + 3*S34 + 19*T + T14 + 2*(T24 + U))*pow(m1,4) + 12*pow(m1,6) - 2*T*((T - T14)*T14 + S34*(2*T + T14 + T24) + 2*pow(S34,2)) + pow(m1,2)*(5*T*T14 + 2*S12*(2*S34 + 3*T + T14) + 4*T*T24 - 2*T14*T24 + S34*(18*T - T14 - 2*T24 - U) + 5*T*U - 3*T14*U + 4*pow(T,2) - pow(T14,2)))*pow(T24 - pow(m1,2),-1) + 2*((S34 + T - T14)*(S12*(S34 + 2*T) + T*(2*(S34 + T) + T24)) - 2*(S34*T14 + T*(S12 + 2*T14 + T24))*U - (17*S12 + 10*S34 + 20*T + 14*T14 + 7*T24 + 15*U)*pow(m1,4) + 24*pow(m1,6) - 2*T*pow(U,2) + pow(m1,2)*(-3*S34*T + 5*S34*T14 + 12*T*T14 + 2*S34*T24 + 6*T*T24 + T14*T24 + (5*S34 + 15*T + 4*T14 + 2*T24)*U + S12*(S34 + 5*T + 4*(T14 + T24) + 6*U) + 4*pow(S12,2) + pow(S34,2) - 4*pow(T,2) + 2*pow(T14,2) + 2*pow(U,2)))*pow(U - pow(m1,2),-1))) + 2*pow(S12,-1)*pow(S34,-1)*(-(pow(S34,-1)*(-4*T*(S34 + T)*(T + T14) - 2*(S34*(T + T14) + T*(3*T + 2*T14))*U - 2*(3*S12 + 2*S34 + 19*T + 3*T14 - T24 + 7*U)*pow(m1,4) + 16*pow(m1,6) - (T + T14)*pow(S12,2) + S12*(-(S34*(T + T14)) - T*(4*(T + T14) + T24) - (2*T + 3*T14)*U + pow(S34,2)) - 2*T24*pow(T,2) - 2*(T + T14)*pow(U,2) + pow(m1,2)*(S12*(-2*S34 + 11*T + 3*T14 + T24 + 5*U) + 2*pow(S12,2) + 2*(-(T14*T24) + 3*T14*U + S34*(5*T + T14 + 2*U) + T*(6*T14 + T24 + 9*U) + 12*pow(T,2) + 2*pow(U,2))))*pow(S12 + T24 + U - 2*pow(m1,2),-1)) - pow(S12 + T + T14 - 2*pow(m1,2),-1)*(-5*S12*(S12 + S34 + 2*(T + U) - 4*pow(m1,2)) + 9*pow(S34,-1)*(-4*T*(S34 + T)*(T + T14) - 2*(S34*(T + T14) + T*(3*T + 2*T14))*U - (23*S12 + 5*S34 + 39*T + 7*T14 - T24 + 15*U)*pow(m1,4) + 20*pow(m1,6) + (-2*T + T24 - U)*pow(S12,2) - 2*T24*pow(T,2) + S12*(-2*S34*T - 2*T*T14 + S34*T24 + T*T24 - (S34 + 6*T + T14 - 2*T24)*U + 2*pow(S34,2) - 6*pow(T,2) - 2*pow(U,2)) - 2*(T + T14)*pow(U,2) + pow(m1,2)*(S12*(-2*S34 + 19*T + 3*T14 + T24 + 13*U) + 6*pow(S12,2) + 2*(-(T14*T24) + 3*T14*U + S34*(5*T + T14 + 2*U) + T*(6*T14 + T24 + 9*U) + 12*pow(T,2) + 2*pow(U,2)))) + 2*((S12 + 2*(S34 + T + U))*(S12*S34 - T14*(2*T + U)) + (9*S12 - 2*(S34 + 7*T + 4*T14 - T24 + 4*U))*pow(m1,4) + 8*pow(m1,6) + pow(m1,2)*(S12*(-10*S34 - 4*T + T14 - 3*U) - 2*pow(S12,2) + 2*(7*T*T14 - T14*T24 + 3*T*U + 4*T14*U + S34*(2*T + T14 + U) + 2*pow(T,2) + pow(U,2))))*pow(T24 - pow(m1,2),-1) + 7*(-((S12 + 2*T)*(S12*(S34 + 2*T) + T*(2*(S34 + T) + T24))) - 2*(S12*(S34 + 2*T) + T*(S34 + 3*T - T14))*U - (31*S12 + 2*S34 + 24*T - 2*T14 + 6*U)*pow(m1,4) + 8*pow(m1,6) - 2*T*pow(U,2) + pow(m1,2)*(S12*(8*S34 + 23*T + 2*T14 + 5*T24 + 12*U) + 8*pow(S12,2) + 2*(T*(3*S34 + 10*T - T14 + T24) + (S34 + 6*T - T14)*U + pow(U,2))))*pow(U - pow(m1,2),-1))) + 2*pow(S12 + T + T14 - 2*pow(m1,2),-2)*(-7*S12 - 9*pow(S34,-2)*(16*(S12 + 2*S34 - T24 + U)*pow(m1,4) + S12*(-4*T*(T24 - 3*U) - 4*pow(S34,2) + 8*pow(T,2) + pow(T24 - U,2)) - 2*pow(m1,2)*(6*S34*(T + T14) - T24*(9*T + 3*T14 + T24) + (25*T + 3*T14 - 2*T24)*U + 4*S34*(T24 + 2*U) + 2*S12*(T24 + 3*U) + 4*pow(S12,2) + 4*pow(S34,2) + 16*pow(T,2) + 3*pow(U,2)) + 2*(U*(2*S34*(3*T + T14) - T14*T24 + 3*T*(2*T14 + T24) + 10*pow(T,2)) + T*(4*(S34 + T)*(T + T14) + 2*(T - T14)*T24 - pow(T24,2)) + (6*T + T14)*pow(U,2))) + pow(T24 - pow(m1,2),-1)*(-7*(S12*S34 - T14*(2*T + U) + (-2*S12 + T14 + U)*pow(m1,2) + pow(m1,4)) + 32*pow(m1,2)*(T*(S34 + T + U) + (-3*T + T24)*pow(m1,2))*pow(T24 - pow(m1,2),-1)) - (-7*(-(S12*S34) + T*(2*T14 + T24) - (2*S34 + 3*T + 2*T14 + 3*T24 + 2*U)*pow(m1,2) + 7*pow(m1,4)) + 2*(-S34 + 2*pow(m1,2))*(S12*(S34 + 2*T) + T*(2*(S34 + T) + T24) + (2*T - T14)*U - (3*S12 + S34 + 8*T)*pow(m1,2) + 4*pow(m1,4))*pow(T24 - pow(m1,2),-1))*pow(U - pow(m1,2),-1) + 9*pow(S34,-1)*((4*S34*T*(S34 + T) + 2*T*(S34 + T14)*T24 + (4*S34*T - 2*S34*T14 - 2*T*T14 + T14*T24)*U + S12*S34*(2*S34 + 4*T - T24 + U) + (8*S12 + 4*T - 4*T14 + 3*T24 - 3*U)*pow(m1,4) - T14*pow(U,2) + pow(m1,2)*(-4*T*T24 - T14*T24 + 5*T14*U - T24*U - 2*S12*(4*S34 + 4*T - T24 + U) + 2*S34*(-7*T + T14 + T24 + U) + pow(U,2)))*pow(T24 - pow(m1,2),-1) + (T*(2*S34 + T24)*(2*(S34 + T) + T24) + S12*(S34 + 2*T)*(2*S34 + T24 - U) + (2*S34*(T - T14) + T*(-2*T + T24))*U + (8*S12 - 4*T + 4*T14 + 5*T24 + 11*U)*pow(m1,4) - 2*T*pow(U,2) - pow(m1,2)*(14*S34*T - 2*S34*T14 + 5*T*T24 + 2*T14*T24 + 4*S12*(2*(S34 + T14) + T24) - 5*T*U + 6*T14*U + 5*T24*U + pow(T24,2) + 2*pow(U,2)))*pow(U - pow(m1,2),-1)) + 32*pow(m1,2)*(-(T14*(S12 + T + U)) + (T14 + U)*pow(m1,2))*pow(-U + pow(m1,2),-2))))/9.;
}