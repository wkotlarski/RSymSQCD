double MSSM::matrixVirt_ddbar_suLsuLdagger(const double alphas, const double s, const double t,
   const double FiniteGs, const double Dminus4, const int divergence, const double mu) const {
   const double u = 2*pow<2>(MassSq) - s - t;
   const double MassGlu2 = pow<2>(MassGlu);
   const double MassSq2 = pow<2>(MassSq);
   const double MassTop2 = pow<2>(MassTop);
   const double mu2 = pow<2>(mu);
   const double alphas3 = pow<3>(alphas);
   setmudim(mu2);
   setuvdiv(0.);
   setlambda(divergence);
   std::complex<double> matrix = (0.296296296296296296*alphas3*pi*(Power(MassSq,4) - 1.*t*u)*(-72.*A0(MassGlu2) + 144.*A0(MassSq2) - 24.*A0(MassTop2) + 28.*s*B0i(bb0,MassSq2,0,MassSq2) - 154.*s*B0i(bb0,s,0,0) + 2.*Dminus4*s*B0i(bb0,s,0,0) - 72.*s*B0i(bb0,s,MassGlu2,MassGlu2) + 96.*B0i(bb00,s,0,0) - 72.*Dminus4*B0i(bb00,s,0,0) + 144.*B0i(bb00,s,MassGlu2,MassGlu2) - 288.*B0i(bb00,s,MassSq2,MassSq2) + 48.*B0i(bb00,s,MassTop2,MassTop2) + 14.*s*B0i(bb1,MassSq2,0,MassSq2) - 156.*s*B0i(bb1,s,0,0) - 72.*s*B0i(bb1,s,MassGlu2,MassGlu2) - 24.*s*B0i(bb1,s,MassTop2,MassTop2) - 4.*(s*s)*C0i(cc0,0,s,0,0,0,0) + 8.*MassSq2*s*C0i(cc0,MassSq2,s,MassSq2,0,MassSq2,MassSq2) - 4.*(s*s)*C0i(cc0,MassSq2,s,MassSq2,0,MassSq2,MassSq2) + 8.*MassGlu2*s*C0i(cc0,MassSq2,s,MassSq2,MassGlu2,0,0) - 36.*(s*s)*C0i(cc0,s,0,0,0,0,0) + 36.*MassGlu2*s*C0i(cc0,s,0,0,MassGlu2,MassGlu2,MassSq2) - 36.*MassSq2*s*C0i(cc0,s,0,0,MassGlu2,MassGlu2,MassSq2) - 72.*(s*s)*C0i(cc0,s,MassSq2,MassSq2,0,0,MassSq2) - 18.*(s*s)*C0i(cc0,s,MassSq2,MassSq2,MassGlu2,MassGlu2,0) - 8.*s*C0i(cc00,0,s,0,0,0,0) - 4.*Dminus4*s*C0i(cc00,0,s,0,0,0,0) + 8.*s*C0i(cc00,0,s,0,MassGlu2,MassSq2,MassSq2) + 8.*s*C0i(cc00,MassSq2,s,MassSq2,0,MassSq2,MassSq2) - 72.*s*C0i(cc00,s,0,0,0,0,0) - 36.*Dminus4*s*C0i(cc00,s,0,0,0,0,0) + 72.*s*C0i(cc00,s,0,0,MassGlu2,MassGlu2,MassSq2) + 36.*s*C0i(cc00,s,MassSq2,MassSq2,0,0,MassSq2) - 4.*(s*s)*C0i(cc1,0,s,0,0,0,0) - 36.*MassGlu2*s*C0i(cc1,MassSq2,s,MassSq2,0,MassGlu2,MassGlu2) - 36.*MassSq2*s*C0i(cc1,MassSq2,s,MassSq2,0,MassGlu2,MassGlu2) + 16.*MassSq2*s*C0i(cc1,MassSq2,s,MassSq2,0,MassSq2,MassSq2) - 6.*(s*s)*C0i(cc1,MassSq2,s,MassSq2,0,MassSq2,MassSq2) + 4.*MassGlu2*s*C0i(cc1,MassSq2,s,MassSq2,MassGlu2,0,0) + 4.*MassSq2*s*C0i(cc1,MassSq2,s,MassSq2,MassGlu2,0,0) + 8.*MassSq2*s*C0i(cc11,MassSq2,s,MassSq2,0,MassSq2,MassSq2) - 2.*(s*s)*C0i(cc11,MassSq2,s,MassSq2,0,MassSq2,MassSq2) + 16.*MassSq2*s*C0i(cc12,MassSq2,s,MassSq2,0,MassSq2,MassSq2) - 4.*(s*s)*C0i(cc12,MassSq2,s,MassSq2,0,MassSq2,MassSq2) - 4.*(s*s)*C0i(cc2,0,s,0,0,0,0) - 36.*MassGlu2*s*C0i(cc2,MassSq2,s,MassSq2,0,MassGlu2,MassGlu2) - 36.*MassSq2*s*C0i(cc2,MassSq2,s,MassSq2,0,MassGlu2,MassGlu2) + 16.*MassSq2*s*C0i(cc2,MassSq2,s,MassSq2,0,MassSq2,MassSq2) - 6.*(s*s)*C0i(cc2,MassSq2,s,MassSq2,0,MassSq2,MassSq2) + 4.*MassGlu2*s*C0i(cc2,MassSq2,s,MassSq2,MassGlu2,0,0) + 4.*MassSq2*s*C0i(cc2,MassSq2,s,MassSq2,MassGlu2,0,0) - 36.*(s*s)*C0i(cc2,s,0,0,0,0,0) + 72.*MassSq2*s*C0i(cc2,s,MassSq2,MassSq2,0,0,MassSq2) - 45.*(s*s)*C0i(cc2,s,MassSq2,MassSq2,0,0,MassSq2) + 18.*(s*s)*C0i(cc2,s,MassSq2,MassSq2,MassGlu2,MassGlu2,0) + 8.*MassSq2*s*C0i(cc22,MassSq2,s,MassSq2,0,MassSq2,MassSq2) - 2.*(s*s)*C0i(cc22,MassSq2,s,MassSq2,0,MassSq2,MassSq2) + 36.*MassSq2*s*C0i(cc22,s,MassSq2,MassSq2,0,0,MassSq2) - 9.*(s*s)*C0i(cc22,s,MassSq2,MassSq2,0,0,MassSq2) - 28.*MassSq2*(s*s)*D0i(dd0,s,0,t,MassSq2,0,MassSq2,0,0,0,MassSq2) + 28.*(s*s)*t*D0i(dd0,s,0,t,MassSq2,0,MassSq2,0,0,0,MassSq2) - 8.*MassSq2*(s*s)*D0i(dd0,s,0,u,MassSq2,0,MassSq2,0,0,0,MassSq2) + 8.*(s*s)*u*D0i(dd0,s,0,u,MassSq2,0,MassSq2,0,0,0,MassSq2) + 14.*MassGlu2*(s*s)*D0i(dd0,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) - 14.*MassSq2*(s*s)*D0i(dd0,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 4.*MassGlu2*(s*s)*D0i(dd0,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) - 4.*MassSq2*(s*s)*D0i(dd0,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 28.*(s*s)*D0i(dd00,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 8.*(s*s)*D0i(dd00,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) - 28.*MassSq2*(s*s)*D0i(dd2,s,0,t,MassSq2,0,MassSq2,0,0,0,MassSq2) + 28.*(s*s)*t*D0i(dd2,s,0,t,MassSq2,0,MassSq2,0,0,0,MassSq2) - 8.*MassSq2*(s*s)*D0i(dd2,s,0,u,MassSq2,0,MassSq2,0,0,0,MassSq2) + 8.*(s*s)*u*D0i(dd2,s,0,u,MassSq2,0,MassSq2,0,0,0,MassSq2) + 14.*MassGlu2*(s*s)*D0i(dd2,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 28.*MassSq2*(s*s)*D0i(dd2,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) - 14.*(s*s)*t*D0i(dd2,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 4.*MassGlu2*(s*s)*D0i(dd2,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 8.*MassSq2*(s*s)*D0i(dd2,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) - 4.*(s*s)*u*D0i(dd2,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 14.*MassSq2*(s*s)*D0i(dd22,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 14.*(s*s)*t*D0i(dd22,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 4.*MassSq2*(s*s)*D0i(dd22,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 4.*(s*s)*u*D0i(dd22,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 14.*MassSq2*(s*s)*D0i(dd23,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) - 14.*(s*s)*t*D0i(dd23,s,MassSq2,t,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) + 4.*MassSq2*(s*s)*D0i(dd23,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) - 4.*(s*s)*u*D0i(dd23,s,MassSq2,u,0,MassSq2,0,MassGlu2,MassGlu2,0,MassSq2) - 28.*MassSq2*(s*s)*D0i(dd3,s,0,t,MassSq2,0,MassSq2,0,0,0,MassSq2) - 8.*MassSq2*(s*s)*D0i(dd3,s,0,u,MassSq2,0,MassSq2,0,0,0,MassSq2) + 24.*FiniteGs*s*log(MassGlu2/mu2) + 24.*FiniteGs*s*log(MassSq2/mu2) + 8.*FiniteGs*s*log(MassTop2/mu2) + 32.*s*Re(B0i(bb0,0,0,0)) + 16.*Dminus4*s*Re(B0i(bb0,0,0,0)) + 64.*s*Re(B0i(bb0,MassSq2,0,MassGlu2)) - 48.*s*Re(B0i(bb0,MassSq2,0,MassSq2)) + 32.*s*Re(B0i(bb1,0,0,0)) + 16.*Dminus4*s*Re(B0i(bb1,0,0,0)) - 32.*s*Re(B0i(bb1,0,MassGlu2,MassSq2)) - 32.*s*Re(B0i(bb1,MassSq2,0,MassSq2)) + 64.*s*Re(B0i(bb1,MassSq2,MassGlu2,0)) + 64.*MassSq2*s*Re(B0i(dbb0,MassSq2,0,MassGlu2)) - 64.*MassSq2*s*Re(B0i(dbb0,MassSq2,0,MassSq2)) - 32.*MassSq2*s*Re(B0i(dbb1,MassSq2,0,MassSq2)) + 64.*MassSq2*s*Re(B0i(dbb1,MassSq2,MassGlu2,0))))/Power(s,3);
   return matrix.real();
}