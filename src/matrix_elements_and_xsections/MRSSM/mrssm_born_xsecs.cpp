double Process::sigmaMRSSMTree_uubar_suLsuLdagger(double alphas, double s12) {
	const double Alfas2 = Sqr(alphas);
   return (Alfas2*pi*(8*(-3*Sqr(MassSq) + MassGlu*MassGlu - 2*s12)*Sqrt(s12*(-4*Sqr(MassSq) + s12)) +
       8*(2*Power4(MassSq) + 2*Power4(MassGlu) - 4*Sqr(MassGlu)*s12 - 3*(s12*s12) + Sqr(MassSq)*(-4*Sqr(MassGlu) + 6*s12))*
        atanh(Sqrt(s12*(-4*Sqr(MassSq) + s12))/(2*Sqr(MassSq) - 2*Sqr(MassGlu) - s12))))/(54.*Power3(s12));
}

// checked against MG with SUSYQCD model
double Process::sigmaMRSSMTree_ddbar_suLsuLdagger(double alphas, double s12 ) {
	const double Alfas2 = Sqr(alphas);
   return (2*Alfas2*pi*pow(-4*Sqr(MassSq) + s12,1.5))/(27.*pow(s12,2.5));
}

double Process::sigmaMRSSMTree_gg_suLsuLdagger(double alphas, double s12 ) {
	const double Alfas2 = Sqr(alphas);
   return (Alfas2*pi*(Sqrt(s12*(-4*Sqr(MassSq) + s12))*(62*Sqr(MassSq) + 5*s12) - 16*Sqr(MassSq)*(Sqr(MassSq) + 4*s12)*atanh(Sqrt(1 - (4*Sqr(MassSq))/s12))))/(48.*Power3(s12));
}

/*
 *    checked with MadGraph
 */
double Process::sigmaMRSSMTree_uu_suLsuR(double alphas, double s ) {
   double MGl2 = pow(MassGlu, 2);
   double a = Sqr(alphas);
   return (-4.*a*pi*(2.*sqrt(s*(-4*Sqr(MassSq) + s)) + (2.*Sqr(MassSq) - 2.*MGl2 - s)*
      log((4.*MGl2 + Sqr(1. + sqrt(1. - (4.*Sqr(MassSq))/s))*s)/
      (4.*MGl2 + Sqr(-1. + sqrt(1. - (4.*Sqr(MassSq))/s))*s))))/(9.*Sqr(s));
}

/*
double Process::sigmaMRSSMTree_uubar_OO(double alphas, double s ) {
   double a2 = Sqr(alphas);
   double m1 = 1000;
   double b = sqrt( 1. - 4. * m1*m1/s );
   return 2 * a2 * pi * pow(b, 3)/(9. * s);
}
*/
