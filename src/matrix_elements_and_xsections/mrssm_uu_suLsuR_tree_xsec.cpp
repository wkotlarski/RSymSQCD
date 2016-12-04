/*
 *    checked with MadGraph
 */

double Process::sigmaMSSMTree_uu_suLsuR( double s ) {
   double MGl2 = pow(MassGlu, 2);
   double a = pdf->alphasQ( mu_r );
   return (-4.*pow(a, 2)*pi*(2.*sqrt(s*(-4*pow(m1,2) + s)) + (2.*pow(m1,2) - 2.*MGl2 - s)*
      log((4.*MGl2 + pow(1. + sqrt(1. - (4.*pow(m1,2))/s),2.)*s)/
      (4.*MGl2 + pow(-1. + sqrt(1. - (4.*pow(m1,2))/s),2)*s))))/(9.*pow(s,2));
}
