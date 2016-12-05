inline double Process::matrixSgluonTree_gg_OO( double s ) {
   double b = sqrt( 1. - 4. * m1*m1/s );
   double a = pdf->alphasQ( mu_r );
   return 3. * pi * pow(a, 2)/(32. * s) 
      * (27.*b - 17*pow(b, 3) + 6.*(-3. + 2. * pow(b, 2) + pow(b, 4) ) * atanh(b) ) ;
}
