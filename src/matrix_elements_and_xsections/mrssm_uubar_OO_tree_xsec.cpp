inline double Process::matrixSgluonTree_qqbar_OO( double s ) {
   double a = pdf->alphasQ( mu_r );
   return 2. * pi * pow(a, 2)/(9. * s) * pow(1. - 4. * m1*m1/s, 3./2. );
}
