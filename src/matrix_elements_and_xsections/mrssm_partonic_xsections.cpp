double Process::sigmaMRSSMTree_uubar_suLsuLdagger( double s12 ) {
	double Alfas2 = pow( pdf->alphasQ( mu_r ), 2);
   return (Alfas2*pi*(8*(-3*(m1*m1) + MassGlu*MassGlu - 2*s12)*Sqrt(s12*(-4*(m1*m1) + s12)) + 
       8*(2*Power(m1,4) + 2*Power(MassGlu,4) - 4*(MassGlu*MassGlu)*s12 - 3*(s12*s12) + m1*m1*(-4*(MassGlu*MassGlu) + 6*s12))*
        atanh(Sqrt(s12*(-4*(m1*m1) + s12))/(2*(m1*m1) - 2*(MassGlu*MassGlu) - s12))))/(54.*Power(s12,3));
} 

double Process::sigmaMRSSMTree_ddbar_suLsuLdagger( double s12 ) {
   double Alfas2 = pow( pdf->alphasQ( mu_r ), 2);
   return (2*Alfas2*pi*Power(-4*(m1*m1) + s12,1.5))/(27.*Power(s12,2.5));
}

double Process::sigmaMRSSMTree_gg_suLsuLdagger( double s12 ) {
   double Alfas2 = pow( pdf->alphasQ( mu_r ), 2);
   return (Alfas2*pi*(Sqrt(s12*(-4*(m1*m1) + s12))*(62*(m1*m1) + 5*s12) - 16*(m1*m1)*(m1*m1 + 4*s12)*atanh(Sqrt(1 - (4*(m1*m1))/s12))))/(48.*Power(s12,3));
}
