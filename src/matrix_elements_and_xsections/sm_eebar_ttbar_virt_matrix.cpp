//
// Created by wojciech on 05.08.17.
//

#include "dilog.hpp"
double Process::matrixSMVirt_eebar_ttbar(double S, double T,
                                         const double FiniteGs, const double Dminus4, int divergence) {
   using dilogarithm::dilog;
   double alphaS = 0.12;
   double Alfa2 = pow (137., -2);
   double MB2 = pow(m1, 2);
   double U = pow(m1, 2) + pow(m2, 2) - S - T;
   double b = sqrt(1 - 4*m1*m1/S);
   double CF = 4/3.;
   const double born = 4.*(8*Alfa2*(pi*pi)*(S*S + 2*(MB2*MB2 - 2*MB2*T + T*(S + T))))/(3.*(S*S));
   // eq. 3.24 of Harris & Owens '02
   const double A1 = -2 * CF * (1 - (1+b*b)/(2*b) * log((1+b)/(1-b)));
   switch (divergence) {
      case -2:
         return 0.;
      case -1:
         return born * alphaS/(2.*pi) * A1;
      case  0:
         return born * alphaS/(2.*pi) * (
            A1*log(91.188*91.188/S)
            + CF*(-2*(1-(1+b*b)/(2*b)*log((1+b)/(1-b)))*log(S/MB2)
               + 3*b*log((1+b)/(1-b)) - 4 + (1+b*b)/b*(-0.5*pow(log((1-b)/(1+b)), 2)
                  + 2*log((1-b)/(1+b))*log(2*b/(1+b)) + 2*dilog((1-b)/(1+b)) + 2/3.*pi*pi)
            )
         )
         // eq. 3.25
         +  8*pi*pi*Alfa2/S * pow(1/3.,2) * alphaS/(2*pi)
            *(4*3*CF*(b*b-1)/b*log((1-b)/(1+b)) * (MB2/S -(T-MB2)*(U-MB2)/(S*S)));
   }
}
