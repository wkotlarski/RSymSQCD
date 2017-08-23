//
// Created by wojciech on 05.08.17.
//

#include "dilog.hpp"
double Process::matrixSMVirt_eebar_ttbar(double S, double T,
                                         const double FiniteGs, const double Dminus4, int divergence) {
   using dilogarithm::dilog;
   ltini(); // for LoopTools
   setmudim(pow(mu_r, 2));
   setlambda(divergence);
   double alphaS = pdf->alphasQ(mu_r);
   double Alfa2 = pow (137., -2);
   double mu = mu_r;
   double MT2 = pow(5, 2);
   double MB2 = pow(5, 2);
   double U = pow(m1, 2) + pow(m2, 2) - S - T;
   double Finite = 0.;
   int Divergence;    // UV divergence from gauge coupling
   if (divergence == 0 || divergence == -2) {
      Divergence = 0;
   } else if (divergence == -1) {
      Divergence = 1;
   }
   double b = sqrt(1 - 4*MT2/S);
   double CF = 4/3.;
   const double born = 4.*(8*Alfa2*(pi*pi)*(S*S + 2*(MB2*MB2 - 2*MB2*T + T*(S + T))))/(3.*(S*S));
   std::cout << "virt: " << born * alphaS/(2.*pi) * (-2*CF * (1 - (1+b*b)/(2*b)*log((1+b)/(1-b)) )) << std::endl;
   double matrix = 2.*3*pow(Alfa2*4*pi,2)*pow(1/3.,2)*( (T+U-2*MT2)/S + 2*MT2/S) *
      4./3.*(-2*(1-(1+b*b)/(2*b)*log((1+b)/(1-b)))*log(S/MT2) + 3*b*log((1+b)/(1-b))-4
      + (1 + b*b)/b*(-0.5*pow(log((1-b)/(1+b)),2) + 2*log((1-b)/(1+b))*log(2*b/(1+b))
           +2*dilog((1-b)/(1+b)) + 2/3. * pi*pi)
   )
   // eq. 3.25
   +  8*pi*pi*Alfa2/S*pow(1/3.,2)*alphaS/(2*pi)*(4*3*4/3.*(b*b-1)/b*log((1-b)/(b+1))*(MT2/S- (T-MT2)*(U-MT2)/(S*S))
   )
   ;
   return matrix;
}
