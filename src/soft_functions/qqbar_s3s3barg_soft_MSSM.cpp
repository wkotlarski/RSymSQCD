#include "mathematica_wrapper.hpp"

std::complex<double> qqbar_s3s3barg_soft_MSSM (double Alfas, double s12, double b, double MGl, double muR, double t, double dS) {
   
   double MGl2 = MGl * MGl;
   double b2 = b*b;
   double Alfas2 = Alfas * Alfas;

   return 0.166666666666666667*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*pow(s12,-2.)*sin(t) - 
   0.166666666666666667*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*pow(s12,-2.)*sin(t) + 
   0.333333333333333333*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*pow(s12,-2.)*sin(t) + 
   0.111111111111111111*Alfas*Alfas2*b*(-1. + b2)*pow(s12,-1.)*sin(t) - 1.11111111111111111*Alfas*Alfas2*b*(-1. + b*cos(t))*pow(s12,-1.)*sin(t) + 
   1.55555555555555556*Alfas*Alfas2*b*(1. + b*cos(t))*pow(s12,-1.)*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*(1. + b*cos(t))*pow(s12,-1.)*sin(t) - 
   0.0740740740740740741*Alfas*Alfas2*b*(-1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*sin(t) + 
   0.148148148148148148*Alfas*Alfas2*b*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*sin(t) - 
   0.0185185185185185185*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*(1. + b*cos(t))*
    (1. + 2.*log(dS) + log(s12) + log(pow(s12,-1.)) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*(1. - 1.*b2)*(-2.*log(dS) + log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*sin(t) + 
   0.0740740740740740741*Alfas*Alfas2*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*sin(t) - 
   0.148148148148148148*Alfas*Alfas2*(-1. + b*cos(t))*(1. + b*cos(t))*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b2,-1.)*pow(s12,-1.)*sin(t) + 0.0115740740740740741*Alfas*Alfas2*b*(1. + b*cos(t))*(-2.*log(dS) + log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*
    pow(1. - 1.*b*cos(t),2.)*sin(t) + 0.00694444444444444444*Alfas*Alfas2*b*(-2.*log(dS) + log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*pow(1. - 1.*b*cos(t),3.)*
    sin(t) - 0.0185185185185185185*Alfas*Alfas2*b*(-1. + b2)*pow(s12,-1.)*pow(-1. + b*cos(t),2.)*sin(t) - 
   0.0462962962962962963*Alfas*Alfas2*b*(1. + b*cos(t))*pow(s12,-1.)*pow(-1. + b*cos(t),2.)*sin(t) - 
   0.00231481481481481481*Alfas*Alfas2*b*(-1. + b2)*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*pow(-1. + b*cos(t),2.)*sin(t) + 
   0.0462962962962962963*Alfas*Alfas2*b*pow(s12,-1.)*pow(-1. + b*cos(t),3.)*sin(t) - 
   0.0648148148148148148*Alfas*Alfas2*b*(-1. + b*cos(t))*pow(s12,-1.)*pow(1. + b*cos(t),2.)*sin(t) - 
   0.00231481481481481481*Alfas*Alfas2*b*(-1. + b2)*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*pow(1. + b*cos(t),2.)*sin(t) - 
   0.00694444444444444444*Alfas*Alfas2*b*(-1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*pow(1. + b*cos(t),2.)*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*(-1. + b2)*(1. + 2.*log(dS) + log(s12) + log(pow(s12,-1.)) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*
    pow(1. + b*cos(t),2.)*sin(t) + 0.0648148148148148148*Alfas*Alfas2*b*pow(s12,-1.)*pow(1. + b*cos(t),3.)*sin(t) - 
   0.0115740740740740741*Alfas*Alfas2*b*(-2.*log(dS) + log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,-1.)*pow(1. + b*cos(t),3.)*sin(t) + 
   170.666666666666667*Alfas*Alfas2*b*(-1. + b2)*s12*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   156.444444444444444*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   170.666666666666667*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   85.3333333333333333*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   270.222222222222222*Alfas*Alfas2*b*MGl2*s12*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   7.11111111111111111*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   85.3333333333333333*Alfas*Alfas2*b*(-1. + b2)*s12*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   132.740740740740741*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   244.148148148148148*Alfas*Alfas2*b*MGl2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   7.11111111111111111*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   82.962962962962963*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    sin(t) - 165.925925925925926*Alfas*Alfas2*b*MGl2*s12*(-1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 91.8518518518518519*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   107.851851851851852*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 215.703703703703704*Alfas*Alfas2*b*MGl2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   9.48148148148148148*Alfas*Alfas2*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   9.48148148148148148*Alfas*Alfas2*MGl2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   29.037037037037037*Alfas*Alfas2*b*MGl2*s12*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   4.74074074074074074*Alfas*Alfas2*b*MGl2*s12*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   10.5185185185185185*Alfas*Alfas2*b*MGl2*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 4.74074074074074074*Alfas*Alfas2*MGl2*s12*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 0.592592592592592593*Alfas*Alfas2*b*MGl2*s12*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 227.555555555555556*Alfas*Alfas2*b*s12*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 106.666666666666667*Alfas*Alfas2*b*(-1. + b2)*s12*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 4.74074074074074074*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 265.481481481481481*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 9.48148148148148148*Alfas*Alfas2*b*(-1. + b2)*s12*(1. + b*cos(t))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 4.74074074074074074*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 165.925925925925926*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   25.4814814814814815*Alfas*Alfas2*b*(-1. + b2)*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 150.518518518518519*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   215.703703703703704*Alfas*Alfas2*b*s12*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 9.48148148148148148*Alfas*Alfas2*(-1. + b2)*s12*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    sin(t) + 18.962962962962963*Alfas*Alfas2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 2.37037037037037037*Alfas*Alfas2*b*s12*pow(-1. + b2,2.)*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 129.580246913580247*Alfas*Alfas2*b*pow(MGl2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 3.16049382716049383*Alfas*Alfas2*b*(-1. + b2)*pow(MGl2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 11.8518518518518519*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 56.8888888888888889*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 66.3703703703703704*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*(1. + b*cos(t))*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 41.4814814814814815*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   53.9259259259259259*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 1.18518518518518519*Alfas*Alfas2*(-1. + b2)*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    sin(t) - 14.2222222222222222*Alfas*Alfas2*b*pow(-1. + b2,2.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   2.07407407407407407*Alfas*Alfas2*b*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   16.5925925925925926*Alfas*Alfas2*b*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   2.07407407407407407*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   10.3703703703703704*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 13.5555555555555556*Alfas*Alfas2*b*(-1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   13.4814814814814815*Alfas*Alfas2*b*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 2.37037037037037037*Alfas*Alfas2*(-1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 2.61728395061728395*Alfas*Alfas2*b*pow(-1. + b2,3.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 0.592592592592592593*Alfas*Alfas2*b*(1. + b*cos(t))*pow(-1. + b2,3.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 1.22222222222222222*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b2,3.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   0.592592592592592593*Alfas*Alfas2*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,3.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 0.049382716049382716*Alfas*Alfas2*b*pow(-1. + b2,4.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 149.333333333333333*Alfas*Alfas2*b*s12*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 123.259259259259259*Alfas*Alfas2*b*MGl2*s12*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 2.37037037037037037*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 10.6666666666666667*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 4.74074074074074074*Alfas*Alfas2*b*MGl2*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 141.62962962962963*Alfas*Alfas2*b*MGl2*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   18.962962962962963*Alfas*Alfas2*MGl2*s12*(1. + b*cos(t))*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b2,-1.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   56.8888888888888889*Alfas*Alfas2*b*pow(s12,3.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   23.1111111111111111*Alfas*Alfas2*b*(-1. + b2)*pow(s12,3.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   66.3703703703703704*Alfas*Alfas2*b*(1. + b*cos(t))*pow(s12,3.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   1.18518518518518519*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*pow(s12,3.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   41.4814814814814815*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,3.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 39.5555555555555556*Alfas*Alfas2*b*(-1. + b2)*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(s12,3.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   53.9259259259259259*Alfas*Alfas2*b*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,3.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 2.37037037037037037*Alfas*Alfas2*(-1. + b2)*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,3.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 4.74074074074074074*Alfas*Alfas2*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,3.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 0.592592592592592593*Alfas*Alfas2*b*pow(-1. + b2,2.)*pow(s12,3.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 32.*Alfas*Alfas2*b*s12*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    sin(t) - 33.5802469135802469*Alfas*Alfas2*b*pow(s12,3.)*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) - 
   0.790123456790123457*Alfas*Alfas2*b*(-1. + b2)*pow(s12,3.)*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   1.18518518518518519*Alfas*Alfas2*b*(1. + b*cos(t))*pow(s12,3.)*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 
   34.6666666666666667*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(s12,3.)*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 4.74074074074074074*Alfas*Alfas2*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*pow(s12,3.)*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 42.6666666666666667*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*sin(t) + 2.61728395061728395*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 
   3.0617283950617284*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 6.71604938271604938*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 
   24.7901234567901235*Alfas*Alfas2*b*MGl2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 1.17283950617283951*Alfas*Alfas2*b*MGl2*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 
   0.185185185185185185*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 0.512345679012345679*Alfas*Alfas2*b*s12*(1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 
   0.0956790123456790123*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    sin(t) - 4.79012345679012346*Alfas*Alfas2*b*(-1. + b2)*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 7.50617283950617284*Alfas*Alfas2*b*(-1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 
   4.04938271604938272*Alfas*Alfas2*b*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 6.51851851851851852*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,3.)*
    pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 
   0.407407407407407407*Alfas*Alfas2*b*(-1. + b2)*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 0.543209876543209877*Alfas*Alfas2*b*MGl2*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 
   11.382716049382716*Alfas*Alfas2*b*s12*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 0.790123456790123457*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 
   5.7037037037037037*Alfas*Alfas2*b*(-1. + b2)*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) + 22.8148148148148148*Alfas*Alfas2*b*MGl2*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 
   11.4074074074074074*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*sin(t) - 0.444444444444444444*Alfas*Alfas2*b*(-1. + b2)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    sin(t) - 10.8641975308641975*Alfas*Alfas2*b*(-1. + b2)*MGl2*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.444444444444444444*Alfas*Alfas2*b*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   13.7283950617283951*Alfas*Alfas2*b*MGl2*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   3.67901234567901235*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   9.58024691358024691*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   2.39506172839506173*Alfas*Alfas2*b*(-1. + b2)*s12*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.222222222222222222*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   6.4691358024691358*Alfas*Alfas2*b*(-1. + b2)*MGl2*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   2.32098765432098765*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.148148148148148148*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.432098765432098765*Alfas*Alfas2*b*(-1. + b2)*MGl2*(-1. + b*cos(t))*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.259259259259259259*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.543209876543209877*Alfas*Alfas2*b*(-1. + b2)*MGl2*(-1. + b*cos(t))*(1. + b*cos(t))*
    (1. + 2.*log(dS) + log(s12) + log(pow(s12,-1.)) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.197530864197530864*Alfas*Alfas2*(-1. + b2)*MGl2*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 1.18518518518518519*Alfas*Alfas2*(-1. + b2)*MGl2*(-1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.0987654320987654321*Alfas*Alfas2*(-1. + b2)*s12*(-1. + b*cos(t))*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 1.58024691358024691*Alfas*Alfas2*MGl2*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.592592592592592593*Alfas*Alfas2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.24691358024691358*Alfas*Alfas2*(-1. + b2)*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   1.41975308641975309*Alfas*Alfas2*b*s12*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.148148148148148148*Alfas*Alfas2*b*MGl2*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.444444444444444444*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.932098765432098765*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.0570987654320987654*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.296296296296296296*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*
    (1. + 2.*log(dS) + log(s12) + log(pow(s12,-1.)) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    sin(t) - 0.592592592592592593*Alfas*Alfas2*MGl2*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.049382716049382716*Alfas*Alfas2*s12*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.296296296296296296*Alfas*Alfas2*s12*(-1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.135802469135802469*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*
    (s12 + 2.*s12*log(dS) + s12*log(s12) + s12*log(pow(s12,-1.)) - 1.*s12*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.0185185185185185185*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.0555555555555555556*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.0740740740740740741*Alfas*Alfas2*s12*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.037037037037037037*Alfas*Alfas2*b*(1. + b*cos(t))*
    (s12 + 2.*s12*log(dS) + s12*log(s12) + s12*log(pow(s12,-1.)) - 1.*s12*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 20.7407407407407407*Alfas*Alfas2*b*pow(MGl2,2.)*pow(s12,-1.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.296296296296296296*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*pow(MGl2,2.)*pow(s12,-1.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.888888888888888889*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*pow(MGl2,2.)*pow(s12,-1.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 10.962962962962963*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*pow(MGl2,2.)*pow(s12,-1.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.814814814814814815*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.592592592592592593*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*(1. + 2.*log(dS) + log(s12) + log(pow(s12,-1.)) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(MGl2,2.)*pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   1.18518518518518519*Alfas*Alfas2*(-1. + b2)*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*
    pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   2.37037037037037037*Alfas*Alfas2*(-1. + b*cos(t))*(1. + b*cos(t))*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b2,-1.)*pow(MGl2,2.)*pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.222222222222222222*Alfas*Alfas2*b*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   1.45679012345679012*Alfas*Alfas2*b*(-1. + b2)*MGl2*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   3.45679012345679012*Alfas*Alfas2*b*s12*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   11.2345679012345679*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.444444444444444444*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   3.08024691358024691*Alfas*Alfas2*b*(-1. + b2)*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.00617283950617283951*Alfas*Alfas2*b*(-1. + b2)*MGl2*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.12962962962962963*Alfas*Alfas2*b*(-1. + b2)*s12*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.882716049382716049*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.777777777777777778*Alfas*Alfas2*b*s12*(1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.233024691358024691*Alfas*Alfas2*b*(-1. + b2)*s12*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.296296296296296296*Alfas*Alfas2*(-1. + b2)*s12*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.691358024691358025*Alfas*Alfas2*s12*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.12345679012345679*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*
    (s12 + 2.*s12*log(dS) + s12*log(s12) + s12*log(pow(s12,-1.)) - 1.*s12*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 1.97530864197530864*Alfas*Alfas2*MGl2*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.592592592592592593*Alfas*Alfas2*s12*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.225308641975308642*Alfas*Alfas2*b*s12*pow(-1. + b2,2.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.000771604938271604938*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b2,2.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   2.22222222222222222*Alfas*Alfas2*b*pow(MGl2,2.)*pow(s12,-1.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.037037037037037037*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*pow(s12,-1.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 2.49382716049382716*Alfas*Alfas2*b*MGl2*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.444444444444444444*Alfas*Alfas2*b*s12*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.722222222222222222*Alfas*Alfas2*b*(-1. + b2)*s12*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 2.87654320987654321*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.00617283950617283951*Alfas*Alfas2*b*MGl2*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.518518518518518519*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.00771604938271604938*Alfas*Alfas2*b*(-1. + b2)*s12*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.237654320987654321*Alfas*Alfas2*b*s12*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.395061728395061728*Alfas*Alfas2*s12*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.691358024691358025*Alfas*Alfas2*b*s12*pow(-1. + b*cos(t),4.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.012345679012345679*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b*cos(t),4.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   4.71604938271604938*Alfas*Alfas2*b*(-1. + b2)*MGl2*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   5.92592592592592593*Alfas*Alfas2*b*s12*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.148148148148148148*Alfas*Alfas2*b*(-1. + b2)*s12*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   8.96296296296296296*Alfas*Alfas2*b*MGl2*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.592592592592592593*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   2.2654320987654321*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.586419753086419753*Alfas*Alfas2*b*(-1. + b2)*MGl2*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.12962962962962963*Alfas*Alfas2*b*(-1. + b2)*s12*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   1.10493827160493827*Alfas*Alfas2*b*MGl2*(-1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.518518518518518519*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.282407407407407407*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 0.049382716049382716*Alfas*Alfas2*b*(-1. + b2)*MGl2*
    (1. + 2.*log(dS) + log(s12) + log(pow(s12,-1.)) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    sin(t) - 0.0987654320987654321*Alfas*Alfas2*s12*(-1. + b*cos(t))*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.024691358024691358*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*
    (s12 + 2.*s12*log(dS) + s12*log(s12) + s12*log(pow(s12,-1.)) - 1.*s12*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.395061728395061728*Alfas*Alfas2*MGl2*(-1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.592592592592592593*Alfas*Alfas2*s12*(-1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.595679012345679012*Alfas*Alfas2*b*s12*pow(-1. + b2,2.)*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.074845679012345679*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b2,2.)*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.012345679012345679*Alfas*Alfas2*b*(s12 + 2.*s12*log(dS) + s12*log(s12) + s12*log(pow(s12,-1.)) - 1.*s12*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,2.)*
    pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   9.33333333333333333*Alfas*Alfas2*b*pow(MGl2,2.)*pow(s12,-1.)*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   1.14814814814814815*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(MGl2,2.)*pow(s12,-1.)*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 2.14814814814814815*Alfas*Alfas2*b*s12*pow(-1. + b*cos(t),2.)*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.265432098765432099*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b*cos(t),2.)*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.197530864197530864*Alfas*Alfas2*s12*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*
    pow(-1. + b*cos(t),2.)*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.37037037037037037*Alfas*Alfas2*b*MGl2*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.592592592592592593*Alfas*Alfas2*b*s12*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.0925925925925925926*Alfas*Alfas2*b*(-1. + b2)*s12*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.185185185185185185*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.0679012345679012346*Alfas*Alfas2*b*MGl2*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.777777777777777778*Alfas*Alfas2*b*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) + 
   0.0169753086419753086*Alfas*Alfas2*b*(-1. + b2)*s12*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 0.0339506172839506173*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*sin(t) - 
   0.111111111111111111*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) + 
   0.111111111111111111*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) + 
   1.0617283950617284*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) - 
   0.00308641975308641975*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) + 0.00308641975308641975*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*
    (2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) + 
   0.0925925925925925926*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) - 0.024691358024691358*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*
    (1. + 2.*log(dS) + log(s12) + log(pow(s12,-1.)) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) - 
   0.049382716049382716*Alfas*Alfas2*(-1. + b2)*(2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) + 0.197530864197530864*Alfas*Alfas2*(-1. + b*cos(t))*(1. + b*cos(t))*
    (2.*b*log(dS) - 1.*log((1. + b)*pow(1. - 1.*b,-1.)) - 1.*b*log(pow(muR,2.)*pow(s12,-1.)))*pow(-1. + b2,-1.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) - 0.234567901234567901*Alfas*Alfas2*b*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) - 0.00308641975308641975*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*
    pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) - 
   0.87654320987654321*Alfas*Alfas2*b*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*sin(t) - 
   0.114197530864197531*Alfas*Alfas2*b*(2.*log(dS) - 1.*log(pow(muR,2.)*pow(s12,-1.)))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*
    sin(t) + 0.00115740740740740741*Alfas*Alfas2*(-1. + b2)*pow(s12,-1.)*pow(-1. + b*cos(t),2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.00115740740740740741*Alfas*Alfas2*(1. + b*cos(t))*pow(s12,-1.)*pow(-1. + b*cos(t),2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00115740740740740741*Alfas*Alfas2*pow(s12,-1.)*pow(-1. + b*cos(t),3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.00115740740740740741*Alfas*Alfas2*(-1. + b2)*pow(s12,-1.)*pow(1. + b*cos(t),2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00115740740740740741*Alfas*Alfas2*(-1. + b*cos(t))*pow(s12,-1.)*pow(1. + b*cos(t),2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.00115740740740740741*Alfas*Alfas2*pow(s12,-1.)*pow(1. + b*cos(t),3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 9.77777777777777778*Alfas*Alfas2*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 2.07407407407407407*Alfas*Alfas2*(-1. + b2)*MGl2*s12*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 4.14814814814814815*Alfas*Alfas2*MGl2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 4.*Alfas*Alfas2*MGl2*s12*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 1.77777777777777778*Alfas*Alfas2*MGl2*s12*(-1. + b*cos(t))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.962962962962962963*Alfas*Alfas2*MGl2*s12*pow(-1. + b2,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 7.7037037037037037*Alfas*Alfas2*(-1. + b2)*s12*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 2.96296296296296296*Alfas*Alfas2*s12*(-1. + b*cos(t))*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.592592592592592593*Alfas*Alfas2*(-1. + b2)*s12*(-1. + b*cos(t))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 4.14814814814814815*Alfas*Alfas2*s12*(1. + b*cos(t))*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 1.77777777777777778*Alfas*Alfas2*s12*pow(-1. + b2,2.)*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.790123456790123457*Alfas*Alfas2*pow(MGl2,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.395061728395061728*Alfas*Alfas2*(-1. + b2)*pow(MGl2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 1.03703703703703704*Alfas*Alfas2*(-1. + b2)*(-1. + b*cos(t))*(1. + b*cos(t))*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 2.25925925925925926*Alfas*Alfas2*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.259259259259259259*Alfas*Alfas2*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.50617283950617284*Alfas*Alfas2*pow(-1. + b2,3.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.481481481481481481*Alfas*Alfas2*(-1. + b*cos(t))*pow(-1. + b2,3.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.12345679012345679*Alfas*Alfas2*pow(-1. + b2,4.)*pow(s12,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 3.55555555555555556*Alfas*Alfas2*MGl2*s12*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.296296296296296296*Alfas*Alfas2*(-1. + b2)*MGl2*s12*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 2.96296296296296296*Alfas*Alfas2*(-1. + b2)*pow(s12,3.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 1.03703703703703704*Alfas*Alfas2*(1. + b*cos(t))*pow(s12,3.)*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.444444444444444444*Alfas*Alfas2*pow(-1. + b2,2.)*pow(s12,3.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.938271604938271605*Alfas*Alfas2*pow(s12,3.)*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.049382716049382716*Alfas*Alfas2*(-1. + b2)*pow(s12,3.)*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.271604938271604938*Alfas*Alfas2*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.351851851851851852*Alfas*Alfas2*MGl2*s12*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.0802469135802469136*Alfas*Alfas2*MGl2*s12*(1. + b*cos(t))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.0401234567901234568*Alfas*Alfas2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.0802469135802469136*Alfas*Alfas2*MGl2*s12*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.0478395061728395062*Alfas*Alfas2*s12*(-1. + b*cos(t))*pow(-1. + b2,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00771604938271604938*Alfas*Alfas2*s12*(1. + b*cos(t))*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.00771604938271604938*Alfas*Alfas2*s12*pow(-1. + b2,4.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.543209876543209877*Alfas*Alfas2*(-1. + b2)*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.888888888888888889*Alfas*Alfas2*(-1. + b*cos(t))*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.444444444444444444*Alfas*Alfas2*(-1. + b2)*(-1. + b*cos(t))*pow(MGl2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.271604938271604938*Alfas*Alfas2*pow(-1. + b2,2.)*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.592592592592592593*Alfas*Alfas2*pow(MGl2,3.)*pow(s12,-1.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.296296296296296296*Alfas*Alfas2*(-1. + b2)*pow(MGl2,3.)*pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.75308641975308642*Alfas*Alfas2*(-1. + b2)*MGl2*s12*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.132716049382716049*Alfas*Alfas2*(-1. + b2)*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.481481481481481481*Alfas*Alfas2*MGl2*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.172839506172839506*Alfas*Alfas2*s12*pow(-1. + b2,2.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.299382716049382716*Alfas*Alfas2*(-1. + b2)*s12*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.481481481481481481*Alfas*Alfas2*MGl2*s12*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.166666666666666667*Alfas*Alfas2*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.166666666666666667*Alfas*Alfas2*s12*pow(-1. + b*cos(t),4.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.259259259259259259*Alfas*Alfas2*(-1. + b2)*MGl2*s12*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.194444444444444444*Alfas*Alfas2*(-1. + b2)*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.259259259259259259*Alfas*Alfas2*MGl2*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.0648148148148148148*Alfas*Alfas2*s12*pow(-1. + b2,2.)*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.12962962962962963*Alfas*Alfas2*s12*pow(-1. + b*cos(t),2.)*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.0648148148148148148*Alfas*Alfas2*(-1. + b2)*s12*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.259259259259259259*Alfas*Alfas2*MGl2*s12*pow(1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.12962962962962963*Alfas*Alfas2*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.037037037037037037*Alfas*Alfas2*(-1. + b2)*MGl2*(-1. + b*cos(t))*(1. + b*cos(t))*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.0555555555555555556*Alfas*Alfas2*MGl2*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.0185185185185185185*Alfas*Alfas2*MGl2*(1. + b*cos(t))*pow(-1. + b2,2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.00925925925925925926*Alfas*Alfas2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.0185185185185185185*Alfas*Alfas2*MGl2*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.0115740740740740741*Alfas*Alfas2*s12*(-1. + b*cos(t))*pow(-1. + b2,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00231481481481481481*Alfas*Alfas2*s12*(1. + b*cos(t))*pow(-1. + b2,3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.00231481481481481481*Alfas*Alfas2*s12*pow(-1. + b2,4.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.0740740740740740741*Alfas*Alfas2*(-1. + b2)*pow(MGl2,2.)*pow(s12,-1.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.037037037037037037*Alfas*Alfas2*pow(-1. + b2,2.)*pow(MGl2,2.)*pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.0308641975308641975*Alfas*Alfas2*(-1. + b2)*MGl2*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.00617283950617283951*Alfas*Alfas2*MGl2*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00771604938271604938*Alfas*Alfas2*(-1. + b2)*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.0169753086419753086*Alfas*Alfas2*s12*pow(-1. + b2,2.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00617283950617283951*Alfas*Alfas2*MGl2*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.00462962962962962963*Alfas*Alfas2*(-1. + b2)*s12*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00308641975308641975*Alfas*Alfas2*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.00308641975308641975*Alfas*Alfas2*s12*pow(-1. + b*cos(t),4.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.00617283950617283951*Alfas*Alfas2*(-1. + b2)*MGl2*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.00617283950617283951*Alfas*Alfas2*MGl2*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00462962962962962963*Alfas*Alfas2*(-1. + b2)*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.00154320987654320988*Alfas*Alfas2*s12*pow(-1. + b2,2.)*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.00308641975308641975*Alfas*Alfas2*s12*pow(-1. + b*cos(t),2.)*pow(1. + b*cos(t),2.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) + 0.00617283950617283951*Alfas*Alfas2*MGl2*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) + 0.00154320987654320988*Alfas*Alfas2*(-1. + b2)*s12*pow(1. + b*cos(t),3.)*
    pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*(4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 
      2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*
    sin(t) - 0.00308641975308641975*Alfas*Alfas2*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (4.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1.)) - 2.*log((1. + b)*pow(1. - 1.*b,-1.))*log(pow(muR,2.)*pow(s12,-1.)) - 4.*PolyLog(2.,2.*b*pow(1. + b,-1.)) - 
      1.*pow(log((1. + b)*pow(1. - 1.*b,-1.)),2.))*sin(t) - 0.00115740740740740741*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*(1. + b*cos(t))*pow(s12,-1.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*(1. + b*cos(t))*pow(s12,-1.)*pow(-1. + b*cos(t),2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(-1. + b2)*pow(s12,-1.)*pow(1. + b*cos(t),2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.00231481481481481481*Alfas*Alfas2*b*(-1. + b*cos(t))*pow(s12,-1.)*pow(1. + b*cos(t),2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   4.14814814814814815*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   2.07407407407407407*Alfas*Alfas2*b*MGl2*s12*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   4.14814814814814815*Alfas*Alfas2*b*(-1. + b2)*s12*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   1.03703703703703704*Alfas*Alfas2*b*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.259259259259259259*Alfas*Alfas2*b*pow(-1. + b2,3.)*pow(s12,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   1.03703703703703704*Alfas*Alfas2*b*(-1. + b2)*pow(s12,3.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-4.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.617283950617283951*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.777777777777777778*Alfas*Alfas2*b*(-1. + b2)*MGl2*s12*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.907407407407407407*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   2.07407407407407407*Alfas*Alfas2*b*MGl2*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.469135802469135802*Alfas*Alfas2*b*MGl2*s12*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.104938271604938272*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.194444444444444444*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.12345679012345679*Alfas*Alfas2*b*s12*pow(-1. + b2,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.0987654320987654321*Alfas*Alfas2*b*(-1. + b2)*pow(MGl2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.283950617283950617*Alfas*Alfas2*b*(-1. + b2)*s12*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   1.03703703703703704*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-3.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.152777777777777778*Alfas*Alfas2*b*(-1. + b2)*MGl2*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.0972222222222222222*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.0216049382716049383*Alfas*Alfas2*b*MGl2*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.0648148148148148148*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.0401234567901234568*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.0648148148148148148*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.0304783950617283951*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*(1. + b*cos(t))*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.00308641975308641975*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(-1. + b2,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.00540123456790123457*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b2,3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.037037037037037037*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*pow(MGl2,2.)*pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.0740740740740740741*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*pow(MGl2,2.)*pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.148148148148148148*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*pow(MGl2,2.)*pow(s12,-1.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.0416666666666666667*Alfas*Alfas2*b*(-1. + b2)*MGl2*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.0648148148148148148*Alfas*Alfas2*b*(-1. + b2)*s12*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.135802469135802469*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.12962962962962963*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.0547839506172839506*Alfas*Alfas2*b*(-1. + b2)*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.0119598765432098765*Alfas*Alfas2*b*s12*pow(-1. + b2,2.)*pow(-1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.0115740740740740741*Alfas*Alfas2*b*(-1. + b2)*s12*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.0308641975308641975*Alfas*Alfas2*b*s12*(1. + b*cos(t))*pow(-1. + b*cos(t),3.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.0324074074074074074*Alfas*Alfas2*b*(-1. + b2)*s12*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.012345679012345679*Alfas*Alfas2*b*MGl2*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.12962962962962963*Alfas*Alfas2*b*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.00308641975308641975*Alfas*Alfas2*b*(-1. + b2)*s12*(-1. + b*cos(t))*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.00617283950617283951*Alfas*Alfas2*b*s12*pow(-1. + b*cos(t),2.)*pow(1. + b*cos(t),2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-2.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.00385802469135802469*Alfas*Alfas2*b*(-1. + b2)*(-1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) - 
   0.00617283950617283951*Alfas*Alfas2*b*(-1. + b2)*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.012345679012345679*Alfas*Alfas2*b*(-1. + b*cos(t))*(1. + b*cos(t))*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.00154320987654320988*Alfas*Alfas2*b*pow(-1. + b2,2.)*pow(4.*MGl2 + s12 + b2*s12 - 2.*b*s12*cos(t),-1.)*
    (-4.*log(dS)*log(pow(muR,2.)*pow(s12,-1.)) + 4.*pow(log(dS),2.) + pow(log(pow(muR,2.)*pow(s12,-1.)),2.))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000434027777777777778*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,3.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*pow(1. + b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00520833333333333333*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00130208333333333333*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*pow(1. + b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. + b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000868055555555555556*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),3.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0324074074074074074*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00289351851851851852*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00202546296296296296*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. + b*cos(t))*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0115740740740740741*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*(1. + b*cos(t))*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*(1. + b*cos(t))*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00173611111111111111*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.072964891975308642*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0180844907407407407*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0362172067901234568*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00202546296296296296*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00607638888888888889*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,3.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0485628858024691358*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0496720679012345679*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0409915123456790123*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0567611882716049383*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. + b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.057195216049382716*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. + b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0140817901234567901*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00144675925925925926*Alfas*Alfas2*b*pow(1. + b*cos(t),3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00072337962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),3.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0115740740740740741*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.348572530864197531*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. - 1.*b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00520833333333333333*Alfas*Alfas2*b*(1. - 1.*b2)*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.18132716049382716*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0460069444444444444*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.176697530864197531*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00217013888888888889*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.214409722222222222*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0653452932098765432*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0234375*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0250289351851851852*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00578703703703703704*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0450424382716049383*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00072337962962962963*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(1. - 1.*b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0425347222222222222*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,3.)*pow(1. - 1.*b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00385802469135802469*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,3.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00354456018518518519*Alfas*Alfas2*b*pow(1. - 1.*b2,4.)*pow(1. - 1.*b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00405092592592592593*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.073591820987654321*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.176311728395061728*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0291280864197530864*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0293209876543209877*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),3.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. + b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00810185185185185185*Alfas*Alfas2*b*MGl2*pow(1. + b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00636574074074074074*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*pow(1. + b*cos(t),3.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),3.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0648148148148148148*Alfas*Alfas2*b*(log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0162037037037037037*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0324074074074074074*Alfas*Alfas2*b*(1. - 1.*b2)*(1. + b*cos(t))*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0405092592592592593*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*(1. + b*cos(t))*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0324074074074074074*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),-1.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00810185185185185185*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0243055555555555556*Alfas*Alfas2*b*pow(1. + b*cos(t),2.)*(log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0162037037037037037*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.056712962962962963*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0324074074074074074*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.162037037037037037*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.194444444444444444*Alfas*Alfas2*b*(1. + b*cos(t))*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0972222222222222222*Alfas*Alfas2*b*(1. - 1.*b2)*(1. + b*cos(t))*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.226851851851851852*Alfas*Alfas2*b*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.453703703703703704*Alfas*Alfas2*b*pow(MGl2,3.)*pow(1. - 1.*b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.226851851851851852*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,3.)*pow(1. - 1.*b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0104166666666666667*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*(1. + b*cos(t))*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000289351851851851852*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00694444444444444444*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.000241126543209876543*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00144675925925925926*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(1. - 1.*b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000144675925925925926*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,3.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0000361689814814814815*Alfas*Alfas2*b*pow(1. - 1.*b2,4.)*pow(1. - 1.*b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000192901234567901235*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*(1. + b*cos(t))*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0138888888888888889*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. - 1.*b*cos(t))*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000771604938271604938*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00260416666666666667*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0277777777777777778*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.049382716049382716*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.015625*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00308641975308641975*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0138888888888888889*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00617283950617283951*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),3.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00173611111111111111*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,3.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0138888888888888889*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),3.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00347222222222222222*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),4.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00173611111111111111*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. + b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 0.5*s12*pow(log(s12),2.) + 
      s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00144675925925925926*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00578703703703703704*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0208333333333333333*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00347222222222222222*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(1. + b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0182291666666666667*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00434027777777777778*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,3.)*pow(1. + b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0451388888888888889*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0173611111111111111*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(MGl2,2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.025462962962962963*Alfas*Alfas2*b*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00810185185185185185*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),3.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00607638888888888889*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00800540123456790123*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00144675925925925926*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0138888888888888889*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00453317901234567901*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0152391975308641975*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0113811728395061728*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000868055555555555556*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,3.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00144675925925925926*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),3.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),4.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00530478395061728395*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. + b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0106095679012345679*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. + b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0447530864197530864*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. - 1.*b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0258487654320987654*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0117669753086419753*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00964506172839506173*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00202546296296296296*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0142746913580246914*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0270061728395061728*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0111882716049382716*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0177469135802469136*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. - 1.*b*cos(t))*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00771604938271604938*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0044367283950617284*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.000964506172839506173*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000964506172839506173*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000192901234567901235*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,3.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000144675925925925926*Alfas*Alfas2*b*pow(1. - 1.*b2,4.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0050154320987654321*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0158179012345679012*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.0033757716049382716*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00655864197530864198*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),3.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00694444444444444444*Alfas*Alfas2*b*MGl2*pow(1. + b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00347222222222222222*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. + b*cos(t),2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (s12*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 1.*(-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-1.*s12*log(s12) + s12*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + s12*log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 0.5*s12*pow(log(s12),2.) - 
      1.*s12*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 0.5*s12*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      s12*(2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) + 
      0.5*pow(s12,-1.)*pow(log(s12),2.) + pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) + 
      0.5*pow(s12,-1.)*pow(log(s12),2.) + pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0740740740740740741*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) + 
      0.5*pow(s12,-1.)*pow(log(s12),2.) + pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0277777777777777778*Alfas*Alfas2*b*(1. + b*cos(t))*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) + 
      0.5*pow(s12,-1.)*pow(log(s12),2.) + pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) + 
      0.5*pow(s12,-1.)*pow(log(s12),2.) + pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,3.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) + 
      0.5*pow(s12,-1.)*pow(log(s12),2.) + pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.025462962962962963*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*(1. + b*cos(t))*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0347222222222222222*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),3.)*pow(1. + b*cos(t),-1.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0115740740740740741*Alfas*Alfas2*b*pow(1. + b*cos(t),2.)*(log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0138888888888888889*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.0185185185185185185*Alfas*Alfas2*b*pow(MGl2,3.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,-1.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,-1.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,-1.)) - 
      0.5*pow(s12,-1.)*pow(log(s12),2.) - 1.*pow(s12,-1.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,-1.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,-1.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(1. + b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + 
      log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) - 
      1.*log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) + 
      0.5*pow(s12,3.)*pow(log(s12),2.) + pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) + 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 1.*pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00462962962962962963*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00159143518518518519*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000361689814814814815*Alfas*Alfas2*b*pow(1. - 1.*b2,3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00202546296296296296*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00289351851851851852*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00209780092592592593*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,3.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000361689814814814815*Alfas*Alfas2*b*pow(1. - 1.*b2,4.)*pow(1. + b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00419560185185185185*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),3.)*pow(1. + b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) - 
   0.00318287037037037037*Alfas*Alfas2*b*(1. - 1.*b2)*pow(1. - 1.*b*cos(t),3.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.000578703703703703704*Alfas*Alfas2*b*pow(1. - 1.*b*cos(t),4.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) - 
      1.*log(pow(muR,2.)*pow(s12,-1.))*(-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.)) + 
      log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.))*
       (-1.*log(s12)*pow(s12,3.) + (-2.*log(dS) - 1.*log(pow(s12,-1.)))*pow(s12,3.) + log(pow(muR,2.)*pow(s12,-1.))*pow(s12,3.)) - 
      0.5*pow(s12,3.)*pow(log(s12),2.) - 1.*pow(s12,3.)*(2.*log(dS)*log(pow(s12,-1.)) + 2.*pow(log(dS),2.) + 0.5*pow(log(pow(s12,-1.)),2.)) - 
      0.5*pow(s12,3.)*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + pow(s12,3.)*
       (2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
         1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
         0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.)))*sin(t) + 
   0.00231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,3.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + 2.*log(dS)*log(pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*pow(log(dS),2.) + 0.5*pow(log(s12),2.) + 0.5*pow(log(pow(s12,-1.)),2.) + 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) - 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00154320987654320988*Alfas*Alfas2*b*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-1.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + 2.*log(dS)*log(pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*pow(log(dS),2.) + 0.5*pow(log(s12),2.) + 0.5*pow(log(pow(s12,-1.)),2.) + 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) - 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00115740740740740741*Alfas*Alfas2*b*(1. - 1.*b2)*(1. + b*cos(t))*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-1.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + 2.*log(dS)*log(pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*pow(log(dS),2.) + 0.5*pow(log(s12),2.) + 0.5*pow(log(pow(s12,-1.)),2.) + 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) - 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00694444444444444444*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + 2.*log(dS)*log(pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*pow(log(dS),2.) + 0.5*pow(log(s12),2.) + 0.5*pow(log(pow(s12,-1.)),2.) + 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) - 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00347222222222222222*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),3.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + 2.*log(dS)*log(pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*pow(log(dS),2.) + 0.5*pow(log(s12),2.) + 0.5*pow(log(pow(s12,-1.)),2.) + 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) - 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.024691358024691358*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + 2.*log(dS)*log(pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*pow(log(dS),2.) + 0.5*pow(log(s12),2.) + 0.5*pow(log(pow(s12,-1.)),2.) + 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) - 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00154320987654320988*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + 2.*log(dS)*log(pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*pow(log(dS),2.) + 0.5*pow(log(s12),2.) + 0.5*pow(log(pow(s12,-1.)),2.) + 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) - 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00925925925925925926*Alfas*Alfas2*b*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (-1.*log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) + 2.*log(dS)*log(pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) - 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) + 
      2.*pow(log(dS),2.) + 0.5*pow(log(s12),2.) + 0.5*pow(log(pow(s12,-1.)),2.) + 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) + 
      pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) - 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00154320987654320988*Alfas*Alfas2*b*pow(MGl2,3.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.112847222222222222*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.193479938271604938*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.196566358024691358*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0162037037037037037*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0486111111111111111*Alfas*Alfas2*b*MGl2*(1. + b*cos(t))*pow(1. - 1.*b2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.163001543209876543*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.227816358024691358*Alfas*Alfas2*b*MGl2*pow(1. + b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.113040123456790123*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00289351851851851852*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),3.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.000771604938271604938*Alfas*Alfas2*b*(1. - 1.*b2)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-1.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0115740740740740741*Alfas*Alfas2*b*pow(MGl2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.513888888888888889*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.348765432098765432*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(MGl2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.350308641975308642*Alfas*Alfas2*b*(1. + b*cos(t))*pow(MGl2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.0115740740740740741*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0115740740740740741*Alfas*Alfas2*b*(1. + b*cos(t))*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.175154320987654321*Alfas*Alfas2*b*(1. - 1.*b2)*(1. + b*cos(t))*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.170138888888888889*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.0115740740740740741*Alfas*Alfas2*b*pow(MGl2,2.)*pow(1. - 1.*b*cos(t),-1.)*pow(1. + b*cos(t),2.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. - 1.*b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. - 1.*b*cos(t),-1.)*(-1.*b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(-1.*b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. - 1.*b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0462962962962962963*Alfas*Alfas2*b*pow(MGl2,3.)*pow(1. + b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.0231481481481481481*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,3.)*pow(1. + b*cos(t),-1.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-4.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0505401234567901235*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0100308641975308642*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0165895061728395062*Alfas*Alfas2*b*MGl2*(1. - 1.*b*cos(t))*(1. + b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.0115740740740740741*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b2,2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.0142746913580246914*Alfas*Alfas2*b*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0127314814814814815*Alfas*Alfas2*b*(1. - 1.*b2)*MGl2*pow(1. - 1.*b*cos(t),2.)*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0212191358024691358*Alfas*Alfas2*b*MGl2*pow(1. + b*cos(t),2.)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-2.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.00771604938271604938*Alfas*Alfas2*b*(1. - 1.*b2)*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-1.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0925925925925925926*Alfas*Alfas2*b*(1. - 1.*b*cos(t))*pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-1.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00385802469135802469*Alfas*Alfas2*b*(1. - 1.*b2)*(1. - 1.*b*cos(t))*pow(1. + b*cos(t),-1.)*
    pow(-1.*MGl2 + 0.25*(1. - 1.*b2)*s12 - 0.5*s12*(1. - 1.*b*cos(t)),-1.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0277777777777777778*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) - 
   0.0540123456790123457*Alfas*Alfas2*b*(1. + b*cos(t))*pow(MGl2,2.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.0154320987654320988*Alfas*Alfas2*b*(1. - 1.*b2)*pow(MGl2,2.)*pow(1. + b*cos(t),-1.)*pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t) + 
   0.00694444444444444444*Alfas*Alfas2*b*pow(1. - 1.*b2,2.)*pow(MGl2,2.)*pow(1. + b*cos(t),-1.)*
    pow(MGl2 - 0.25*(1. - 1.*b2)*s12 + 0.5*s12*(1. - 1.*b*cos(t)),-3.)*
    (log(s12)*(-2.*log(dS) - 1.*log(pow(s12,-1.))) - 2.*log(dS)*log(pow(s12,-1.)) - 
      1.*(-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)))*log(pow(muR,2.)*pow(s12,-1.)) + 
      (-2.*log(dS) - 1.*log(s12) - 1.*log(pow(s12,-1.)) + log(pow(muR,2.)*pow(s12,-1.)))*
       log(pow(1. + b*cos(t),2.)*pow(1. - 1.*b2*pow(cos(t),2.) - 1.*b2*pow(sin(t),2.),-1.)) + 
      2.*PolyLog(2.,pow(1. + b*cos(t),-1.)*(b*cos(t) - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*PolyLog(2.,-1.*pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(b*cos(t) + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))) - 
      2.*pow(log(dS),2.) - 0.5*pow(log(s12),2.) - 0.5*pow(log(pow(s12,-1.)),2.) - 0.5*pow(log(pow(muR,2.)*pow(s12,-1.)),2.) - 
      1.*pow(log(pow(1. + b*cos(t),-1.)*(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.) + 
      0.5*pow(log(pow(1. - 1.*pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5),-1.)*(1. + pow(b2*pow(cos(t),2.) + b2*pow(sin(t),2.),0.5))),2.))*sin(t);
}
