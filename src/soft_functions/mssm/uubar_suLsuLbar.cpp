#include "mathematica_wrapper.hpp"

std::complex<double> uubar_suLsuLbarg_soft_MSSM(double Alfas, double s12, 
        double b, double MGl2, double muR, double t, double dS) {
   
   double Alfas2 = Alfas * Alfas;
   
return -0.00077160493827160493827*Alfas*Alfas2*pow(b,2)*pow(s12,-2)*
  pow(4.*MGl2 + s12 - 2.*b*s12*cos(t) + s12*pow(b,2),-4)*
  (110592.*b*s12*pow(MGl2,3) - 221184.*s12*cos(t)*pow(b,2)*pow(MGl2,3) + 
    110592.*s12*pow(b,3)*pow(MGl2,3) + 110592.*b*pow(MGl2,4) + 
    49152.*b*s12*log(dS)*pow(MGl2,4) - 
    24576.*s12*log((1. + b)*pow(1. - 1.*b,-1))*pow(MGl2,4) + 
    3072.*s12*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(MGl2,4) - 
    24576.*b*s12*log(pow(muR,2)*pow(s12,-1))*pow(MGl2,4) + 
    49152.*b*s12*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(MGl2,4) - 
    1536.*s12*log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(MGl2,4) - 43008.*b*s12*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,4) + 
    21504.*b*s12*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,4) - 
    12288.*b*s12*log(dS)*log(-1.*pow(1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(MGl2,4) + 
    6144.*b*s12*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,4) - 
    43008.*b*s12*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(MGl2,4) + 43008.*b*s12*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(MGl2,4) + 
    12288.*b*s12*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*
     pow(MGl2,4) + 3072.*s12*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*
     pow(b,2)*pow(MGl2,4) - 1536.*s12*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(MGl2,4) + 
    50688.*b*pow(MGl2,2)*pow(s12,2) - 
    165888.*cos(t)*pow(b,2)*pow(MGl2,2)*pow(s12,2) + 
    165888.*pow(b,3)*pow(MGl2,2)*pow(s12,2) + 
    82944.*cos(t)*pow(b,3)*pow(MGl2,2)*pow(s12,2) - 
    165888.*cos(t)*pow(b,4)*pow(MGl2,2)*pow(s12,2) + 
    41472.*pow(b,5)*pow(MGl2,2)*pow(s12,2) + 
    65536.*b*log(dS)*pow(MGl2,3)*pow(s12,2) - 
    32768.*log((1. + b)*pow(1. - 1.*b,-1))*pow(MGl2,3)*pow(s12,2) + 
    49152.*b*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(MGl2,3)*
     pow(s12,2) + 4096.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*
     pow(MGl2,3)*pow(s12,2) - 6144.*b*cos(t)*log(dS)*
     log((1. + b)*pow(1. - 1.*b,-1))*pow(MGl2,3)*pow(s12,2) - 
    32768.*b*log(pow(muR,2)*pow(s12,-1))*pow(MGl2,3)*pow(s12,2) + 
    65536.*b*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(MGl2,3)*pow(s12,2) - 
    2048.*log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(MGl2,3)*pow(s12,2) + 3072.*b*cos(t)*
     log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(MGl2,3)*pow(s12,2) - 40960.*b*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,3)*
     pow(s12,2) + 20480.*b*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,3)*
     pow(s12,2) - 32768.*b*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,3)*
     pow(s12,2) + 16384.*b*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,3)*
     pow(s12,2) - 40960.*b*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(MGl2,3)*pow(s12,2) + 40960.*b*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(MGl2,3)*
     pow(s12,2) + 32768.*b*PolyLog(2.,
      b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(MGl2,3)*pow(s12,2) - 
    98304.*cos(t)*log(dS)*pow(b,2)*pow(MGl2,3)*pow(s12,2) - 
    24576.*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*pow(MGl2,3)*
     pow(s12,2) + 7168.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*
     pow(MGl2,3)*pow(s12,2) + 49152.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,2)*pow(MGl2,3)*pow(s12,2) - 
    98304.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(MGl2,3)*
     pow(s12,2) - 3584.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(MGl2,3)*pow(s12,2) + 
    86016.*cos(t)*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,2)*pow(MGl2,3)*pow(s12,2) - 
    43008.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(MGl2,3)*pow(s12,2) + 24576.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(MGl2,3)*pow(s12,2) - 12288.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(MGl2,3)*pow(s12,2) + 86016.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,2)*pow(MGl2,3)*
     pow(s12,2) - 86016.*cos(t)*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,2)*
     pow(MGl2,3)*pow(s12,2) - 24576.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,2)*
     pow(MGl2,3)*pow(s12,2) + 49152.*log(dS)*pow(b,3)*pow(MGl2,3)*
     pow(s12,2) - 6144.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*
     pow(b,3)*pow(MGl2,3)*pow(s12,2) - 
    24576.*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(MGl2,3)*pow(s12,2) + 
    49152.*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(MGl2,3)*
     pow(s12,2) + 3072.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(MGl2,3)*pow(s12,2) - 
    43008.*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*
     pow(b,3)*pow(MGl2,3)*pow(s12,2) + 
    21504.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,3)*pow(s12,2) - 12288.*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,3)*pow(s12,2) + 6144.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,3)*pow(s12,2) - 43008.*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,3)*pow(MGl2,3)*
     pow(s12,2) + 43008.*PolyLog(2.,
      b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,3)*pow(MGl2,3)*
     pow(s12,2) + 12288.*PolyLog(2.,
      b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,3)*pow(MGl2,3)*
     pow(s12,2) + 3072.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*
     pow(MGl2,3)*pow(s12,2) - 1536.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,4)*pow(MGl2,3)*pow(s12,2) + 
    66816.*b*MGl2*pow(s12,3) - 50688.*MGl2*cos(t)*pow(b,2)*pow(s12,3) + 
    66816.*MGl2*pow(b,3)*pow(s12,3) + 
    41472.*MGl2*cos(t)*pow(b,3)*pow(s12,3) - 
    138240.*MGl2*cos(t)*pow(b,4)*pow(s12,3) + 
    62208.*MGl2*pow(b,5)*pow(s12,3) + 
    41472.*MGl2*cos(t)*pow(b,5)*pow(s12,3) - 
    41472.*MGl2*cos(t)*pow(b,6)*pow(s12,3) + 
    6912.*MGl2*pow(b,7)*pow(s12,3) + 
    55296.*b*log(dS)*pow(MGl2,2)*pow(s12,3) - 
    27648.*log((1. + b)*pow(1. - 1.*b,-1))*pow(MGl2,2)*pow(s12,3) + 
    49152.*b*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(MGl2,2)*
     pow(s12,3) - 8832.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*
     pow(MGl2,2)*pow(s12,3) - 6144.*b*cos(t)*log(dS)*
     log((1. + b)*pow(1. - 1.*b,-1))*pow(MGl2,2)*pow(s12,3) - 
    27648.*b*log(pow(muR,2)*pow(s12,-1))*pow(MGl2,2)*pow(s12,3) + 
    55296.*b*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(MGl2,2)*pow(s12,3) + 
    4416.*log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(MGl2,2)*pow(s12,3) + 3072.*b*cos(t)*
     log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(MGl2,2)*pow(s12,3) - 11520.*b*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,2)*
     pow(s12,3) + 5760.*b*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,2)*
     pow(s12,3) - 26112.*b*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,2)*
     pow(s12,3) + 13056.*b*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(MGl2,2)*
     pow(s12,3) - 11520.*b*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(MGl2,2)*pow(s12,3) + 11520.*b*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(MGl2,2)*
     pow(s12,3) + 26112.*b*PolyLog(2.,
      b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(MGl2,2)*pow(s12,3) - 
    98304.*cos(t)*log(dS)*pow(b,2)*pow(MGl2,2)*pow(s12,3) - 
    43008.*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*pow(MGl2,2)*
     pow(s12,3) - 18432.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*
     pow(MGl2,2)*pow(s12,3) - 3456.*log(dS)*
     log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*pow(MGl2,2)*pow(s12,3) + 
    2304.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*
     pow(MGl2,2)*pow(s12,3) + 49152.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,2)*pow(MGl2,2)*pow(s12,3) - 
    98304.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(MGl2,2)*
     pow(s12,3) + 1728.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(MGl2,2)*pow(s12,3) - 
    1152.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(MGl2,2)*pow(s12,3) + 
    61440.*cos(t)*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,2)*pow(MGl2,2)*pow(s12,3) - 
    30720.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(MGl2,2)*pow(s12,3) + 49152.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(MGl2,2)*pow(s12,3) - 24576.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(MGl2,2)*pow(s12,3) + 61440.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,2)*pow(MGl2,2)*
     pow(s12,3) - 61440.*cos(t)*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,2)*
     pow(MGl2,2)*pow(s12,3) - 49152.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,2)*
     pow(MGl2,2)*pow(s12,3) + 86016.*log(dS)*pow(b,3)*pow(MGl2,2)*
     pow(s12,3) + 36864.*cos(t)*log(dS)*pow(b,3)*pow(MGl2,2)*pow(s12,3) + 
    36864.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,3)*pow(MGl2,2)*
     pow(s12,3) - 10752.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*
     pow(b,3)*pow(MGl2,2)*pow(s12,3) - 
    43008.*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(MGl2,2)*pow(s12,3) - 
    18432.*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(MGl2,2)*
     pow(s12,3) + 86016.*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) + 36864.*cos(t)*log(dS)*
     log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(MGl2,2)*pow(s12,3) + 
    5376.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(MGl2,2)*pow(s12,3) - 
    62976.*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*
     pow(b,3)*pow(MGl2,2)*pow(s12,3) - 
    32256.*cos(t)*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,3)*pow(MGl2,2)*pow(s12,3) + 
    31488.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) + 16128.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) - 33792.*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) - 9216.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) + 16896.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) + 4608.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) - 62976.*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,3)*pow(MGl2,2)*
     pow(s12,3) - 32256.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,3)*pow(MGl2,2)*
     pow(s12,3) + 62976.*PolyLog(2.,
      b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,3)*pow(MGl2,2)*
     pow(s12,3) + 32256.*cos(t)*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) + 33792.*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) + 9216.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,3)*
     pow(MGl2,2)*pow(s12,3) - 73728.*cos(t)*log(dS)*pow(b,4)*pow(MGl2,2)*
     pow(s12,3) - 9216.*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*
     pow(MGl2,2)*pow(s12,3) + 6528.*log(dS)*
     log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*pow(MGl2,2)*pow(s12,3) + 
    2304.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*
     pow(MGl2,2)*pow(s12,3) + 36864.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,4)*pow(MGl2,2)*pow(s12,3) - 
    73728.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,4)*pow(MGl2,2)*
     pow(s12,3) - 3264.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,4)*pow(MGl2,2)*pow(s12,3) - 
    1152.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,4)*pow(MGl2,2)*pow(s12,3) + 
    64512.*cos(t)*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,4)*pow(MGl2,2)*pow(s12,3) - 
    32256.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(MGl2,2)*pow(s12,3) + 18432.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(MGl2,2)*pow(s12,3) - 9216.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(MGl2,2)*pow(s12,3) + 64512.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,4)*pow(MGl2,2)*
     pow(s12,3) - 64512.*cos(t)*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,4)*
     pow(MGl2,2)*pow(s12,3) - 18432.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,4)*
     pow(MGl2,2)*pow(s12,3) + 18432.*log(dS)*pow(b,5)*pow(MGl2,2)*
     pow(s12,3) - 4608.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*
     pow(b,5)*pow(MGl2,2)*pow(s12,3) - 
    9216.*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*pow(MGl2,2)*pow(s12,3) + 
    18432.*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*pow(MGl2,2)*
     pow(s12,3) + 2304.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,5)*pow(MGl2,2)*pow(s12,3) - 
    16128.*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*
     pow(b,5)*pow(MGl2,2)*pow(s12,3) + 
    8064.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(MGl2,2)*pow(s12,3) - 4608.*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(MGl2,2)*pow(s12,3) + 2304.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(MGl2,2)*pow(s12,3) - 16128.*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,5)*pow(MGl2,2)*
     pow(s12,3) + 16128.*PolyLog(2.,
      b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,5)*pow(MGl2,2)*
     pow(s12,3) + 4608.*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*
     pow(b,5)*pow(MGl2,2)*pow(s12,3) + 
    1152.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,6)*pow(MGl2,2)*
     pow(s12,3) - 576.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,6)*pow(MGl2,2)*pow(s12,3) + 
    1008.*b*pow(s12,4) + 18432.*b*MGl2*log(dS)*pow(s12,4) - 
    9216.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*pow(s12,4) + 
    27648.*b*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(s12,4) - 
    4992.*MGl2*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(s12,4) + 
    8832.*b*MGl2*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*
     pow(s12,4) - 9216.*b*MGl2*log(pow(muR,2)*pow(s12,-1))*pow(s12,4) + 
    18432.*b*MGl2*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(s12,4) + 
    2496.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(s12,4) - 4416.*b*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(s12,4) - 
    768.*b*MGl2*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(s12,4) + 
    384.*b*MGl2*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(s12,4) - 
    7680.*b*MGl2*log(dS)*log(-1.*pow(1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(s12,4) + 
    3840.*b*MGl2*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(s12,4) - 
    768.*b*MGl2*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(s12,4) + 
    768.*b*MGl2*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(s12,4) + 7680.*b*MGl2*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(s12,4) - 
    5760.*cos(t)*pow(b,2)*pow(s12,4) - 
    55296.*MGl2*cos(t)*log(dS)*pow(b,2)*pow(s12,4) - 
    26112.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*pow(s12,4) - 
    12288.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*
     pow(s12,4) - 7872.*MGl2*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*
     pow(b,2)*pow(s12,4) + 1536.*MGl2*cos(t)*log(dS)*
     log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*pow(s12,4) + 
    27648.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(s12,4) - 
    55296.*MGl2*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,2)*
     pow(s12,4) + 3936.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(s12,4) - 
    768.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(s12,4) + 
    11520.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(s12,4) - 5760.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(s12,4) + 26112.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(s12,4) - 13056.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(s12,4) + 11520.*MGl2*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,2)*pow(s12,4) - 
    11520.*MGl2*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,2)*pow(s12,4) - 26112.*MGl2*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,2)*pow(s12,4) \
+ 9216.*pow(b,3)*pow(s12,4) + 6336.*cos(t)*pow(b,3)*pow(s12,4) + 
    52224.*MGl2*log(dS)*pow(b,3)*pow(s12,4) + 
    24576.*MGl2*cos(t)*log(dS)*pow(b,3)*pow(s12,4) + 
    36864.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,3)*
     pow(s12,4) + 4224.*MGl2*cos(t)*log(dS)*
     log((1. + b)*pow(1. - 1.*b,-1))*pow(b,3)*pow(s12,4) - 
    26112.*MGl2*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,4) - 
    12288.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,4) + 
    52224.*MGl2*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,4) + 
    24576.*MGl2*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*
     pow(s12,4) - 2112.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,4) - 
    21120.*MGl2*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,3)*pow(s12,4) - 
    15360.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,4) + 10560.*MGl2*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,4) + 7680.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,4) - 25344.*MGl2*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,4) - 12288.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,4) + 12672.*MGl2*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,4) + 6144.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,4) - 21120.*MGl2*PolyLog(2.,
      (b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,3)*pow(s12,4) - 
    15360.*MGl2*cos(t)*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(b,3)*pow(s12,4) + 21120.*MGl2*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,3)*pow(s12,4) \
+ 15360.*MGl2*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,3)*pow(s12,4) + 25344.*MGl2*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,3)*pow(s12,4) \
+ 12288.*MGl2*cos(t)*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*
     pow(b,3)*pow(s12,4) - 26496.*cos(t)*pow(b,4)*pow(s12,4) - 
    73728.*MGl2*cos(t)*log(dS)*pow(b,4)*pow(s12,4) - 
    15360.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*pow(s12,4) - 
    9216.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*pow(s12,4) - 
    960.*MGl2*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*pow(s12,4) + 
    2688.*MGl2*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*
     pow(s12,4) + 36864.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,4)*
     pow(s12,4) - 73728.*MGl2*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,4)*pow(s12,4) + 480.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,4)*pow(s12,4) - 
    1344.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,4)*pow(s12,4) + 
    52224.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(s12,4) - 26112.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(s12,4) + 30720.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(s12,4) - 15360.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(s12,4) + 52224.*MGl2*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,4)*pow(s12,4) - 
    52224.*MGl2*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,4)*pow(s12,4) - 30720.*MGl2*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,4)*pow(s12,4) \
+ 16128.*pow(b,5)*pow(s12,4) + 14688.*cos(t)*pow(b,5)*pow(s12,4) + 
    30720.*MGl2*log(dS)*pow(b,5)*pow(s12,4) + 
    18432.*MGl2*cos(t)*log(dS)*pow(b,5)*pow(s12,4) + 
    9216.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,5)*pow(s12,4) - 
    5760.*MGl2*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,5)*
     pow(s12,4) - 15360.*MGl2*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*
     pow(s12,4) - 9216.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*
     pow(s12,4) + 30720.*MGl2*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*
     pow(s12,4) + 18432.*MGl2*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,5)*pow(s12,4) + 2880.*MGl2*cos(t)*
     log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*
     pow(s12,4) - 23808.*MGl2*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,4) - 16128.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,4) + 11904.*MGl2*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,4) + 8064.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,4) - 10752.*MGl2*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,4) - 4608.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,4) + 5376.*MGl2*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,4) + 2304.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,4) - 23808.*MGl2*PolyLog(2.,
      (b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,5)*pow(s12,4) - 
    16128.*MGl2*cos(t)*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(b,5)*pow(s12,4) + 23808.*MGl2*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,5)*pow(s12,4) \
+ 16128.*MGl2*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,5)*pow(s12,4) + 10752.*MGl2*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,5)*pow(s12,4) \
+ 4608.*MGl2*cos(t)*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*
     pow(b,5)*pow(s12,4) - 24192.*cos(t)*pow(b,6)*pow(s12,4) - 
    18432.*MGl2*cos(t)*log(dS)*pow(b,6)*pow(s12,4) - 
    1536.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,6)*pow(s12,4) + 
    2112.*MGl2*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,6)*
     pow(s12,4) + 1152.*MGl2*cos(t)*log(dS)*
     log((1. + b)*pow(1. - 1.*b,-1))*pow(b,6)*pow(s12,4) + 
    9216.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,6)*pow(s12,4) - 
    18432.*MGl2*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,6)*
     pow(s12,4) - 1056.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,6)*pow(s12,4) - 
    576.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,6)*pow(s12,4) + 
    16128.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,6)*
     pow(s12,4) - 8064.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,6)*
     pow(s12,4) + 4608.*MGl2*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,6)*
     pow(s12,4) - 2304.*MGl2*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,6)*
     pow(s12,4) + 16128.*MGl2*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,6)*pow(s12,4) - 
    16128.*MGl2*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,6)*pow(s12,4) - 4608.*MGl2*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,6)*pow(s12,4) \
+ 6912.*pow(b,7)*pow(s12,4) + 5184.*cos(t)*pow(b,7)*pow(s12,4) + 
    3072.*MGl2*log(dS)*pow(b,7)*pow(s12,4) - 
    1152.*MGl2*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,7)*
     pow(s12,4) - 1536.*MGl2*log(pow(muR,2)*pow(s12,-1))*pow(b,7)*
     pow(s12,4) + 3072.*MGl2*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,7)*
     pow(s12,4) + 576.*MGl2*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,7)*pow(s12,4) - 
    2688.*MGl2*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,7)*pow(s12,4) + 
    1344.*MGl2*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,4) - 768.*MGl2*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,4) + 384.*MGl2*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,4) - 2688.*MGl2*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(b,7)*pow(s12,4) + 2688.*MGl2*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,7)*pow(s12,4) \
+ 768.*MGl2*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,7)*
     pow(s12,4) - 3456.*cos(t)*pow(b,8)*pow(s12,4) + 
    192.*MGl2*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,8)*pow(s12,4) - 
    96.*MGl2*log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(b,8)*pow(s12,4) + 432.*pow(b,9)*pow(s12,4) + 
    1984.*b*log(dS)*pow(s12,5) - 
    992.*log((1. + b)*pow(1. - 1.*b,-1))*pow(s12,5) + 
    4608.*b*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(s12,5) - 
    644.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(s12,5) + 
    2496.*b*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(s12,5) - 
    992.*b*log(pow(muR,2)*pow(s12,-1))*pow(s12,5) + 
    1984.*b*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(s12,5) + 
    322.*log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(s12,5) - 1248.*b*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(s12,5) + 
    56.*b*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*
     pow(s12,5) - 28.*b*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(s12,5) - 
    752.*b*log(dS)*log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*
     pow(s12,5) + 376.*b*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(s12,5) + 
    56.*b*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(s12,5) - 
    56.*b*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(s12,5) + 
    752.*b*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(s12,5) - 
    9216.*cos(t)*log(dS)*pow(b,2)*pow(s12,5) - 
    5760.*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*pow(s12,5) - 
    3456.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*pow(s12,5) - 
    2996.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*pow(s12,5) - 
    1104.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,2)*
     pow(s12,5) + 4608.*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,2)*
     pow(s12,5) - 9216.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,2)*pow(s12,5) + 1498.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(s12,5) + 
    552.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,2)*pow(s12,5) + 
    384.*cos(t)*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,2)*pow(s12,5) - 
    192.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(s12,5) + 3840.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(s12,5) - 1920.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,2)*
     pow(s12,5) + 384.*cos(t)*PolyLog(2.,
      (b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,2)*pow(s12,5) - 
    384.*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,2)*pow(s12,5) - 3840.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,2)*pow(s12,5) \
+ 11520.*log(dS)*pow(b,3)*pow(s12,5) + 
    6912.*cos(t)*log(dS)*pow(b,3)*pow(s12,5) + 
    11008.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,3)*pow(s12,5) + 
    4192.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,3)*
     pow(s12,5) - 5760.*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,5) - 
    3456.*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,5) + 
    11520.*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,5) + 
    6912.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,5) - 
    2096.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,3)*pow(s12,5) - 
    1632.*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*
     pow(b,3)*pow(s12,5) - 1440.*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,5) + 816.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,5) + 720.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,5) - 5184.*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,5) - 3264.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,5) + 2592.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,5) + 1632.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,3)*
     pow(s12,5) - 1632.*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(b,3)*pow(s12,5) - 1440.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,3)*pow(s12,5) + 
    1632.*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,3)*
     pow(s12,5) + 1440.*cos(t)*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,3)*pow(s12,5) \
+ 5184.*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,3)*
     pow(s12,5) + 3264.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,3)*pow(s12,5) \
- 22016.*cos(t)*log(dS)*pow(b,4)*pow(s12,5) - 
    5376.*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*pow(s12,5) - 
    4032.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*pow(s12,5) - 
    2448.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*pow(s12,5) - 
    600.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,4)*
     pow(s12,5) + 11008.*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,4)*
     pow(s12,5) - 22016.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,4)*pow(s12,5) + 1224.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,4)*pow(s12,5) + 
    300.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,4)*pow(s12,5) + 
    8000.*cos(t)*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,4)*pow(s12,5) - 
    4000.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(s12,5) + 10624.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(s12,5) - 5312.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,4)*
     pow(s12,5) + 8000.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,4)*pow(s12,5) - 
    8000.*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,4)*pow(s12,5) - 10624.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,4)*pow(s12,5) \
+ 10752.*log(dS)*pow(b,5)*pow(s12,5) + 
    8064.*cos(t)*log(dS)*pow(b,5)*pow(s12,5) + 
    6144.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,5)*pow(s12,5) + 
    928.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,5)*
     pow(s12,5) - 5376.*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*pow(s12,5) - 
    4032.*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*pow(s12,5) + 
    10752.*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*pow(s12,5) + 
    8064.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,5)*pow(s12,5) - 
    464.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,5)*pow(s12,5) - 
    5568.*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*
     pow(b,5)*pow(s12,5) - 5520.*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,5) + 2784.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,5) + 2760.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,5) - 4992.*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,5) - 3552.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,5) + 2496.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,5) + 1776.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,5)*
     pow(s12,5) - 5568.*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(b,5)*pow(s12,5) - 5520.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,5)*pow(s12,5) + 
    5568.*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,5)*
     pow(s12,5) + 5520.*cos(t)*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,5)*pow(s12,5) \
+ 4992.*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,5)*
     pow(s12,5) + 3552.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,5)*pow(s12,5) \
- 12288.*cos(t)*log(dS)*pow(b,6)*pow(s12,5) - 
    1664.*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,6)*pow(s12,5) - 
    1152.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,6)*pow(s12,5) + 
    112.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,6)*pow(s12,5) + 
    648.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,6)*
     pow(s12,5) + 6144.*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,6)*
     pow(s12,5) - 12288.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,6)*pow(s12,5) - 56.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,6)*pow(s12,5) - 
    324.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,6)*pow(s12,5) + 
    9216.*cos(t)*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,6)*pow(s12,5) - 
    4608.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,6)*
     pow(s12,5) + 4608.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,6)*
     pow(s12,5) - 2304.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,6)*
     pow(s12,5) + 9216.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,6)*pow(s12,5) - 
    9216.*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,6)*pow(s12,5) - 4608.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,6)*pow(s12,5) \
+ 3328.*log(dS)*pow(b,7)*pow(s12,5) + 
    2304.*cos(t)*log(dS)*pow(b,7)*pow(s12,5) + 
    768.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,7)*pow(s12,5) - 
    864.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,7)*
     pow(s12,5) - 1664.*log(pow(muR,2)*pow(s12,-1))*pow(b,7)*pow(s12,5) - 
    1152.*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,7)*pow(s12,5) + 
    3328.*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,7)*pow(s12,5) + 
    2304.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,7)*pow(s12,5) + 
    432.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,7)*pow(s12,5) - 
    2656.*log(dS)*log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*
     pow(b,7)*pow(s12,5) - 2016.*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,5) + 1328.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,5) + 1008.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,5) - 1088.*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,5) - 576.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,5) + 544.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,5) + 288.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,7)*
     pow(s12,5) - 2656.*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(b,7)*pow(s12,5) - 2016.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,7)*pow(s12,5) + 
    2656.*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,7)*
     pow(s12,5) + 2016.*cos(t)*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,7)*pow(s12,5) \
+ 1088.*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,7)*
     pow(s12,5) + 576.*cos(t)*PolyLog(2.,
      b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,7)*pow(s12,5) - 
    1536.*cos(t)*log(dS)*pow(b,8)*pow(s12,5) - 
    96.*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,8)*pow(s12,5) + 
    220.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,8)*pow(s12,5) + 
    144.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,8)*
     pow(s12,5) + 768.*cos(t)*log(pow(muR,2)*pow(s12,-1))*pow(b,8)*
     pow(s12,5) - 1536.*cos(t)*log(dS)*log(pow(muR,2)*pow(s12,-1))*
     pow(b,8)*pow(s12,5) - 110.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,8)*pow(s12,5) - 
    72.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(b,8)*pow(s12,5) + 1344.*cos(t)*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,8)*
     pow(s12,5) - 672.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,8)*
     pow(s12,5) + 384.*cos(t)*log(dS)*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,8)*
     pow(s12,5) - 192.*cos(t)*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,8)*
     pow(s12,5) + 1344.*cos(t)*
     PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*pow(b,8)*pow(s12,5) - 
    1344.*cos(t)*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*
     pow(b,8)*pow(s12,5) - 384.*cos(t)*
     PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,8)*pow(s12,5) \
+ 192.*log(dS)*pow(b,9)*pow(s12,5) - 
    96.*cos(t)*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,9)*
     pow(s12,5) - 96.*log(pow(muR,2)*pow(s12,-1))*pow(b,9)*pow(s12,5) + 
    192.*log(dS)*log(pow(muR,2)*pow(s12,-1))*pow(b,9)*pow(s12,5) + 
    48.*cos(t)*log((1. + b)*pow(1. - 1.*b,-1))*log(pow(muR,2)*pow(s12,-1))*
     pow(b,9)*pow(s12,5) - 168.*log(dS)*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,9)*
     pow(s12,5) + 84.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(-1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,9)*
     pow(s12,5) - 48.*log(dS)*log(-1.*pow(1. + b*cos(t),2)*
       pow(-1. + pow(b,2),-1))*pow(b,9)*pow(s12,5) + 
    24.*log(pow(muR,2)*pow(s12,-1))*
     log(-1.*pow(1. + b*cos(t),2)*pow(-1. + pow(b,2),-1))*pow(b,9)*
     pow(s12,5) - 168.*PolyLog(2.,(b - 1.*b*cos(t))*pow(-1. + b,-1))*
     pow(b,9)*pow(s12,5) + 168.*
     PolyLog(2.,b*(1. + cos(t))*pow(-1. + b*cos(t),-1))*pow(b,9)*pow(s12,5) \
+ 48.*PolyLog(2.,b*(-1. + cos(t))*pow(1. + b*cos(t),-1))*pow(b,9)*
     pow(s12,5) + 12.*log(dS)*log((1. + b)*pow(1. - 1.*b,-1))*pow(b,10)*
     pow(s12,5) - 6.*log((1. + b)*pow(1. - 1.*b,-1))*
     log(pow(muR,2)*pow(s12,-1))*pow(b,10)*pow(s12,5) - 
    49152.*b*s12*pow(MGl2,4)*pow(log(dS),2) - 
    65536.*b*pow(MGl2,3)*pow(s12,2)*pow(log(dS),2) + 
    98304.*cos(t)*pow(b,2)*pow(MGl2,3)*pow(s12,2)*pow(log(dS),2) - 
    49152.*pow(b,3)*pow(MGl2,3)*pow(s12,2)*pow(log(dS),2) - 
    55296.*b*pow(MGl2,2)*pow(s12,3)*pow(log(dS),2) + 
    98304.*cos(t)*pow(b,2)*pow(MGl2,2)*pow(s12,3)*pow(log(dS),2) - 
    86016.*pow(b,3)*pow(MGl2,2)*pow(s12,3)*pow(log(dS),2) - 
    36864.*cos(t)*pow(b,3)*pow(MGl2,2)*pow(s12,3)*pow(log(dS),2) + 
    73728.*cos(t)*pow(b,4)*pow(MGl2,2)*pow(s12,3)*pow(log(dS),2) - 
    18432.*pow(b,5)*pow(MGl2,2)*pow(s12,3)*pow(log(dS),2) - 
    18432.*b*MGl2*pow(s12,4)*pow(log(dS),2) + 
    55296.*MGl2*cos(t)*pow(b,2)*pow(s12,4)*pow(log(dS),2) - 
    52224.*MGl2*pow(b,3)*pow(s12,4)*pow(log(dS),2) - 
    24576.*MGl2*cos(t)*pow(b,3)*pow(s12,4)*pow(log(dS),2) + 
    73728.*MGl2*cos(t)*pow(b,4)*pow(s12,4)*pow(log(dS),2) - 
    30720.*MGl2*pow(b,5)*pow(s12,4)*pow(log(dS),2) - 
    18432.*MGl2*cos(t)*pow(b,5)*pow(s12,4)*pow(log(dS),2) + 
    18432.*MGl2*cos(t)*pow(b,6)*pow(s12,4)*pow(log(dS),2) - 
    3072.*MGl2*pow(b,7)*pow(s12,4)*pow(log(dS),2) - 
    1984.*b*pow(s12,5)*pow(log(dS),2) + 
    9216.*cos(t)*pow(b,2)*pow(s12,5)*pow(log(dS),2) - 
    11520.*pow(b,3)*pow(s12,5)*pow(log(dS),2) - 
    6912.*cos(t)*pow(b,3)*pow(s12,5)*pow(log(dS),2) + 
    22016.*cos(t)*pow(b,4)*pow(s12,5)*pow(log(dS),2) - 
    10752.*pow(b,5)*pow(s12,5)*pow(log(dS),2) - 
    8064.*cos(t)*pow(b,5)*pow(s12,5)*pow(log(dS),2) + 
    12288.*cos(t)*pow(b,6)*pow(s12,5)*pow(log(dS),2) - 
    3328.*pow(b,7)*pow(s12,5)*pow(log(dS),2) - 
    2304.*cos(t)*pow(b,7)*pow(s12,5)*pow(log(dS),2) + 
    1536.*cos(t)*pow(b,8)*pow(s12,5)*pow(log(dS),2) - 
    192.*pow(b,9)*pow(s12,5)*pow(log(dS),2) - 
    768.*s12*pow(MGl2,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    13824.*b*s12*pow(MGl2,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    768.*s12*pow(b,2)*pow(MGl2,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    1024.*pow(MGl2,3)*pow(s12,2)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    18432.*b*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    1536.*b*cos(t)*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    1792.*pow(b,2)*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    27648.*cos(t)*pow(b,2)*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    13824.*pow(b,3)*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    1536.*cos(t)*pow(b,3)*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    768.*pow(b,4)*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    2208.*pow(MGl2,2)*pow(s12,3)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    9408.*b*pow(MGl2,2)*pow(s12,3)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    1536.*b*cos(t)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    864.*pow(b,2)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    28224.*cos(t)*pow(b,2)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    24192.*pow(b,3)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    13056.*cos(t)*pow(b,3)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    1632.*pow(b,4)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    21312.*cos(t)*pow(b,4)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    5184.*pow(b,5)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    1152.*cos(t)*pow(b,5)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    288.*pow(b,6)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    1248.*MGl2*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    2112.*b*MGl2*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    2208.*b*MGl2*cos(t)*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    1968.*MGl2*pow(b,2)*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    9792.*MGl2*cos(t)*pow(b,2)*pow(s12,4)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    11616.*MGl2*pow(b,3)*pow(s12,4)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    5856.*MGl2*cos(t)*pow(b,3)*pow(s12,4)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    240.*MGl2*pow(b,4)*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    21408.*MGl2*cos(t)*pow(b,4)*pow(s12,4)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    8640.*MGl2*pow(b,5)*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    6624.*MGl2*cos(t)*pow(b,5)*pow(s12,4)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    528.*MGl2*pow(b,6)*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    5472.*MGl2*cos(t)*pow(b,6)*pow(s12,4)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    864.*MGl2*pow(b,7)*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    288.*MGl2*cos(t)*pow(b,7)*pow(s12,4)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    48.*MGl2*pow(b,8)*pow(s12,4)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    161.*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    174.*b*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    624.*b*cos(t)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    749.*pow(b,2)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    780.*cos(t)*pow(b,2)*pow(s12,5)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    1704.*pow(b,3)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    128.*cos(t)*pow(b,3)*pow(s12,5)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    612.*pow(b,4)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    4506.*cos(t)*pow(b,4)*pow(s12,5)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    2640.*pow(b,5)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    2036.*cos(t)*pow(b,5)*pow(s12,5)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    28.*pow(b,6)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    3618.*cos(t)*pow(b,6)*pow(s12,5)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    936.*pow(b,7)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    864.*cos(t)*pow(b,7)*pow(s12,5)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    55.*pow(b,8)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    468.*cos(t)*pow(b,8)*pow(s12,5)*
     pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    54.*pow(b,9)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) + 
    24.*cos(t)*pow(b,9)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    3.*pow(b,10)*pow(s12,5)*pow(log((1. + b)*pow(1. - 1.*b,-1)),2) - 
    12288.*b*s12*pow(MGl2,4)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    16384.*b*pow(MGl2,3)*pow(s12,2)*pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    24576.*cos(t)*pow(b,2)*pow(MGl2,3)*pow(s12,2)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    12288.*pow(b,3)*pow(MGl2,3)*pow(s12,2)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    13824.*b*pow(MGl2,2)*pow(s12,3)*pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    24576.*cos(t)*pow(b,2)*pow(MGl2,2)*pow(s12,3)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    21504.*pow(b,3)*pow(MGl2,2)*pow(s12,3)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    9216.*cos(t)*pow(b,3)*pow(MGl2,2)*pow(s12,3)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    18432.*cos(t)*pow(b,4)*pow(MGl2,2)*pow(s12,3)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    4608.*pow(b,5)*pow(MGl2,2)*pow(s12,3)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    4608.*b*MGl2*pow(s12,4)*pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    13824.*MGl2*cos(t)*pow(b,2)*pow(s12,4)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    13056.*MGl2*pow(b,3)*pow(s12,4)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    6144.*MGl2*cos(t)*pow(b,3)*pow(s12,4)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    18432.*MGl2*cos(t)*pow(b,4)*pow(s12,4)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    7680.*MGl2*pow(b,5)*pow(s12,4)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    4608.*MGl2*cos(t)*pow(b,5)*pow(s12,4)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    4608.*MGl2*cos(t)*pow(b,6)*pow(s12,4)*
     pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    768.*MGl2*pow(b,7)*pow(s12,4)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    496.*b*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    2304.*cos(t)*pow(b,2)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    2880.*pow(b,3)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    1728.*cos(t)*pow(b,3)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    5504.*cos(t)*pow(b,4)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    2688.*pow(b,5)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    2016.*cos(t)*pow(b,5)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    3072.*cos(t)*pow(b,6)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    832.*pow(b,7)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    576.*cos(t)*pow(b,7)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) + 
    384.*cos(t)*pow(b,8)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    48.*pow(b,9)*pow(s12,5)*pow(log(pow(muR,2)*pow(s12,-1)),2) - 
    21504.*b*s12*pow(MGl2,4)*pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    20480.*b*pow(MGl2,3)*pow(s12,2)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    43008.*cos(t)*pow(b,2)*pow(MGl2,3)*pow(s12,2)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    21504.*pow(b,3)*pow(MGl2,3)*pow(s12,2)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    5760.*b*pow(MGl2,2)*pow(s12,3)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    30720.*cos(t)*pow(b,2)*pow(MGl2,2)*pow(s12,3)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    31488.*pow(b,3)*pow(MGl2,2)*pow(s12,3)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    16128.*cos(t)*pow(b,3)*pow(MGl2,2)*pow(s12,3)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    32256.*cos(t)*pow(b,4)*pow(MGl2,2)*pow(s12,3)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    8064.*pow(b,5)*pow(MGl2,2)*pow(s12,3)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    384.*b*MGl2*pow(s12,4)*pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    5760.*MGl2*cos(t)*pow(b,2)*pow(s12,4)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    10560.*MGl2*pow(b,3)*pow(s12,4)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    7680.*MGl2*cos(t)*pow(b,3)*pow(s12,4)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    26112.*MGl2*cos(t)*pow(b,4)*pow(s12,4)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    11904.*MGl2*pow(b,5)*pow(s12,4)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    8064.*MGl2*cos(t)*pow(b,5)*pow(s12,4)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    8064.*MGl2*cos(t)*pow(b,6)*pow(s12,4)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    1344.*MGl2*pow(b,7)*pow(s12,4)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    28.*b*pow(s12,5)*pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    192.*cos(t)*pow(b,2)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    816.*pow(b,3)*pow(s12,5)*pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    720.*cos(t)*pow(b,3)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    4000.*cos(t)*pow(b,4)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    2784.*pow(b,5)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    2760.*cos(t)*pow(b,5)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    4608.*cos(t)*pow(b,6)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    1328.*pow(b,7)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    1008.*cos(t)*pow(b,7)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) + 
    672.*cos(t)*pow(b,8)*pow(s12,5)*
     pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    84.*pow(b,9)*pow(s12,5)*pow(log((-1. + b)*pow(-1. + b*cos(t),-1)),2) - 
    6144.*b*s12*pow(MGl2,4)*pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    16384.*b*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    12288.*cos(t)*pow(b,2)*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    6144.*pow(b,3)*pow(MGl2,3)*pow(s12,2)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    13056.*b*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    24576.*cos(t)*pow(b,2)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    16896.*pow(b,3)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    4608.*cos(t)*pow(b,3)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    9216.*cos(t)*pow(b,4)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    2304.*pow(b,5)*pow(MGl2,2)*pow(s12,3)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    3840.*b*MGl2*pow(s12,4)*pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    13056.*MGl2*cos(t)*pow(b,2)*pow(s12,4)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    12672.*MGl2*pow(b,3)*pow(s12,4)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    6144.*MGl2*cos(t)*pow(b,3)*pow(s12,4)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    15360.*MGl2*cos(t)*pow(b,4)*pow(s12,4)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    5376.*MGl2*pow(b,5)*pow(s12,4)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    2304.*MGl2*cos(t)*pow(b,5)*pow(s12,4)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    2304.*MGl2*cos(t)*pow(b,6)*pow(s12,4)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    384.*MGl2*pow(b,7)*pow(s12,4)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    376.*b*pow(s12,5)*pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    1920.*cos(t)*pow(b,2)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    2592.*pow(b,3)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    1632.*cos(t)*pow(b,3)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    5312.*cos(t)*pow(b,4)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    2496.*pow(b,5)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    1776.*cos(t)*pow(b,5)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    2304.*cos(t)*pow(b,6)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    544.*pow(b,7)*pow(s12,5)*pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),
      2) - 288.*cos(t)*pow(b,7)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) + 
    192.*cos(t)*pow(b,8)*pow(s12,5)*
     pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    24.*pow(b,9)*pow(s12,5)*pow(log((1. - 1.*b)*pow(1. + b*cos(t),-1)),2) - 
    4.*s12*PolyLog(2.,2.*b*pow(1. + b,-1))*(1. + pow(b,2))*
     (40.*MGl2*s12 + 24.*MGl2*s12*pow(b,2) - 
       4.*b*s12*cos(t)*(12.*MGl2 + 5.*s12 + 3.*s12*pow(b,2)) + 
       48.*pow(MGl2,2) - 161.*pow(s12,2) + 16.*pow(b,2)*pow(s12,2) + 
       6.*cos(t)*pow(b,2)*pow(s12,2) + 3.*pow(b,4)*pow(s12,2))*
     pow(4.*MGl2 + s12 - 2.*b*s12*cos(t) + s12*pow(b,2),2) - 
    16.*b*s12*PolyLog(2.,b*(1. + cos(t))*pow(-1. + b,-1))*
     (104.*MGl2*s12 + 24.*MGl2*s12*pow(b,2) - 
       4.*b*s12*cos(t)*(12.*MGl2 + 13.*s12 + 3.*s12*pow(b,2)) + 
       48.*pow(MGl2,2) + 47.*pow(s12,2) + 32.*pow(b,2)*pow(s12,2) + 
       6.*cos(t)*pow(b,2)*pow(s12,2) + 3.*pow(b,4)*pow(s12,2))*
     pow(4.*MGl2 + s12 - 2.*b*s12*cos(t) + s12*pow(b,2),2))*pow(sin(t),3);

}
