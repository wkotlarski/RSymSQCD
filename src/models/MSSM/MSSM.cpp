#include "models/MSSM.hpp"
#include "constants.hpp"

#include <boost/math/special_functions/pow.hpp>
#include "Li2.hpp"

#include "clooptools.h"

using boost::math::pow;

double inline Power(double x, int i) {
   return pow(x, i);
}

MSSM::MSSM(MSSMParameters const& params)
   : MassTop {params.MassTop},
     MassGlu {params.MassGlu},
     MassSq {params.MassSq},
     eta {{std::sqrt(1.+pow<2>(params.delta)), 0., params.delta, params.eta_sign}},
     WidthGlu {params.WidthGlu}
   {};

#include <iostream>
double MSSM::matrixTree_uubar_suLsuLdagger(const double alphas, const double s, const double t) const {
   return (2*pow<2>(alphas)*pi*(6 + (s*(2*pow<2>(MassGlu) + 3*s - 2*t))/pow<2>(MassGlu*MassGlu - t))*(-pow<2>(pow<2>(MassSq) - t) - s*t))/(27.*pow<4>(s));
}
double MSSM::matrixTree_ddbar_suLsuLdagger(const double alphas, const double s, const double t) const {
   return (-4*(alphas*alphas)*pi*(Power(MassSq*MassSq - t,2) + s*t))/(9.*Power(s,4));
}
double MSSM::matrixTree_gg_suLsuLdagger(const double alphas, const double s, const double t) const {
}
double MSSM::sigmaTree_uubar_suLsuLdagger(double alphas, double s, double Dminus4) const {
   return (-2*(alphas*alphas)*pi*(2*sqrt(s*(-4*(MassSq*MassSq) + s))*(-(MassGlu*MassGlu) + 3*(MassSq*MassSq) + 2*s) + 2*(2*pow<4>(MassGlu) + 2*pow<4>(MassSq) + 6*(MassSq*MassSq)*s - 3*(s*s) - 4*(MassGlu*MassGlu)*(MassSq*MassSq + s))*atanh(sqrt(s*(-4*(MassSq*MassSq) + s))/(2*(MassGlu*MassGlu) - 2*(MassSq*MassSq) + s))))/(27.*pow<3>(s));
}
double MSSM::sigmaTree_ddbar_suLsuLdagger(double alphas, double s, double Dminus4) const {
   return (2*(alphas*alphas)*pi*Power(s*(-4*(MassSq*MassSq) + s),1.5))/(27.*Power(s,4));
}
double MSSM::sigmaTree_gg_suLsuLdagger(double alphas, double s, double Dminus4) const {
}

#include "mssm_uu_suLsuR_virt_matrix.cpp"
#include "mssm_uu_suLsuRg_soft.cpp"
#include "mssm_uu_suLsuRg_hard.cpp"

#include "mssm_uubar_suLsuLdagger_virt_matrix.cpp"
#include "mssm_uubar_suLsuLdaggerg_soft.cpp"
#include "mssm_uubar_suLsuLdaggerg_hard.cpp"

#include "mssm_ddbar_suLsuLdagger_virt_matrix.cpp"
#include "mssm_ddbar_suLsuLdaggerg_soft.cpp"
#include "mssm_ddbar_suLsuLdaggerg_hard.cpp"
//#include "mssm_ddbar_suLsuLdagger_virt_matrix.cpp"
//#include "mssm_gg_suLsuLdagger_virt_matrix.cpp"
