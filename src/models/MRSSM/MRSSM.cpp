#include "models/MRSSM.hpp"
#include "constants.hpp"

#include <boost/math/special_functions/pow.hpp>
#include "spdlog/spdlog.h"

#include "clooptools.h"
#include "Li2.hpp"

#include <cmath>

using boost::math::pow;

MRSSM::MRSSM(MRSSMParameters const& params)
   : MassTop {params.MassTop},
     MassGlu {params.MassGlu},
     MasssigmaO {params.MasssigmaO},
     MassphiO {std::sqrt(pow<2>(params.MasssigmaO) + 4*pow<2>(params.MassGlu))},
     MassSq {params.MassSq},
     eta {{std::sqrt(1.+pow<2>(params.delta)), 0., params.delta, params.eta_sign}},
     WidthGlu {params.WidthGlu}
   {};

double MRSSM::matrixTree_udbar_suLsdLdagger(double alphas, double S, double T, int Dminus4Power) const {
   double U = 2*pow<2>(MassSq) - S - T;
   return  (315.82734083485946*pow<2>(alphas)*(-pow<4>(MassSq) + T*U))/pow<2>(-(pow<2>(MassGlu)) + T);
}

#include "mrssm_uu_suLsuR_born_matrix.cpp"
#include "mrssm_uubar_suLsuLdagger_born_matrix.cpp"
#include "mrssm_ddbar_suLsuLdagger_born_matrix.cpp"
#include "mrssm_gg_suLsuLdagger_born_matrix.cpp"

double MRSSM::matrixTree_uubar_glglbar(double alphaS, double S, double T, int Dminus4Power) const
{
   const double U = 2*pow<2>(MassGlu) - S - T;
   const double MsquaredReal = 0.5*16.*52.63789013914324*(alphaS*alphaS)*(36*pow<2>(MassGlu)/S + 18*pow<2>(pow<2>(MassGlu) - T)/pow<2>(S) + (4*pow<2>(pow<2>(MassGlu) - T))/pow<2>(pow<2>(MassSq) - T) + (9.*(pow<2>(MassGlu)))/(-(pow<2>(MassSq)) + T) + (9.*pow<2>(pow<2>(MassGlu) - T))/(S*(-(pow<2>(MassSq)) + T)) + (18.*pow<2>(pow<2>(MassGlu) - U))/(pow<2>(S)) + (4.*pow<2>(pow<2>(MassGlu) - U))/pow<2>(pow<2>(MassSq) - U) + (9.*(pow<2>(MassGlu)))/(-(pow<2>(MassSq)) + U) + (9.*pow<2>(pow<2>(MassGlu) - U))/(S*(-(pow<2>(MassSq)) + U)))
;
   return MsquaredReal/18.;
}

double MRSSM::matrixTree_gg_glglbar(double alphaS, double S, double T, int Dminus4Power) const
{
   const double U = 2*pow<2>(MassGlu) - S - T;
   const double MsquaredReal = 8.*pow<2>(4.*pi*alphaS)*3.*24.*(1.-(U - pow<2>(MassGlu))*(T - pow<2>(MassGlu))/(pow<2>(S)))*(pow<2>(S)/((T - pow<2>(MassGlu))*(U - pow<2>(MassGlu))) - 2. + 4.*pow<2>(MassGlu)*S/((T - pow<2>(MassGlu))*(U - pow<2>(MassGlu)))*(1. - pow<2>(MassGlu)*S/((T - pow<2>(MassGlu))*(U - pow<2>(MassGlu)))));

return MsquaredReal/256.;
}

double MRSSM::sigmaTree_uubar_suLsuLdagger(double alphas, double s12, double Dminus4) const {
   const double Alfas2 = pow<2>(alphas);
   return (1+0.5*Dminus4)*(Alfas2*pi*(8*(-3*pow<2>(MassSq) + MassGlu*MassGlu - 2*s12)*sqrt(s12*(-4*pow<2>(MassSq) + s12)) +
       8*(2*pow<4>(MassSq) + 2*pow<4>(MassGlu) - 4*pow<2>(MassGlu)*s12 - 3*(s12*s12) + pow<2>(MassSq)*(-4*pow<2>(MassGlu) + 6*s12))*
        atanh(sqrt(s12*(-4*pow<2>(MassSq) + s12))/(2*pow<2>(MassSq) - 2*pow<2>(MassGlu) - s12))))/(54.*pow<3>(s12));
}

// checked against MG with SUSYQCD model
double MRSSM::sigmaTree_ddbar_suLsuLdagger(double alphas, double s12, double Dminus4 ) const {
   const double Alfas2 = pow<2>(alphas);
   return (1+0.5*Dminus4)*(2*Alfas2*pi*pow(-4*pow<2>(MassSq) + s12,1.5))/(27.*pow(s12,2.5));
}

double MRSSM::sigmaTree_gg_suLsuLdagger(double alphas, double s, double Dminus4) const {
   const double MassSq2 = pow<2>(MassSq);
   const double coeff1 =
      (alphas*alphas*pi*(62*(MassSq*MassSq)*sqrt(s*(-4*(MassSq*MassSq) + s)) +
                         5*sqrt(pow<3>(s)*(-4*(MassSq*MassSq) + s)) -
                         16*(MassSq*MassSq)*(MassSq*MassSq +
                            4*s)*atanh(sqrt(1 -
                                  (4*(MassSq*MassSq))/s))))/(48.*pow<3>(s));
   const double coeffDm4 =
      -1/96.*(alphas*alphas*pi*(130*MassSq2*sqrt(s*(-4*MassSq2 + s)) + 5*sqrt(pow<3>(s)*(-4*MassSq2 + s)) - 32*MassSq2*(MassSq2 + 4*s)*atanh(sqrt(1 - (4*MassSq2)/s))))/pow<3>(s);
   return (1+0.5*Dminus4)*coeff1 + Dminus4*coeffDm4;
}

/*
 *    checked with MadGraph
 */
double MRSSM::sigmaTree_uu_suLsuR(double alphas, double s, double Dminus4) const {
   const double MGl2 = pow<2>(MassGlu);
   const double a = pow<2>(alphas);
   return (1+0.5*Dminus4)*(-4.*a*pi*(2.*sqrt(s*(-4*pow<2>(MassSq) + s)) + (2.*pow<2>(MassSq) - 2.*MGl2 - s)*
      log((4.*MGl2 + pow<2>(1. + sqrt(1. - (4.*pow<2>(MassSq))/s))*s)/
      (4.*MGl2 + pow<2>(-1. + sqrt(1. - (4.*pow<2>(MassSq))/s))*s))))/(9.*pow<2>(s));
}

/*
 *  pp -> suL suR
 */

#include "mrssm_uu_suLsuR_virt_matrix.cpp"
#include "mrssm_uu_suLsuRg_soft_dp.cpp"
#include "mrssm_uu_suLsuRg_soft_sp.cpp"
#include "mrssm_uu_suLsuRg_soft_finite.cpp"
#include "mrssm_uu_suLsuRg_hard.cpp"
#include "mrssm_gu_suLsuRubar_hard.cpp"
#include "mrssm_gu_suLsuRubar_hard-DR.cpp"
#include "mrssm_gu_suLsuRubar_hard-DS.cpp"

/*
 *  pp -> suL suL*
 */
#include "mrssm_uubar_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_uubar_suLsuLdaggerg_soft_dp.cpp"
#include "mrssm_uubar_suLsuLdaggerg_soft_sp.cpp"
#include "mrssm_uubar_suLsuLdaggerg_soft_finite.cpp"
#include "mrssm_uubar_suLsuLdaggerg_hard.cpp"

#include "mrssm_ddbar_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_ddbar_suLsuLdaggerg_soft_dp.cpp"
#include "mrssm_ddbar_suLsuLdaggerg_soft_sp.cpp"
#include "mrssm_ddbar_suLsuLdaggerg_soft_finite.cpp"
#include "mrssm_ddbar_suLsuLdaggerg_hard.cpp"

#include "mrssm_gg_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_gg_suLsuLdaggerg_soft_dp.cpp"
#include "mrssm_gg_suLsuLdaggerg_soft_sp.cpp"
#include "mrssm_gg_suLsuLdaggerg_soft_finite.cpp"
#include "mrssm_gg_suLsuLdaggerg_hard.cpp"

#include "mrssm_gu_suLsuLdaggeru_hard.cpp"
#include "mrssm_gu_suLsuLdaggeru_hard-DS.cpp"
#include "mrssm_gu_suLsuLdaggeru_hard-DR.cpp"

#include "mrssm_gubar_suLsuLdaggerubar_hard.cpp"
#include "mrssm_gubar_suLsuLdaggerubar_hard-DR.cpp"

#include "mrssm_gd_suLsuLdaggerd_hard.cpp"

#include "mrssm_gdbar_suLsuLdaggerdbar_hard.cpp"

