#include "models/MRSSM.hpp"
#include "constants.hpp"

#include <boost/math/special_functions/pow.hpp>

#include <Eigen/Dense>
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

double MRSSM::matrixTree_GG_suLsuLdagger(double alphas, double S, double T) const { // agrees with Philip
   const double U = 2*pow<2>(MassSq) - S - T;
   static constexpr double k = 2.*2*8*8;
   static constexpr double h = 1.;
   return  h/k*(-157.91367041742973*pow<2>(alphas)*((96.*(pow<4>(MassSq) + T*U - (pow<2>(MassSq))*(T + U)))/(pow<2>(S)) + pow<2>(MassSq)*(37.333333333333336/(-(pow<2>(MassSq)) + U) - 85.33333333333333*(T/pow<2>(-(pow<2>(MassSq)) + T) + U/pow<2>(-(pow<2>(MassSq)) + U))) + (37.333333333333336*(pow<2>(MassSq)) + (5.333333333333333*(9.*pow<4>(MassSq) - 3.*(pow<2>(MassSq))*(2.*(pow<2>(MassSq)) + S) + (S + T)*(S + U)))/(-(pow<2>(MassSq)) + U))/(-(pow<2>(MassSq)) + T) - (48.*((5.*pow<4>(MassSq) + pow<2>(MassSq)*U - T*(-(pow<2>(MassSq)) + U) - (pow<2>(MassSq))*(3.*S + 4.*T + 2.*U))/(-(pow<2>(MassSq)) + U) + (5.*pow<4>(MassSq) + pow<2>(MassSq)*U - T*(-(pow<2>(MassSq)) + U) - (pow<2>(MassSq))*(3.*S + 2.*T + 4.*U))/(-(pow<2>(MassSq)) + T)))/S));
}

double MRSSM::matrixTree_ddbar_suLsuLdagger(double alphas, double S, double T ) const { // agrees with Philip
   const double U = 2*pow<2>(MassSq) - S - T;
   static constexpr double k = 2.*2*3*3;
   static constexpr double h = 2.*2;
   return  h/k*(-631.6546816697189*pow<2>(alphas)*(pow<4>(MassSq) - T*U))/pow<2>(S);
}

double MRSSM::matrixTree_uubar_suLsuLdagger(double alphas, double S, double T) const {
   const double U = 2*pow<2>(MassSq) - S - T;
   static constexpr double k = 2.*2*3*3;
   static constexpr double h = 2.*2;
   return  h/k*((315.82734083485946*pow<2>(alphas)*(-pow<4>(MassSq) + T*U))/(pow<2>(S)) + 355.3057584392169*pow<2>(alphas)*pow(0.3333333333333333/S - 1./(-(pow<2>(MassGlu)) + T),2)*(-pow<4>(MassSq) + T*U) - 236.8705056261446*pow<2>(alphas)*(0.3333333333333333/S - 1./(-(pow<2>(MassGlu)) + T))*(1/S - 0.3333333333333333/(-(pow<2>(MassGlu)) + T))*(-pow<4>(MassSq) + T*U) + 355.3057584392169*pow<2>(alphas)*pow(1/S - 0.3333333333333333/(-(pow<2>(MassGlu)) + T),2)*(-pow<4>(MassSq) + T*U));
}

double MRSSM::matrixTree_udbar_suLsdLdagger(double alphas, double S, double T) const {
   double U = 2*pow<2>(MassSq) - S - T;
   return  (315.82734083485946*pow<2>(alphas)*(-pow<4>(MassSq) + T*U))/pow<2>(-(pow<2>(MassGlu)) + T);
}

double MRSSM::matrixTree_uu_suLsuR(double alphas, double S, double T) const {
   static constexpr double h = 2.*2;
   static constexpr double k = 2.*2*3*3;
   const double U = 2*pow<2>(MassSq) - S - T;
   return h/k*((315.82734083485946*pow<2>(alphas)*(-pow<4>(MassSq) + T*U))/pow<2>(-pow<2>(MassGlu) + T) +
           (315.82734083485946*pow<2>(alphas)*(-pow<4>(MassSq) + T*U))/pow<2>(-pow<2>(MassGlu) + U));
}

double MRSSM::matrixTree_uubar_glglbar(double alphaS, double S, double T) const
{
   const double U = 2*pow<2>(MassGlu) - S - T;
   const double MsquaredReal = 0.5*16.*52.63789013914324*(alphaS*alphaS)*(36*pow<2>(MassGlu)/S + 18*pow<2>(pow<2>(MassGlu) - T)/pow<2>(S) + (4*pow<2>(pow<2>(MassGlu) - T))/pow<2>(pow<2>(MassSq) - T) + (9.*(pow<2>(MassGlu)))/(-(pow<2>(MassSq)) + T) + (9.*pow<2>(pow<2>(MassGlu) - T))/(S*(-(pow<2>(MassSq)) + T)) + (18.*pow<2>(pow<2>(MassGlu) - U))/(pow<2>(S)) + (4.*pow<2>(pow<2>(MassGlu) - U))/pow<2>(pow<2>(MassSq) - U) + (9.*(pow<2>(MassGlu)))/(-(pow<2>(MassSq)) + U) + (9.*pow<2>(pow<2>(MassGlu) - U))/(S*(-(pow<2>(MassSq)) + U)))
;
   return MsquaredReal/18.;
}

double MRSSM::matrixTree_gg_glglbar(double alphaS, double S, double T) const
{
   const double U = 2*pow<2>(MassGlu) - S - T;
   const double MsquaredReal = 8.*pow<2>(4.*pi*alphaS)*3.*24.*(1.-(U - pow<2>(MassGlu))*(T - pow<2>(MassGlu))/(pow<2>(S)))*(pow<2>(S)/((T - pow<2>(MassGlu))*(U - pow<2>(MassGlu))) - 2. + 4.*pow<2>(MassGlu)*S/((T - pow<2>(MassGlu))*(U - pow<2>(MassGlu)))*(1. - pow<2>(MassGlu)*S/((T - pow<2>(MassGlu))*(U - pow<2>(MassGlu)))));

return MsquaredReal/256.;
}

double MRSSM::sigmaTree_uubar_suLsuLdagger(double alphas, double s12) const {
   const double Alfas2 = pow<2>(alphas);
   return (Alfas2*pi*(8*(-3*pow<2>(MassSq) + MassGlu*MassGlu - 2*s12)*sqrt(s12*(-4*pow<2>(MassSq) + s12)) +
       8*(2*pow<4>(MassSq) + 2*pow<4>(MassGlu) - 4*pow<2>(MassGlu)*s12 - 3*(s12*s12) + pow<2>(MassSq)*(-4*pow<2>(MassGlu) + 6*s12))*
        atanh(sqrt(s12*(-4*pow<2>(MassSq) + s12))/(2*pow<2>(MassSq) - 2*pow<2>(MassGlu) - s12))))/(54.*pow<3>(s12));
}

// checked against MG with SUSYQCD model
double MRSSM::sigmaTree_ddbar_suLsuLdagger(double alphas, double s12 ) const {
   const double Alfas2 = pow<2>(alphas);
   return (2*Alfas2*pi*pow(-4*pow<2>(MassSq) + s12,1.5))/(27.*pow(s12,2.5));
}

double MRSSM::sigmaTree_gg_suLsuLdagger(double alphas, double s12 ) const {
   const double Alfas2 = pow<2>(alphas);
   return (Alfas2*pi*(sqrt(s12*(-4*pow<2>(MassSq) + s12))*(62*pow<2>(MassSq) + 5*s12) - 16*pow<2>(MassSq)*(pow<2>(MassSq) + 4*s12)*atanh(sqrt(1 - (4*pow<2>(MassSq))/s12))))/(48.*pow<3>(s12));
}

/*
 *    checked with MadGraph
 */
double MRSSM::sigmaTree_uu_suLsuR(double alphas, double s) const {
   const double MGl2 = pow<2>(MassGlu);
   const double a = pow<2>(alphas);
   return (-4.*a*pi*(2.*sqrt(s*(-4*pow<2>(MassSq) + s)) + (2.*pow<2>(MassSq) - 2.*MGl2 - s)*
      log((4.*MGl2 + pow<2>(1. + sqrt(1. - (4.*pow<2>(MassSq))/s))*s)/
      (4.*MGl2 + pow<2>(-1. + sqrt(1. - (4.*pow<2>(MassSq))/s))*s))))/(9.*pow<2>(s));
}

#include "mrssm_uu_suLsuR_virt_matrix.cpp"
#include "mrssm_ud_suLsdR_virt_matrix.cpp"

#include "mrssm_uu_suLsuRg_soft.cpp"
#include "mrssm_uu_suLsuRg_hard.cpp"
#include "mrssm_gu_suLsuRubar_hard.cpp"
#include "mrssm_gu_suLsuRubar_hard-DR.cpp"
#include "mrssm_gu_suLsuRubar_hard-DS.cpp"
#include "mrssm_gu_suLsuRubar_hard-DS_unsympli.cpp"

#include "mrssm_gg_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_ddbar_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_uubar_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_uubar_suLsuLdaggerg_soft.cpp"
#include "mrssm_uubar_suLsuLdaggerg_hard.cpp"
#include "mrssm_ddbar_suLsuLdaggerg_soft.cpp"
#include "mrssm_ddbar_suLsuLdaggerg_hard.cpp"
#include "mrssm_gd_suLsuLdaggerd_hard.cpp"
#include "mrssm_gu_suLsuLdaggeru_hard-DS.cpp"
#include "mrssm_gu_suLsuLdaggeru_hard-DR.cpp"
#include "mrssm_gu_suLsuLdaggeru_hard.cpp"
#include "mrssm_gg_suLsuLdaggerg_hard.cpp"
#include "mrssm_gg_suLsuLdaggerg_soft.cpp"

//#include "mrssm_gu_suLsuRubar_hard-DR_wEtaDep.cpp"
// #include "mrssm_gu_suLsuLdaggeru_hard-DR_wEta_noSimplify.cpp"
// wEta_noSimplify is faster than wEta
