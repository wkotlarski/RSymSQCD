#include "models/MRSSM.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include <Eigen/Dense>
#include "clooptools.h"

#include <cmath>
#include <iostream>

MRSSM::MRSSM(boost::property_tree::ptree const& pt) {

   MassTop = pt.get<double>("masses.top");
   MassGlu = pt.get<double>("masses.gluino");
   MasssigmaO = pt.get<double>("masses.pseudoscalar_sgluon");
   MassphiO  = sqrt( pow(MasssigmaO,2) + 4.0 * pow(MassGlu, 2) );
   MassSq = pt.get<double>("masses.squarks");
   double eta_sign = pt.get<double>("technical parameters.eta_sign");
   double delta = pt.get<double>("technical parameters.delta");
   WidthGlu = pt.get<double>("technical parameters.WidthOverMass") * MassGlu;

   // choose a gage vector \eta for DR matrix elements
   eta = {std::sqrt(1.+Sqr(delta)), 0., delta, eta_sign};
}

#include "mrssm_born_mes.cpp"
#include "mrssm_born_xsecs.cpp"

#include "mrssm_uu_suLsuR_virt_matrix.cpp"
#include "mrssm_ud_suLsdR_virt_matrix.cpp"

#include "mrssm_uu_suLsuRg_soft.cpp"
#include "mrssm_uu_suLsuRg_hard.cpp"
#include "mrssm_gu_suLsuRubar_hard.cpp"
#include "mrssm_gu_suLsuRubar_hard-DR.cpp"
#include "mrssm_gu_suLsuRubar_hard-DR_wEtaDep.cpp"
#include "mrssm_gu_suLsuRubar_hard-DS.cpp"

#include "mrssm_gg_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_ddbar_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_uubar_suLsuLdagger_virt_matrix.cpp"
#include "mrssm_uubar_suLsuLdaggerg_soft.cpp"
#include "mrssm_uubar_suLsuLdaggerg_hard.cpp"
#include "mrssm_ddbar_suLsuLdaggerg_soft.cpp"
#include "mrssm_ddbar_suLsuLdaggerg_hard.cpp"
#include "mrssm_gd_suLsuLdaggerd_hard.cpp"
// wEta_noSimplify is faster than wEta
#include "mrssm_gu_suLsuLdaggeru_hard-DS.cpp"
#include "mrssm_gu_suLsuLdaggeru_hard-DR.cpp"
#include "mrssm_gu_suLsuLdaggeru_hard-DR_wEta_noSimplify.cpp"
#include "mrssm_gu_suLsuLdaggeru_hard.cpp"
#include "mrssm_gg_suLsuLdaggerg_hard.cpp"
#include "mrssm_gg_suLsuLdaggerg_soft.cpp"

