#include "models/MRSSM.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include <boost/math/special_functions/pow.hpp>

#include <Eigen/Dense>
#include "clooptools.h"

#include <cmath>

using boost::math::pow;

MRSSM::MRSSM(MRSSMParameters const& params)
   : MassTop {params.MassTop},
     MassGlu {params.MassGlu},
     MasssigmaO {params.MasssigmaO},
     MassphiO {std::sqrt(Sqr(params.MasssigmaO) + 4*Sqr(params.MassGlu))},
     MassSq {params.MassSq},
     eta {{std::sqrt(1.+Sqr(params.delta)), 0., params.delta, params.eta_sign}},
     WidthGlu {params.WidthGlu}
   {};

#include "mrssm_born_mes.cpp"
#include "mrssm_born_xsecs.cpp"

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
