#include "models/Sgluons.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include "clooptools.h"

#include <cmath>

Sgluons::Sgluons(SgluonParameters const& params)
   : mO {params.mO}
   {};

double Sgluons::sigmaSgluonsTree_qqbar_OO(double alphas, double s) {
   double a2 = Sqr(alphas);
   double b = sqrt( 1. - 4. * Sqr(mO)/s );
   return 2 * a2 * pi * Power3(b)/(9. * s);
}

#include "sgluons_gg_OOg_soft.cpp"
#include "sgluons_gg_OOg_hard.cpp"
#include "sgluons_qqbar_OOg_soft.cpp"
#include "sgluons_qqbar_OOg_hard.cpp"
