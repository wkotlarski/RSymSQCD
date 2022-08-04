#include "models/Sgluons.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include <boost/math/special_functions/pow.hpp>
#include "clooptools.h"

#include <cmath>

using boost::math::pow;

Sgluons::Sgluons(SgluonParameters const& params)
   : mO {params.mO}, mt {params.mt}
   {};

// Eq. 7.2 of arXiv:1611.06622
double Sgluons::sigmaSgluonsTree_qqbar_OO(double alphas, double s) {
   const double b = std::sqrt(1. - 4.*pow<2>(mO)/s);
   return 2*pow<2>(alphas)*pi*pow<3>(b)/(9.*s);
}

// Eq. 7.3 of arXiv:1611.06622
double Sgluons::sigmaSgluonsTree_gg_OO(double alphas, double s) {
   const double b = std::sqrt(1. - 4.*pow<2>(mO)/s);
   return 3*pow<2>(alphas)*pi/(32.*s)*(27*b - 17*pow<3>(b) + 6*(-3 + 2*pow<2>(b) + pow<4>(b))*std::atanh(b));
}

#include "sgluons_gg_OOg_soft.cpp"
#include "sgluons_gg_OOg_hard.cpp"
#include "sgluons_qqbar_OOg_soft.cpp"
#include "sgluons_qqbar_OOg_hard.cpp"
