#include "models/Sgluons.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include "clooptools.h"

#include <cmath>

Sgluons::Sgluons(boost::property_tree::ptree const& pt) {
   mO = pt.get<double>("masses.sgluon");
}

/*
double Sgluons::sigmaSgluonsTree_uubar_OO(double alphas, double s ) {
   double a2 = Sqr(alphas);
   double b = sqrt( 1. - 4. * Sqr(mO)/s );
   return 2 * a2 * pi * Cube(b)/(9. * s);
}
*/
