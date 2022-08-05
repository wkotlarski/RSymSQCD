#include "models/MSSM.hpp"

#include <boost/math/special_functions/pow.hpp>

#include "clooptools.h"

using boost::math::pow;

double MSSM::matrixMSSMTree_uubar_suLsuRdagger(double alphas, double S, double T) const { // agrees with Philip
   const double U = 2*MassSq*MassSq - S - T;
   return (-631.6546816697189*pow<2>(alphas)*(pow(MassSq,4) - 1.*T*U))/(S*S);
}

#include "mssm_uu_suLsuR_virt_matrix.cpp"
