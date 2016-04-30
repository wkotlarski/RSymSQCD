/*
 * mathematica_wrapper.hpp
 *
 *  Created on: 6 mar 2016
 */

#ifndef MATHEMATICA_WRAPPER_HPP_
#define MATHEMATICA_WRAPPER_HPP_

#include "gsl/gsl_sf_lambert.h"

template <typename Base, typename Exponent>
  inline Base Power(Base base, Exponent exp) {
     return std::pow(base, exp);
  }

inline double Sqrt(double a) {
   return std::sqrt(a);
}

inline double Abs(double x) {
   return std::abs(x);
}

inline double ProductLog( int i, double x) {
  return gsl_sf_lambert_Wm1( x );
}




#endif /* MATHEMATICA_WRAPPER_HPP_ */
