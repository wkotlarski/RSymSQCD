#ifndef MATHEMATICA_WRAPPER_HPP_
#define MATHEMATICA_WRAPPER_HPP_

#include "dilog.hpp"

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

template <typename Arg> inline double Log (Arg arg) {
    return log(arg);
}

template <typename Arg> inline double Cos (Arg arg) {
    return cos(arg);
}

template <typename Arg> inline double Sin (Arg arg) {
    return sin(arg);
}

inline std::complex<double> PolyLog(int i, double arg) {
   if (i != 2) {
      throw std::invalid_argument("We handle only dilogarithms");
   }
   std::complex<double> z(arg, 0.);
   return dilogarithm::dilog(z);
}

#endif // MATHEMATICA_WRAPPER_HPP_
