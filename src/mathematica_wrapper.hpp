#ifndef MATHEMATICA_WRAPPER_HPP_
#define MATHEMATICA_WRAPPER_HPP_

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

#endif /* MATHEMATICA_WRAPPER_HPP_ */
