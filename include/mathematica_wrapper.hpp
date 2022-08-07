#ifndef MATHEMATICA_WRAPPER_HPP_
#define MATHEMATICA_WRAPPER_HPP_

#include <complex>

template <typename Arg>
inline double Log (Arg arg) {
   return std::log(arg);
}

template <typename Arg>
inline double Cos (Arg arg) {
   return std::cos(arg);
}

template <typename Base, typename Exponent>
Base Power(Base base, Exponent exp) noexcept
{
   return std::pow(base, exp);
}

template <typename T>
constexpr std::complex<T> Sqr(const std::complex<T>& a) noexcept
{
   return a * a;
}

template <typename T, class = std::enable_if_t<std::is_arithmetic<T>::value,T>>
constexpr T Sqr(T a) noexcept
{
   return a * a;
}

template <typename Base>
constexpr Base Power2(Base b) noexcept
{
   return b * b;
}

template <typename Base>
constexpr Base Power3(Base b) noexcept
{
   return b * b * b;
}

template <typename Base>
constexpr Base Power4(Base b) noexcept
{
   return Power2(Power2(b));
}

#endif // MATHEMATICA_WRAPPER_HPP_
