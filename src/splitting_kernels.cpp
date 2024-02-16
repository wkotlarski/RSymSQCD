#include "constants.hpp"
#include "mathematica_wrapper.hpp"
#include "splitting_kernels.hpp"

std::array<double, 2> Pqq(double z) {
   return {CF*(1.+z*z)/(1.-z), -CF*(1.-z)};
}

std::array<double, 2> Pgq(double z) {
   return Pqq(1 - z);
};

std::array<double, 2> Pgg(double z) {
   return {2.*CA*(z/(1-z) + (1.-z)/z + z*(1.-z)), 0.};
}

std::array<double, 2> Pqg(double z) {
   return {0.5*(Sqr(z) + Sqr(1.-z)), -z*(1.-z)};
}

std::array<double, 2> get_sp(SplittingKernel sp, double z) {
   switch(sp) {
      case SplittingKernel::Pqq:
         return Pqq(z);
      case SplittingKernel::Pgg:
         return Pgg(z);
      case SplittingKernel::Pqg:
         return Pqg(z);
      case SplittingKernel::Pgq:
         return Pgq(z);
      default:
         throw("Unknown splitting kernel");
   }
}