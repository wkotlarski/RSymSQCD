#include "XSection_Tree.hpp"

#include "cuba.h"

#include <functional>

namespace {
int forwarder(const int *ndim, const double xx[],
   const int *ncomp, double ff[], void *userdata) {
    return static_cast<XSection_Tree*>(userdata)->integrand(ndim, xx, ncomp, ff, nullptr);
}
}
/*
 * E1 = sqrtS/2 (1 + (m1^2 - m2^2)/s)
 * E2 = sqrtS/2 (1 + (m2^2 - m1^2)/s)
 * b  = sqrt(1 - 2(m1^2 + m2^2)/s + ((m1 + m2)(m1-m2)/s)^2)
 * p  = sqrtS/2 b
 *
 */
int XSection_Tree::integrand(const int *ndim, const double xx[],
    const int *ncomp, double ff[], void *userdata) {

   const double x1min = 4.*Sqr(m1)/S;
   static constexpr double xmax = 1.;
   const double x1 = x1min + (xmax - x1min ) * xx[0];
   const double x2min = 4.*Sqr(m1)/(S*x1);
   const double x2 = x2min + (xmax - x2min) * xx[1];
   const double s = S * x1 * x2;     //partonic

   double pdf_flux = 0.0;
   for (const auto& inner : flav_) {
      pdf_flux += inner.at(2) * pdf->xfxQ(inner.at(0), x1, mu_f) * pdf->xfxQ(inner.at(1), x2, mu_f);
   }
   pdf_flux /= (x1 * x2);

   
   const double alphas = pdf->alphasQ(muR);
   /* integration of |M^B|^2 */
   if (!processID.partonic) {
      const double Tmin = Sqr(m1) - 0.5*s - std::sqrt(0.25*Sqr(s) - Sqr(m1)*s);
      const double Tmax = Sqr(m1) - 0.5*s + std::sqrt(0.25*Sqr(s) - Sqr(m1)*s);
      const double T = Tmin + (Tmax-Tmin)*xx[2];
      const double jacobian = (Tmax-Tmin)*(xmax-x1min)*(xmax-x2min);
      const double squaredM = f(alphas, s, T);
      double dSigmaPart = squaredM/(16.*pi*Sqr(s));

      ff[0] = dSigmaPart * pdf_flux * jacobian * to_fb;
   }
   /* integration of partonic cross section */
   else {
      ff[0] = (processID.*processID.sigmaPartTree1)(alphas, s) * to_fb * pdf_flux *
         pow(-4.*pow(m1, 2) + S, 2)*xx[0] /
         (S*(-4*pow(m1, 2)*(-1 + xx[0]) + S*xx[0]));
   }

   return 0;
}

 std::array<double, 3> XSection_Tree::integrate() {
   constexpr int ndim = 3;
   constexpr int ncomp = 1;
   constexpr int nvec = 1;
   constexpr double accuracy_rel = 1e-6;
   constexpr double accuracy_abs = 1e-12;
   constexpr int eval_min = 1000;
   constexpr int eval_max = 1000000000;
   const int verbose = vm["verbosity-born"].as<int>();        // adjust shown output 0 ... 3
   constexpr int last = 4;
   constexpr int key = 0;
   int nregions, neval, fail;
   double integral[ncomp], error[ncomp], prob[ncomp];

   Cuhre( ndim, ncomp, forwarder, this, nvec,
      accuracy_rel, accuracy_abs, verbose | last,
      eval_min, eval_max, key, NULL, NULL,
      &nregions, &neval, &fail, integral, error, prob );

   std::array <double, 3> result{ integral[0], error[0], prob[0] }; 

   return result;
}
