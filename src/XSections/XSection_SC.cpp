#include "XSections/XSection_SC.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include "cuba.h"
#include "LHAPDF/LHAPDF.h"

namespace {
int forwarder_sc(const int *, const double xx[],
   const int *, double ff[], void *userdata) {
   ff[0] = static_cast<XSection_SC*>(userdata)->integrand_sc(xx);
   return 0.;
}
int forwarder_c1(const int *, const double xx[],
   const int *, double ff[], void *userdata) {
    ff[0] = static_cast<XSection_SC*>(userdata)->integrand_c1(xx);
    return 0;
}
int forwarder_c2(const int *, const double xx[],
   const int *, double ff[], void *userdata) {
    ff[0] = static_cast<XSection_SC*>(userdata)->integrand_c2(xx);
    return 0.;
}
}

std::array<double, 3> XSection_SC::integrate() {

   // integrals dimensions, number of integrands
   static constexpr int ndim  = 3;
   static constexpr int ncomp = 1;
   // accuraccy
   const double accuracy_rel_sc {std::pow(10., -integration_precision_)};
   const double accuracy_rel_c  {std::pow(10., -integration_precision_)};
   static constexpr double accuracy_abs {1e-12};

   static constexpr int neval_min = 10'000;
   long long int neval;
   static constexpr long long int neval_max {1'000'000'000};

   int nregions, fail;

   double integral_sc[ncomp] {0.};
   double error_sc[ncomp]    {0.};
   double prob_sc[ncomp]     {0.};
   if (f_soft_.has_value()) {
      llCuhre(ndim, ncomp, forwarder_sc, this, 1,
         accuracy_rel_sc, accuracy_abs, integration_verbosity_,
         neval_min, neval_max, 1, nullptr, nullptr,
         &nregions, &neval, &fail, integral_sc, error_sc, prob_sc);
   }

   double integral_c1[ncomp] {0.};
   double error_c1[ncomp]    {0.};
   double prob_c1[ncomp]     {0.};
   if (sp_.at(0).first == SplittingKernel::Pqq || sp_.at(0).first == SplittingKernel::Pgg ||
       sp_.at(1).first == SplittingKernel::Pgg || sp_.at(1).first == SplittingKernel::Pgg) {
      llCuhre(ndim, ncomp, forwarder_c1, this, 1,
         accuracy_rel_c, accuracy_abs, integration_verbosity_,
         neval_min, neval_max, 1, nullptr, nullptr,
         &nregions, &neval, &fail, integral_c1, error_c1, prob_c1);
   }

   double integral_c2[ncomp], error_c2[ncomp], prob_c2[ncomp];
   llCuhre(ndim, ncomp, forwarder_c2, this, 1,
      accuracy_rel_c, accuracy_abs, integration_verbosity_,
      neval_min, neval_max, 1, nullptr, nullptr,
      &nregions, &neval, &fail, integral_c2, error_c2, prob_c2);

   std::array <double, 3> result_finite {
      integral_sc[0]  + integral_c1[0] + integral_c2[0],
      std::sqrt(Sqr(error_sc[0]) + Sqr(error_c1[0]) + Sqr(error_c2[0])),
      prob_sc[0]  + prob_c1[0] + prob_c2[0]
   };

  return result_finite;
}

double XSection_SC::integrand_sc(const double xx[]) {

   // integration variables
   const double x1 = 4.*Sqr(m1_/sqrtS_)    + (1-4.*Sqr(m1_/sqrtS_))    * xx[0];
   const double x2 = 4.*Sqr(m1_/sqrtS_)/x1 + (1-4.*Sqr(m1_/sqrtS_)/x1) * xx[1];

   double pdf_flux = 0.0;
   for (const auto& f : flav_) {
      pdf_flux += f.at(2) * pdf_->xfxQ(f.at(0), x1, muF_) * pdf_->xfxQ(f.at(1), x2, muF_);
   }
   pdf_flux /= x1 * x2;

   const double s12 = x1 * x2 * Sqr(sqrtS_);
   const double th = xx[2] * pi;
   const double alphas = pdf_->alphasQ(muR_);
   double result = pdf_flux * f_soft_.value()(alphas, s12, th, dS_, muR_);
   result *= to_fb;

   // jakobian
   result *= pi*Sqr(-4*Sqr(m1_) + Sqr(sqrtS_))*xx[0] /
           (Sqr(sqrtS_)*(-4*Sqr(m1_)*(-1 + xx[0]) + Sqr(sqrtS_)*xx[0]));
   return result;
}

double XSection_SC::integrand_c1(const double xx[]) {

   // integration variables
   const double x1 = 4.*Sqr(m1_)/Sqr(sqrtS_)      + (1-4.*Sqr(m1_)/Sqr(sqrtS_))      * xx[0];
   const double x2 = 4.*Sqr(m1_)/(Sqr(sqrtS_)*x1) + (1-4.*Sqr(m1_)/(Sqr(sqrtS_)*x1)) * xx[1];

   double pdf_flux = 0.0;
   for (const auto& inner : flav_) {
      pdf_flux += inner.at(2) * pdf_->xfxQ(inner.at(0), x1, muF_) * pdf_->xfxQ(inner.at(1), x2, muF_);
   }
   pdf_flux /= x1 * x2;

   const double s12 = x1 * x2 * Sqr(sqrtS_);
   const double alphas = pdf_->alphasQ(muR_);
   double result = 0.;
   for (const auto& [ker, sigma] : sp_) {
      if (!sigma) continue;
      if (ker == SplittingKernel::Pqq) {
         result += CF*(2.*std::log(dS_) + 1.5)*sigma.value()(alphas, s12);
      }
      else if (ker == SplittingKernel::Pgg) {
         static constexpr int Nf = 5;
         result += (2.*CA*std::log(dS_) + (11.*CA + 2.*Nf)/6)*sigma.value()(alphas, s12);
      }
   };
   result *= alphas/two_pi*std::log(Sqr(muR_/muF_));

   return result*(pow(-4.*pow(m1_, 2) + Sqr(sqrtS_), 2)*xx[0] /
          (Sqr(sqrtS_)*(-4*pow(m1_, 2)*(-1 + xx[0]) + Sqr(sqrtS_)*xx[0])))*pdf_flux*to_fb;
}

double XSection_SC::integrand_c2(const double xx[]) {

   // scale integration variables as Cuba works in a unit hipercube
   const double x1 = 4.*Sqr(m1_)/Sqr(sqrtS_)      + (1-4.*Sqr(m1_)/Sqr(sqrtS_))      * xx[0];
   const double x2 = 4.*Sqr(m1_)/(Sqr(sqrtS_)*x1) + (1-4.*Sqr(m1_)/(Sqr(sqrtS_)*x1)) * xx[1];
   const double z = x1 + (1-dS_-x1)*xx[2];

   if (x1/z > 1.) {
      return 0.;
   }

   double result = 0.;
   const double s12 = x1 * x2 * Sqr(sqrtS_);
   const double alphas = pdf_->alphasQ(muR_);
   for (const auto& f : flav_) {
      if (sp_.at(0).second) {
         result += f.at(2) * pdf_->xfxQ(f.at(0), x1/z, muF_)/(x1/z) * pdf_->xfxQ(f.at(1), x2, muF_)/x2
               * ( get_sp(sp_.at(0).first, z).at(0) * std::log(0.5*dC_ * s12/Sqr(muF_) * Sqr(1 - z)/z ) -
               get_sp(sp_.at(0).first, z).at(1)) * sp_.at(0).second.value()(alphas, s12);
      }
      if (sp_.at(1).second) {
         result += f.at(2) * pdf_->xfxQ(f.at(0), x2, muF_)/x2 * pdf_->xfxQ(f.at(1), x1/z, muF_)/(x1/z)
               * (get_sp(sp_.at(1).first, z).at(0) * std::log(0.5*dC_ * s12/Sqr(muF_) * Sqr(1 - z)/z ) -
               get_sp(sp_.at(1).first, z).at(1)) * sp_.at(1).second.value()(alphas, s12);
      }
   }

   result *= alphas/two_pi * 1./z * to_fb;

   // multiply by jakobian of integration variable transformation
   result *= (Sqr(-4*Sqr(m1_) + Sqr(sqrtS_))*xx[0]*(4*Sqr(m1_)*(-1 + xx[0]) - Sqr(sqrtS_)*(-1 + dS_ + xx[0])))/
        (Sqr(Sqr(sqrtS_))*(-4*Sqr(m1_)*(-1 + xx[0]) + Sqr(sqrtS_)*xx[0]));

   return result;
}
