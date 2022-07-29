#include "XSection_SC.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include "LHAPDF/LHAPDF.h"

namespace {
int forwarder_sc(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {
    return static_cast<XSection_SC*>(userdata)->integrand_sc(ndim, xx, ncomp, ff, nullptr);
}
int forwarder_c1(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {
    return static_cast<XSection_SC*>(userdata)->integrand_c1(ndim, xx, ncomp, ff, nullptr);
}
int forwarder_c2(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {
    return static_cast<XSection_SC*>(userdata)->integrand_c2(ndim, xx, ncomp, ff, nullptr);
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

   cubareal integral_sc[ncomp], error_sc[ncomp], prob_sc[ncomp];
   llCuhre(ndim, ncomp, forwarder_sc, this, 1,
      accuracy_rel_sc, accuracy_abs, integration_verbosity_,
      neval_min, neval_max, 1, NULL, NULL,
      &nregions, &neval, &fail, integral_sc, error_sc, prob_sc);

   double integral_c1[ncomp] = {0.};
   double error_c1[ncomp]    = {0.};
   double prob_c1[ncomp]     = {0.};
   if (sp_.at(0).first == SplittingKernel::Pqq || sp_.at(0).first == SplittingKernel::Pgg ||
       sp_.at(1).first == SplittingKernel::Pgg || sp_.at(1).first == SplittingKernel::Pgg) {
      llCuhre(ndim, ncomp, forwarder_c1, this, 1,
         accuracy_rel_c, accuracy_abs, integration_verbosity_,
         neval_min, neval_max, 1, NULL, NULL,
         &nregions, &neval, &fail, integral_c1, error_c1, prob_c1);
   }

   cubareal integral_c2[ncomp], error_c2[ncomp], prob_c2[ncomp];
   llCuhre(ndim, ncomp, forwarder_c2, this, 1,
      accuracy_rel_c, accuracy_abs, integration_verbosity_,
      neval_min, neval_max, 1, NULL, NULL,
      &nregions, &neval, &fail, integral_c2, error_c2, prob_c2);

   std::array <double, 3> result_finite {
      integral_sc[0]  + integral_c1[0] + integral_c2[0],
      std::sqrt(Sqr(error_sc[0]) + Sqr(error_c1[0]) + Sqr(error_c2[0])),
      prob_sc[0]  + prob_c1[0] + prob_c2[0]
   };

  return result_finite;
}

int XSection_SC::integrand_sc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

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
   ff[0] = pdf_flux * f_soft_(alphas, s12, th, dS_, muR_);
   ff[0] *= to_fb;

   // jakobian
   ff[0] *= pi*Sqr(-4*Sqr(m1_) + Sqr(sqrtS_))*xx[0] /
           (Sqr(sqrtS_)*(-4*Sqr(m1_)*(-1 + xx[0]) + Sqr(sqrtS_)*xx[0]));
   return 0;
}

// @todo if mu_r != mu_f one needs these term,
//    otherwise it's 0

int XSection_SC::integrand_c1(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

   // integration variables
   const double x1 = 4.*Sqr(m1_)/Sqr(sqrtS_)      + (1-4.*Sqr(m1_)/Sqr(sqrtS_))      * xx[0];
   const double x2 = 4.*Sqr(m1_)/(Sqr(sqrtS_)*x1) + (1-4.*Sqr(m1_)/(Sqr(sqrtS_)*x1)) * xx[1];
   const double th = pi * xx[2];

   double pdf_flux = 0.0;
   for (const auto& inner : flav_) {
      pdf_flux += inner.at(2) * pdf_->xfxQ(inner.at(0), x1, muF_) * pdf_->xfxQ(inner.at(1), x2, muF_);
   }
   pdf_flux /= x1 * x2;

   const double s12 = x1 * x2 * Sqr(sqrtS_);
   const double alphas = pdf_->alphasQ(muR_);
   double result = 0.;
   for (const auto& el : sp_) {
      if (el.first == SplittingKernel::Pqq) {
         result += CF*(2.*std::log(dS_) + 1.5)*el.second(alphas, s12);
      }
      else if (el.first == SplittingKernel::Pgg) {
         static constexpr int Nf = 5;
         result += (2.*CA*std::log(dS_) + (11.*CA + 2.*Nf)/6)*el.second(alphas, s12);
      }
   };
   result *= alphas/two_pi*std::log(Sqr(muR_/muF_));

   ff[0] = result*(pow(-4.*pow(m1_, 2) + Sqr(sqrtS_), 2)*xx[0] /
          (Sqr(sqrtS_)*(-4*pow(m1_, 2)*(-1 + xx[0]) + Sqr(sqrtS_)*xx[0])))*pdf_flux*to_fb;

   return 0;
}

int XSection_SC::integrand_c2(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

   // scale integration variables as Cuba works in a unit hipercube
   const double x1 = 4.*Sqr(m1_)/Sqr(sqrtS_)      + (1-4.*Sqr(m1_)/Sqr(sqrtS_))      * xx[0];
	const double x2 = 4.*Sqr(m1_)/(Sqr(sqrtS_)*x1) + (1-4.*Sqr(m1_)/(Sqr(sqrtS_)*x1)) * xx[1];
   const double z = x1 + (1-dS_-x1)*xx[2];

   ff[0] = 0.0;

   if (x1/z > 1.) {
      return 0;
   }

   const double s12 = x1 * x2 * Sqr(sqrtS_);
   const double alphas = pdf_->alphasQ(muR_);
   for (const auto& f : flav_) {
      ff[0] += f.at(2) * pdf_->xfxQ(f.at(0), x1/z, muF_)/(x1/z) * pdf_->xfxQ(f.at(1), x2, muF_)/x2
            * ( get_sp(sp_.at(0).first, z).at(0) * std::log( dC_/2. * s12/Sqr(muF_) * Sqr(1 - z)/z ) -
           get_sp(sp_.at(0).first, z).at(1)) * sp_.at(0).second(alphas, s12);
      ff[0] += f.at(2) * pdf_->xfxQ(f.at(0), x2, muF_)/x2 * pdf_->xfxQ(f.at(1), x1/z, muF_)/(x1/z)
           * (get_sp(sp_.at(1).first, z).at(0) * std::log( dC_/2. * s12/Sqr(muF_) * Sqr(1 - z)/z ) -
           get_sp(sp_.at(1).first, z).at(1)) * sp_.at(1).second(alphas, s12);
   }

   ff[0] *= alphas/two_pi * 1./z * to_fb;

   // multiply by jakobian of integration variable transformation
   ff[0] *= (Sqr(-4*Sqr(m1_) + Sqr(sqrtS_))*xx[0]*(4*Sqr(m1_)*(-1 + xx[0]) - Sqr(sqrtS_)*(-1 + dS_ + xx[0])))/
        (Sqr(Sqr(sqrtS_))*(-4*Sqr(m1_)*(-1 + xx[0]) + Sqr(sqrtS_)*xx[0]));

   return 0;
}
