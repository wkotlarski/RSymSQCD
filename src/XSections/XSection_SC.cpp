#include "XSections/XSection_SC.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"
#include "Li2.hpp"

#include <boost/math/special_functions/pow.hpp>

#include "cuba.h"
#include "LHAPDF/LHAPDF.h"

#include <thread>

using boost::math::pow;

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
int forwarder_c2A(const int *, const double xx[],
   const int *, double ff[], void *userdata) {
    ff[0] = static_cast<XSection_SC*>(userdata)->integrand_c2A(xx);
    return 0.;
}
int forwarder_c2B(const int *, const double xx[],
   const int *, double ff[], void *userdata) {
    ff[0] = static_cast<XSection_SC*>(userdata)->integrand_c2B(xx);
    return 0.;
}

} // anonymous

std::array<double, 3> XSection_SC::integrate() {

   // integrals dimensions, number of integrands
   static constexpr int ncomp = 1;
   // accuraccy
   const double accuracy_rel_sc {std::pow(10., -integration_precision_)};
   const double accuracy_rel_c  {std::pow(10., -integration_precision_)};
   static constexpr double accuracy_abs {1e-12};

   static constexpr int neval_min = 10'000;
   long long int neval;
   static constexpr long long int neval_max {1'000'000'000};

   int nregions, fail;

   const char* env_cubacores = std::getenv("CUBACORES");
   const int ncores = env_cubacores ? std::atoi(env_cubacores) : std::thread::hardware_concurrency();
   const int pn = 10'000; // this is Cuba's default, see arXiv:1408.6373
   cubacores(&ncores, &pn);

   double integral_sc[ncomp] {0.};
   double error_sc[ncomp]    {0.};
   double prob_sc[ncomp]     {0.};
   if (f_soft_.has_value()) {
      llCuhre(3, ncomp, forwarder_sc, this, 1,
         accuracy_rel_sc, accuracy_abs, integration_verbosity_,
         neval_min, neval_max, 1, nullptr, nullptr,
         &nregions, &neval, &fail, integral_sc, error_sc, prob_sc);
   }

   double integral_c1[ncomp] {0.};
   double error_c1[ncomp]    {0.};
   double prob_c1[ncomp]     {0.};
   if (sp_.at(0).first == SplittingKernel::Pqq || sp_.at(0).first == SplittingKernel::Pgg ||
       sp_.at(1).first == SplittingKernel::Pgg || sp_.at(1).first == SplittingKernel::Pgg) {
      llCuhre(3, ncomp, forwarder_c1, this, 1,
         accuracy_rel_c, accuracy_abs, integration_verbosity_,
         neval_min, neval_max, 1, nullptr, nullptr,
         &nregions, &neval, &fail, integral_c1, error_c1, prob_c1);
   }

   // c2A and c2B are finite in the dS->0 limit
   double integral_c2A[ncomp], error_c2A[ncomp], prob_c2A[ncomp];
   llCuhre(4, ncomp, forwarder_c2A, this, 1,
      accuracy_rel_c, accuracy_abs, integration_verbosity_,
      neval_min, neval_max, 1, nullptr, nullptr,
      &nregions, &neval, &fail, integral_c2A, error_c2A, prob_c2A);

   double integral_c2B[ncomp], error_c2B[ncomp], prob_c2B[ncomp];
   llCuhre(4, ncomp, forwarder_c2B, this, 1,
      accuracy_rel_c, accuracy_abs, integration_verbosity_,
      neval_min, neval_max, 1, nullptr, nullptr,
      &nregions, &neval, &fail, integral_c2B, error_c2B, prob_c2B);

   std::array <double, 3> result_finite {
      integral_sc[0]  + integral_c1[0] + integral_c2A[0] + integral_c2B[0],
      std::sqrt(pow<2>(error_sc[0]) + pow<2>(error_c1[0]) + pow<2>(error_c2A[0]) + pow<2>(error_c2B[0])),
      prob_sc[0]  + prob_c1[0] + prob_c2A[0] + prob_c2B[0]
   };

  return result_finite;
}

double XSection_SC::integrand_sc(const double xx[]) {

   const double x1min = 4.*pow<2>(m1_/sqrtS_);
   static constexpr double xmax = 1.;
   const double x1 = x1min + (xmax-x1min)*xx[0];
   const double x2min = 4.*pow<2>(m1_/sqrtS_)/x1;
   const double x2 = x2min + (xmax-x2min)*xx[1];

   double pdf_flux = 0.0;
   for (const auto& f : flav_) {
      pdf_flux += f.at(2) * pdf_->xfxQ(f.at(0), x1, muF_) * pdf_->xfxQ(f.at(1), x2, muF_);
   }
   pdf_flux /= x1 * x2;

   const double s12 = x1 * x2 * pow<2>(sqrtS_);
   const double th = xx[2] * pi;
   const double alphas = pdf_->alphasQ(muR_);
   double result = pdf_flux * f_soft_.value()(alphas, s12, th, dS_, muR_);
   result *= to_fb;

   // jacobian
   const double jacobian = pi*(1.-x1min)*(1.-x2min);
   result *= jacobian;
   return result;
}

double XSection_SC::integrand_c1(const double xx[]) {

   const double x1min = 4.*pow<2>(m1_/sqrtS_);
   static constexpr double xmax = 1.;
   const double x1 = x1min + (xmax-x1min)*xx[0];
   const double x2min = 4.*pow<2>(m1_/sqrtS_)/x1;
   const double x2 = x2min + (xmax-x2min)*xx[1];

   const double s12 = x1 * x2 * pow<2>(sqrtS_);

   const double Tmin = pow<2>(m1_) - 0.5*s12 - std::sqrt(0.25*pow<2>(s12) - pow<2>(m1_)*s12);
   const double Tmax = pow<2>(m1_) - 0.5*s12 + std::sqrt(0.25*pow<2>(s12) - pow<2>(m1_)*s12);
   const double T = Tmin + (Tmax-Tmin)*xx[2];

   double pdf_flux = 0.0;
   for (const auto& inner : flav_) {
      pdf_flux += inner.at(2) * pdf_->xfxQ(inner.at(0), x1, muF_) * pdf_->xfxQ(inner.at(1), x2, muF_);
   }
   pdf_flux /= x1 * x2;

   const double alphas = pdf_->alphasQ(muR_);
   const double a = s12/pow<2>(muF_)*0.5*dC_;
   double result = 0.;
   if (sp_.at(0).second) {
      const double xsec_ord0 = sp_.at(0).second.value()(alphas, s12, T, 0);
      const double xsec_ordEps = sp_.at(0).second.value()(alphas, s12, T, 1);
      if (sp_.at(0).first == SplittingKernel::Pqq) {
         result += CF*(2.*std::log(dS_) + 1.5)*(std::log(pow<2>(muR_/muF_))*xsec_ord0 - 2.*xsec_ordEps) +
            2.*CF*(
               -0.5*pow<2>(std::log(1 - dS_)) + 2*std::log(1 - dS_)*std::log(dS_) - pow<2>(std::log(dS_)) + pow<2>(std::log(1 - x1)) + (std::log(x1)*(-2*std::log(pow<2>(-1 + x1)) + std::log(x1)))/2. +
                std::log(a)*(std::log((-1 + dS_)*(-1 + x1)) - std::log(dS_*x1)) + polylogarithm::Li2(dS_) - polylogarithm::Li2(1 - x1)
			    )*xsec_ord0;
      }
      else if (sp_.at(0).first == SplittingKernel::Pgg) {
         static constexpr int Nf = 5;
         result += (2.*CA*std::log(dS_) + (11.*CA - 2.*Nf)/6)*(std::log(pow<2>(muR_/muF_))*xsec_ord0 - 2.*xsec_ordEps) +
            2.*CA*(
               -0.5*pow<2>(std::log(1 - dS_)) + 2*std::log(1 - dS_)*std::log(dS_) - pow<2>(std::log(dS_)) + pow<2>(std::log(1 - x1)) + (std::log(x1)*(-2*std::log(pow<2>(-1 + x1)) + std::log(x1)))/2. +
                std::log(a)*(std::log((-1 + dS_)*(-1 + x1)) - std::log(dS_*x1)) + polylogarithm::Li2(dS_) - polylogarithm::Li2(1 - x1)
			    )*xsec_ord0;
      }
   }
   if (sp_.at(1).second) {
      const double xsec_ord0 = sp_.at(1).second.value()(alphas, s12, T, 0);
      const double xsec_ordEps = sp_.at(1).second.value()(alphas, s12, T, 1);
      if (sp_.at(1).first == SplittingKernel::Pqq) {
         result += CF*(2.*std::log(dS_) + 1.5)*(std::log(pow<2>(muR_/muF_))*xsec_ord0 - 2.*xsec_ordEps) +
            2.*CF*(
            -0.5*pow<2>(std::log(1 - dS_)) + 2*std::log(1 - dS_)*std::log(dS_) - pow<2>(std::log(dS_)) + pow<2>(std::log(1 - x2)) + (std::log(x2)*(-2*std::log(pow<2>(-1 + x2)) + std::log(x2)))/2. +
             std::log(a)*(std::log((-1 + dS_)*(-1 + x2)) - std::log(dS_*x2)) + polylogarithm::Li2(dS_) - polylogarithm::Li2(1 - x2)
			    )*xsec_ord0;
      }
      else if (sp_.at(1).first == SplittingKernel::Pgg) {
         static constexpr int Nf = 5;
         result += (2.*CA*std::log(dS_) + (11.*CA - 2.*Nf)/6)*(std::log(pow<2>(muR_/muF_))*xsec_ord0 - 2.*xsec_ordEps) +
            2*CA*(
            -0.5*pow<2>(std::log(1 - dS_)) + 2*std::log(1 - dS_)*std::log(dS_) - pow<2>(std::log(dS_)) + pow<2>(std::log(1 - x2)) + (std::log(x2)*(-2*std::log(pow<2>(-1 + x2)) + std::log(x2)))/2. +
             std::log(a)*(std::log((-1 + dS_)*(-1 + x2)) - std::log(dS_*x2)) + polylogarithm::Li2(dS_) - polylogarithm::Li2(1 - x2)
			    )*xsec_ord0;
      }
   }

   const double jacobian = (xmax-x1min)*(xmax-x2min)*(Tmax-Tmin);
   result *= alphas/two_pi;

   return result*jacobian*pdf_flux*to_fb/(16.*pi*pow<2>(s12));
}

double XSection_SC::integrand_c2A(const double xx[]) {

   const double x1min = 4.*pow<2>(m1_/sqrtS_);
   static constexpr double xmax = 1.;
   const double x1 = x1min + (xmax-x1min)*xx[0];
   const double x2min = 4.*pow<2>(m1_/sqrtS_)/x1;
   const double x2 = x2min + (xmax-x2min)*xx[1];

   const double s12 = x1 * x2 * pow<2>(sqrtS_);

   const double z = x1 + (1-dS_-x1)*xx[2];

   const double Tmin = pow<2>(m1_) - 0.5*s12 - std::sqrt(0.25*pow<2>(s12) - pow<2>(m1_)*s12);
   const double Tmax = pow<2>(m1_) - 0.5*s12 + std::sqrt(0.25*pow<2>(s12) - pow<2>(m1_)*s12);
   const double T = Tmin + (Tmax-Tmin)*xx[3];

   const double phaseFlux= 1/(16.*pi*pow<2>(s12));
   if (x1/z > 1.) {
      return 0.;
   }

   double result = 0.;
   const double alphas = pdf_->alphasQ(muR_);
   for (const auto& f : flav_) {
      if (sp_.at(0).second) {
         double rem = 0.;
         if (sp_.at(0).first == SplittingKernel::Pqq) {
            rem = 2*CF/(1-z)*pdf_->xfxQ(f.at(0), x1, muF_)/x1;
         }
         else if (sp_.at(0).first == SplittingKernel::Pgg) {
            rem = 2*CA/(1-z)*pdf_->xfxQ(f.at(0), x1, muF_)/x1;
         }

         result += f.at(2) * pdf_->xfxQ(f.at(1), x2, muF_)/x2 * (
               (pdf_->xfxQ(f.at(0), x1/z, muF_)/(x1/z)*get_sp(sp_.at(0).first, z).at(0) -  rem) * std::log(0.5*dC_ * s12/pow<2>(muF_) * pow<2>(1 - z)/z ) -
               pdf_->xfxQ(f.at(0), x1/z, muF_)/(x1/z)*get_sp(sp_.at(0).first, z).at(1)) * sp_.at(0).second.value()(alphas, s12, T, 0)*phaseFlux;
      }
   }

   result *= alphas/two_pi * 1/z * to_fb;

   // multiply by jakobian of integration variable transformation
   const double jacobian = (1.-x1min)*(1.-x2min)*(Tmax-Tmin)*(1-dS_-x1);
   result *= jacobian;

   return result;
}

double XSection_SC::integrand_c2B(const double xx[]) {

   const double x1min = 4.*pow<2>(m1_/sqrtS_);
   static constexpr double xmax = 1.;
   const double x1 = x1min + (xmax-x1min)*xx[0];
   const double x2min = 4.*pow<2>(m1_/sqrtS_)/x1;
   const double x2 = x2min + (xmax-x2min)*xx[1];

   const double s12 = x1 * x2 * pow<2>(sqrtS_);

   const double z = x2 + (1-dS_-x2)*xx[2];

   const double Tmin = pow<2>(m1_) - 0.5*s12 - std::sqrt(0.25*pow<2>(s12) - pow<2>(m1_)*s12);
   const double Tmax = pow<2>(m1_) - 0.5*s12 + std::sqrt(0.25*pow<2>(s12) - pow<2>(m1_)*s12);
   const double T = Tmin + (Tmax-Tmin)*xx[3];

   const double phaseFlux= 1/(16.*pi*pow<2>(s12));
   if (x2/z > 1.) {
      return 0.;
   }

   double result = 0.;
   const double alphas = pdf_->alphasQ(muR_);
   for (const auto& f : flav_) {
      if (sp_.at(1).second) {
         double rem = 0.;
         if (sp_.at(1).first == SplittingKernel::Pqq) {
            rem = 2*CF/(1-z)*pdf_->xfxQ(f.at(1), x2, muF_)/x2;
         }
         else if (sp_.at(1).first == SplittingKernel::Pgg) {
            rem = 2*CA/(1-z)*pdf_->xfxQ(f.at(1), x2, muF_)/x2;
         }

         result += f.at(2) * pdf_->xfxQ(f.at(0), x1, muF_)/x1 * (
               (pdf_->xfxQ(f.at(1), x2/z, muF_)/(x2/z)*get_sp(sp_.at(1).first, z).at(0) -  rem) * std::log(0.5*dC_ * s12/pow<2>(muF_) * pow<2>(1 - z)/z ) -
               pdf_->xfxQ(f.at(1), x2/z, muF_)/(x2/z)*get_sp(sp_.at(1).first, z).at(1)) * sp_.at(1).second.value()(alphas, s12, T, 0)*phaseFlux;
      }
   }

   result *= alphas/two_pi * 1/z * to_fb;

   // multiply by jakobian of integration variable transformation
   const double jacobian = (1.-x1min)*(1.-x2min)*(Tmax-Tmin)*(1-dS_-x2);
   result *= jacobian;

   return result;
}
