#include "XSections/XSection_Virt.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include "clooptools.h"
#include "cuba.h"
#include "LHAPDF/LHAPDF.h"

namespace {

int forwarder(const int *, const double xx[],
   const int *, double ff[], void *userdata) {
   ff[0] = static_cast<XSection_Virt*>(userdata)->integrand(xx);
   return 0.;
}

}

double XSection_Virt::integrand(const double xx[]) {
   const double x1min = 4.*Sqr(m1_/sqrtS_);
   static constexpr double xmax = 1.;
   const double x1 = x1min + (xmax-x1min)*xx[1];
   const double x2min = 4.*Sqr(m1_/sqrtS_)/x1;
   const double x2 = x2min + (xmax-x2min)*xx[2];
   const double s = Sqr(sqrtS_)*x1*x2;
   const double Tmin = Sqr(m1_) - 0.5*s - std::sqrt(0.25*Sqr(s) - Sqr(m1_)*s);
   const double Tmax = Sqr(m1_) - 0.5*s + std::sqrt(0.25*Sqr(s) - Sqr(m1_)*s);
   const double T = Tmin + (Tmax-Tmin)*xx[0];
   const double jacobian = (Tmax-Tmin)*(1.-x1min)*(1.-x2min);

   const double alphas = pdf_->alphasQ(muR_);

   int FiniteGs = 1;
   int Divergence = 0;
   double Dminus4Pow = 0;
   const double dSigmaPart1 = f(alphas, s, T, FiniteGs, Dminus4Pow, Divergence, muR_);

   FiniteGs = 0;

   // contraction with O(eps) from Dminus4
   Divergence = -1;
   Dminus4Pow = 1;
   const double dSigmaPart2 = -2.*f(alphas, s, T, FiniteGs, Dminus4Pow, Divergence, muR_);

   // contraction with O(eps^2) from O(eps^2) mismatch between soft and virt prefactors
   Divergence = -2;
   Dminus4Pow = 0;
   double dSigmaPart3 = pi_sqr_div_six * f(alphas, s, T, FiniteGs, Dminus4Pow, Divergence, muR_);

   // contraction with O(eps^2) from Dminus4
   Divergence = -2;
   Dminus4Pow = 2;
   const double dSigmaPart4 = 4.*f(alphas, s, T, FiniteGs, Dminus4Pow, Divergence, muR_);

   const double dSigmaHad = (dSigmaPart1 + dSigmaPart2 + dSigmaPart3 + dSigmaPart4)/(16.*pi*Sqr(s));

   double pdf_flux = 0.0;
   for (const auto& fv : flav_) {
      pdf_flux += fv.at(2) * pdf_->xfxQ(fv.at(0), x1, muF_) * pdf_->xfxQ(fv.at(1), x2, muF_);
   }
   pdf_flux /= x1*x2;

   // return result in femtobarns
   return dSigmaHad * jacobian * pdf_flux * to_fb;
}

std::array<double, 3> XSection_Virt::integrate() {
   static constexpr int ndim = 3;
   static constexpr int ncomp = 1;
   static constexpr int nvec = 1;
   static constexpr double accuracy_abs = 1e-12;
   static constexpr int eval_min = 1000;
   static constexpr int eval_max = 1000000;
   int nregions, neval, fail;
   double integral[ncomp], error[ncomp], prob[ncomp];
   static constexpr int last = 4;
   static constexpr int key = 0;
   // Divonne specific
   static constexpr int seed = 1;
   static constexpr double border = 1e-6; // don't go lower than 1e-6 for 1e-4 relative accuraccy
   static constexpr int key1 = 47;
   static constexpr int key2 = 1;
   static constexpr int key3 = 1;
   static constexpr int maxpass = 5;
   static constexpr double maxchisq = 10.;
   static constexpr double mindeviation = 0.25;
   static constexpr int ngiven = 0;
   static constexpr int ldxgiven = ndim;
   static constexpr int nextra = 0;
   const double prec_virt = std::pow(10., -integration_precision_);

   static bool looptools_initialized = false;
   if(!looptools_initialized) {
      std::cout << std::endl;
      ltini();
      setuvdiv(1);
      looptools_initialized = true;
   }
   /*
    *    LoopTools fails for phase space points near the border.
    *    That's why we have to use Divonne with option border.
    *    It's not clear if this is just a sign of numerical instability
    *    of LoopTools or something more serious.
    */
   Divonne(ndim, ncomp, forwarder, this, nvec,
        prec_virt, accuracy_abs, integration_verbosity_, seed,
        eval_min, eval_max, key1, key2, key3, maxpass,
        border, maxchisq, mindeviation,
        ngiven, ldxgiven, nullptr, nextra, nullptr,
        nullptr, nullptr,
        &nregions, &neval, &fail, integral, error, prob);

   return {integral[0], error[0], prob[0]};
}
