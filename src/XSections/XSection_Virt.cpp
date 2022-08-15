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
    const double s = Sqr(sqrtS_)*x1*x2;     //partonic
    const double Tmin = Sqr(m1_) - 0.5*s - std::sqrt(0.25*Sqr(s) - Sqr(m1_)*s);
    const double Tmax = Sqr(m1_) - 0.5*s + std::sqrt(0.25*Sqr(s) - Sqr(m1_)*s);
    const double T = Tmin + (Tmax-Tmin)*xx[0];
    const double jacobian = (Tmax-Tmin)*(1.-x1min)*(1.-x2min);

    int FiniteGs = 1;
    double Dminus4 = 0;
    int Divergence = 0;     // O(eps)
    const double alphas = pdf_->alphasQ(muR_);

    double squaredMReal = f(alphas, s, T, FiniteGs, Dminus4, Divergence, muR_);

    double dSigmaPart1 = squaredMReal/(8.*pi*Sqr(s));

    // contraction with O(eps) from Dminus4
    Divergence = -1;           // O(eps)
    FiniteGs = 0;
    squaredMReal = f(alphas, s, T, FiniteGs, Dminus4, Divergence, muR_);

    Dminus4 = -2.;
    const double squaredMRealMinus2 = f(alphas, s, T, FiniteGs, Dminus4, Divergence, muR_);

    const double dSigmaPart3 = (squaredMRealMinus2 - squaredMReal)/(8.*pi*Sqr(s));

    // contraction with O(eps^2) prefactor of loop integral
    // and with product of O(eps) prefactors of phase space and loop integral
    Divergence = -2;
    Dminus4 = 0;
    squaredMReal = f(alphas, s, T, FiniteGs, Dminus4, Divergence, muR_);

   double dSigmaPart4 = squaredMReal*pi/(48.*Sqr(s));

   const double dSigmaHad = (dSigmaPart1 + dSigmaPart3 + dSigmaPart4);

   double pdf_flux = 0.0;
   for (const auto& fv : flav_) {
      pdf_flux += fv.at(2) * pdf_->xfxQ(fv.at(0), x1, muF_) * pdf_->xfxQ(fv.at(1), x2, muF_);
   }
   pdf_flux /= (x1 * x2);

   return dSigmaHad*jacobian*to_fb * pdf_flux;   // in femto barn
}


 std::array<double, 3> XSection_Virt::integrate() {
   constexpr int ndim = 3;
   constexpr int ncomp = 1;
   constexpr int nvec = 1;
   constexpr double accuracy_abs = 1e-12;
   constexpr int eval_min = 1000;
   constexpr int eval_max = 1000000;
   int nregions, neval, fail;
   double integral[ncomp], error[ncomp], prob[ncomp];
   constexpr int last = 4;
   constexpr int key = 0;
   // Divonne specific
   constexpr int seed = 1;
   constexpr double border = 1e-6; // don't go lower than 1e-6 for 1e-4 relative accuraccy
   constexpr int key1 = 47;
   constexpr int key2 = 1;
   constexpr int key3 = 1;
   constexpr int maxpass = 5;
   constexpr double maxchisq = 10.;
   constexpr double mindeviation = .25;
   constexpr int ngiven = 0;
   constexpr int ldxgiven = ndim;
   constexpr int nextra = 0;
   const double prec_virt = std::pow(10., -integration_precision_);

   static bool looptools_initialized = false;
   if(!looptools_initialized) {
      std::cout << std::endl;
      ltini();
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

   std::array<double, 3> result {integral[0], error[0], prob[0]};

   return result;
}
