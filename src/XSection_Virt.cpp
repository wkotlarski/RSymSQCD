#include "XSection_Virt.hpp"

int XSection_Virt::integrand(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

    const double x1min = 4. * pow( m1, 2 )/S;
    static constexpr double xmax = 1.;
    const double x1 = x1min + (xmax - x1min) * xx[1];
    const double x2min = 4. * pow( m1, 2 )/(S*x1);
    const double x2 = x2min + (xmax - x2min) * xx[2];
    const double s = S * x1 * x2;     //partonic
    const double Tmin = pow( m1, 2 ) - s/2. - sqrt( pow(s, 2)/4 -
                  pow( m1, 2 )*s);
    const double Tmax = pow( m1, 2 ) - s/2. + sqrt( pow(s, 2)/4. -
                  pow( m1, 2 )*s);
    const double T = xx[0]*(Tmax-Tmin) + Tmin;
    const double jacobian = (Tmax-Tmin)*(1.-x1min)*(1.-x2min);

    int FiniteGs = 1;
    double Dminus4 = 0;
    int Divergence = 0;     // O(eps)

    //using Func = double (Process::* double)(double, double, int, double, int);
    //Func* f_ptr = processID->*processID->matrixelementVirt;

    double squaredMReal = (processID->*processID->matrixelementVirt)(
      s, T, FiniteGs, Dminus4, Divergence);

    double dSigmaPart1 = 2.*squaredMReal*(processID->h)*pi/(pow(4.*pi,2))/
                         (processID->k)/(pow(s,2));

    // contraction with O(eps) from Dminus4
    Divergence = -1;           // O(eps)
    FiniteGs = 0;
    squaredMReal = (processID->*processID->matrixelementVirt)(
      s, T, FiniteGs, Dminus4, Divergence);

    Dminus4 = -2.;
    double squaredMRealMinus2 = (processID->*processID->matrixelementVirt)(
                         s, T, FiniteGs, Dminus4, Divergence);

    double dSigmaPart3 = 2.*(squaredMRealMinus2 - squaredMReal)*
                         (processID->h)*pi/(pow(4.*pi,2))/
                         (processID->k)/(pow(s,2));

    // contraction with O(eps^2) prefactor of loop integral
    // and with product of O(eps) prefactors of phase space and loop integral
    Divergence = -2;
    Dminus4 = 0;
   squaredMReal = (processID->*processID->matrixelementVirt)(
      s, T, FiniteGs, Dminus4, Divergence);

    double dSigmaPart4 = 2.*squaredMReal*(processID->h)*pi/(pow(4.*pi,2))/
                         (processID->k)/(pow(s,2))
                         *(pow(pi,2.)/6.);

    double dSigmaHad = (dSigmaPart1 + dSigmaPart3 + dSigmaPart4);

   double pdf_flux = 0.0;
   for (const auto& flav : processID->flav) {
      pdf_flux += flav.at(2) * pdf->xfxQ( flav.at(0), x1, mu_f ) * pdf->xfxQ( flav.at(1), x2, mu_f );
   }
   pdf_flux /= (x1 * x2);

    ff[0] = dSigmaHad*jacobian*to_fb * pdf_flux;   // in femto barn
    return 0;
}


 std::array<double, 3> XSection_Virt::integrate() {
   constexpr int ndim = 3;
   constexpr int ncomp = 1;
   constexpr int nvec = 1;
   constexpr double accuracy_abs = 1e-12;
   constexpr int eval_min = 1000;
   constexpr int eval_max = 1000000;
   const int verbose = vm["verbosity-virt"].as<int>();        // adjust output 0 ... 3
   int nregions, neval, fail;
   cubareal integral[ncomp], error[ncomp], prob[ncomp];
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
   const double prec_virt = pow( 10., -vm["precision-virt"].as<int>());

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
   Divonne(ndim, ncomp, integrand, NULL, nvec,
        prec_virt, accuracy_abs, verbose, seed,
        eval_min, eval_max, key1, key2, key3, maxpass,
        border, maxchisq, mindeviation,
        ngiven, ldxgiven, NULL, nextra, NULL,
        NULL, NULL,
        &nregions, &neval, &fail, integral, error, prob);

   std::array <double, 3> result {integral[0], error[0], prob[0]};

   return result;
}
