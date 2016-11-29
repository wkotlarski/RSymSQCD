#include "XSection_Virt.hpp"

int XSection_Virt::integrand(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

    double mu_f = pt.get<double>("collider setup.mu_f");
    double x1min = 4. * pow( m1, 2 )/S;
    double xmax = 1.;
    double x1 = x1min + (xmax - x1min) * xx[1];
    double x2min = 4. * pow( m1, 2 )/(S*x1);
    double x2 = x2min + (xmax - x2min) * xx[2];
    double s = S * x1 * x2;     //partonic 
    double Tmin = pow( m1, 2 ) - s/2. - sqrt( pow(s, 2)/4 -
                  pow( m1, 2 )*s);
    double Tmax = pow( m1, 2 ) - s/2. + sqrt( pow(s, 2)/4. -
                  pow( m1, 2 )*s);
    double T = xx[0]*(Tmax-Tmin) + Tmin;
    double jacobian = (Tmax-Tmin)*(1.-x1min)*(1.-x2min);

    double FiniteGs = 1;
    double Dminus4 = 0;
    int Divergence = 0;     // O(eps) 
    double squaredMReal = (processID->*processID->matrixelementVirt)(
      s, T, FiniteGs, Dminus4, Divergence);
    
    double dSigmaPart1 = 2.*squaredMReal*(processID->h)*M_PI/(pow(4.*M_PI,2))/
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
                         (processID->h)*M_PI/(pow(4.*M_PI,2))/
                         (processID->k)/(pow(s,2));

    // contraction with O(eps^2) prefactor of loop integral
    // and with product of O(eps) prefactors of phase space and loop integral
    Divergence = -2;
    Dminus4 = 0;
   squaredMReal = (processID->*processID->matrixelementVirt)(
      s, T, FiniteGs, Dminus4, Divergence);
    
    double dSigmaPart4 = 2.*squaredMReal*(processID->h)*M_PI/(pow(4.*M_PI,2))/
                         (processID->k)/(pow(s,2))
                         *(pow(M_PI,2.)/6.);

    double dSigmaHad = (dSigmaPart1 + dSigmaPart3 + dSigmaPart4)
                     * pdf->xfxQ( processID->f1, x1, mu_f )/x1
                     * pdf->xfxQ( processID->f2, x2, mu_f )/x2;
//////////////////////////////////////////////////////////////////////
// here ends the part without prefactors
//////////////////////////////////////////////////////////////////////

    ff[0] = dSigmaHad*jacobian*to_fb;   // in femto barn
    return 1;
}


 std::array<double, 3> XSection_Virt::integrate() {
    constexpr int ndim = 3;
    constexpr int ncomp = 1;
    constexpr int nvec = 1;
    //constexpr double accuracy_rel = 1e-3;
    constexpr double accuracy_abs = 1e-12;
    constexpr int eval_min = 1000;
    constexpr int eval_max = 1000000;
    constexpr int verbose = 0;        // adjust shown output 0 ... 3
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];
    // Vegas specific
    //int seed = 1;
    //int nstart = 1000;
    //int nincrease = 500;
    //int nbatch = 1000;
    //int grindo = 0;
    // Cuhre specific
    constexpr int last = 4;
    constexpr int key = 0;

    //  Vegas(ndim, ncomp, integrand, NULL, nvec,
    //  accuracy_rel, accuracy_abs, verbose, seed,
    //  eval_min, eval_max, nstart, nincrease, nbatch, 
    //  grindo, NULL, NULL,
    //  &neval, &fail, integral, error, prob);

      Cuhre(ndim, ncomp, integrand, NULL, nvec,
      prec_virt, accuracy_abs, verbose | last,
      eval_min, eval_max, key, NULL, NULL,
      &nregions, &neval, &fail, integral, error, prob);

    std::array <double, 3> result{integral[0], error[0], prob[0]}; 

    return result;
}
