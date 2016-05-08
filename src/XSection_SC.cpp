#include "XSection_SC.hpp"

std::array<double, 3> XSection_SC::integrate() {
  
  //  integral dimension, number of integrands
  constexpr int ndim { 3 }, ncomp { 1 };
  //  accuraccy
  constexpr double accuracy_rel { 1e-5 }, 
          accuracy_abs { 1e-12 };

  constexpr int neval_min = 10000;
  long long int neval;
  constexpr long long int neval_max { 1000000000 }; 
    // @TODO: read from external source strtoll( "1e+3", NULL, 10 );

  // technical (Vegas specific) stuff
  constexpr int nstart      = 200000;
  constexpr int nincrease   = 100;
  constexpr int nbatch      = 1000;
  constexpr int gridno      = 0;
  const char* state_file    = "";
  int nregions, fail;

  cubareal integral_sc[ncomp], error_sc[ncomp], prob_sc[ncomp];
  llVegas( ndim, ncomp, integrand_sc, NULL, 1,
           accuracy_rel, accuracy_abs, 8 | 1, 0,
           neval_min, neval_max, nstart, nincrease, nbatch,
           gridno, state_file, NULL,
           &neval, &fail, integral_sc, error_sc, prob_sc );

  cubareal integral_c[ncomp], error_c[ncomp], prob_c[ncomp];
  llVegas( ndim, ncomp, integrand_sc, NULL, 1,
           accuracy_rel, accuracy_abs, 8 | 1, 0,
           neval_min, neval_max, nstart, nincrease, nbatch,
           gridno, state_file, NULL,
           &neval, &fail, integral_c, error_c, prob_c );
  
  std::array <double, 3> result_finite { 
      integral_sc[0] + integral_c[0], 
      sqrt( pow( error_sc[0], 2) + pow ( error_c[0], 2) ),
      prob_sc[0] + prob_c[0] 
  };
  
  return result_finite;
}

int XSection_SC::integrand_sc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
    ff[0] = 1;
  return 0;
}

// @todo if muR != muF one needs one more term here

int XSection_SC::integrand_c(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
    ff[0] = 1;
  return 0;
}