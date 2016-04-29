/*
 * XSection_HnonC.cpp
 *
 *  Created on: 22 kwi 2016
 */
#include <string>
#include "XSection_HnonC.h"

XSection_HnonC::XSection_HnonC() {

  // TODO Auto-generated constructor stub

}

XSection_HnonC::~XSection_HnonC() {
  // TODO Auto-generated destructor stub
}

void XSection_HnonC::show_settings() {
  std::cout << S_sqrt << std::endl;
        double x1 = 0.1;
        double m = 1000;
        std::cout << pdf->xfxQ(21, x1, m)/x1 << std::endl;
}


int XSection_HnonC::integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata) {
    ff[0] = 1;
    return 0;
}

double XSection_HnonC::integrate() {
  constexpr int ndim = 7;
  constexpr int ncomp = 1;
  constexpr int accuracy_rel = 1e-3;
  constexpr int accuracy_abs = 1e-12;
  constexpr int eval_min = 0;
  const int eval_max = 1000;
  constexpr int nstart = 200000;
  constexpr int nincrease = 100;
  constexpr int nbatch = 1000;
  constexpr int gridno = 0;
  //string state_file = ""; //"gg_s8ps8pg_vegas.tmp"

  long long int max_eval = strtoll( "1e+3", NULL, 10 );
  long long int neval;
  int comp, nregions, fail;
  cubareal integral[ncomp], error[ncomp], prob[ncomp];
  llVegas(ndim, ncomp, integrand, NULL, 1,
         accuracy_rel, accuracy_abs, 8 | 1, 0,
         eval_min, eval_max, nstart, nincrease, nbatch,
         gridno, "", NULL, &neval, &fail, integral, error, prob);
  return 1.;//((double)integral[0]);
}
