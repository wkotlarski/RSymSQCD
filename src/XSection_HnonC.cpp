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
    ff[0] = 2.;
    return 0;
}

double XSection_HnonC::integrate() {
  const int ndim = 7;
  const int ncomp = 1;
  constexpr double accuracy_rel = 1e-3;
  constexpr double accuracy_abs = 1e-12;
  constexpr int eval_min = 0;
  constexpr int nstart = 0;
  constexpr int nincrease = 100;
  constexpr int nbatch = 1000;
  constexpr int gridno = 0;
  const char* state_file = "";

  long long int eval_max = 1e+5;//strtoll( "1e+3", NULL, 10 );
  long long int neval;
  int nregions, fail;
  cubareal integral[ncomp], error[ncomp], prob[ncomp];
  llVegas(ndim, ncomp, integrand, NULL, 1,
         accuracy_rel, accuracy_abs, 8 | 1, 0,
         0, 1000, 100, nincrease, nbatch,
         gridno, state_file, NULL,
         &neval, &fail, integral, error, prob);
  return (double)integral[0];
}
