#ifndef SRC_XSECTION_SC_H_
#define SRC_XSECTION_SC_H_

#include "cuba.h"

#include "XSection.hpp"

class XSection_SC: public XSection {

  public:
    std::array<double, 3> integrate();
    static double f(double*, double*, double*, double*);

  private:
    static int integrand_sc(const int*, const cubareal[],
              const int*, cubareal[], void*);
    static int integrand_c1(const int*, const cubareal[],
              const int*, cubareal[], void*);
    static int integrand_c2(const int*, const cubareal[],
              const int*, cubareal[], void*);
};

#endif // SRC_XSECTION_SC_H_
