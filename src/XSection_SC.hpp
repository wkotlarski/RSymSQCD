#ifndef SRC_XSECTION_SC_H_
#define SRC_XSECTION_SC_H_

#include "XSection_Real.hpp"
// @todo: we probably cannot redistribute this package
#include "dilog.hpp"
// to get min(x,y) function
#include <algorithm>

class XSection_SC: public XSection_Real {

  public:
    std::array<double, 3> integrate();

  private:
    static int integrand_sc(const int*, const cubareal[],
              const int*, cubareal[], void*);
    static int integrand_c1(const int*, const cubareal[],
              const int*, cubareal[], void*);
    static int integrand_c2(const int*, const cubareal[],
              const int*, cubareal[], void*);
};

#endif /* SRC_XSECTION_SC_H_ */
