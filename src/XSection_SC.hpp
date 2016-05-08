#ifndef SRC_XSECTION_SC_H_
#define SRC_XSECTION_SC_H_

#include "XSection_Real.hpp"
#include "dilog.hpp"

class XSection_SC: public XSection_Real {

  public:
    std::array<double, 3> integrate();

  private:
    static int integrand_sc(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata);
    static int integrand_c(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata);
};

#endif /* SRC_XSECTION_SC_H_ */
