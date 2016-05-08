#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include "XSection_Real.hpp"
#include "constants.hpp"

// neede to do Euler rotation
#include "rk/rk.hh"
#include "rk/geom3.hh"
// hard-non collinear ME
#include "Process_uu_ulurg.h"


class XSection_HnonC : public virtual XSection_Real {

  public:
    std::array<double, 3> integrate();

  private:
    static int integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata);
};
#endif /* XSECTION_HNONC_H_ */
