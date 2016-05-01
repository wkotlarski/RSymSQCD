#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include "XSection.h"

#include "constants.hpp"

// neede to do Euler rotations
#include "rk/rk.hh"
#include "rk/geom3.hh"

class XSection_HnonC : public virtual XSection {

  public:
    static double dS;
    static double dC;
    std::array<double, 3> integrate();

  private:
    static int integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata);
};

#endif /* XSECTION_HNONC_H_ */
