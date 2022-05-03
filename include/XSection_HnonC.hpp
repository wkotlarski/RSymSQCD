#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include "cuba.h"

#include "XSection.hpp"

class XSection_HnonC : public virtual XSection {

  public:
    std::array<double, 3> integrate();

  private:
    static int integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata);
};

#endif // XSECTION_HNONC_H_
