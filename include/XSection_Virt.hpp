#ifndef XSECTION_VIRT_H_
#define XSECTION_VIRT_H_

#include "cuba.h"

#include "XSection.hpp"

class XSection_Virt : public virtual XSection {
   private:
   static int integrand(const int *ndim, const cubareal xx[],
      const int *ncomp, cubareal ff[], void *userdata);

   public:
      std::array<double, 3> integrate();
};

#endif // XSECTION_VIRT_H_
