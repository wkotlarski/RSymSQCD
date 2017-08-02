#ifndef XSECTION_VIRT_H_
#define XSECTION_VIRT_H_

#include "LHAPDF/LHAPDF.h"

#include "XSection.hpp"
#include "CSDipole.hpp"

class XSection_Virt : public virtual XSection {

public:
   std::array<double, 3> integrate();
   static std::vector<CSDipole> cs_dipoles;

private:
   static int integrand(const int *ndim, const cubareal xx[],
      const int *ncomp, cubareal ff[], void *userdata);

};

#endif /* XSECTION_VIRT_H_ */
