#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include "XSection.hpp"
#include "constants.hpp"
#include "CSDipole.hpp"

// neede to do Euler rotation
#include "rk/rk.hh"
#include "rk/geom3.hh"

class XSection_HnonC : public virtual XSection {

public:
   XSection_HnonC() = default;
   std::array<double, 3> integrate();
   static std::vector<CSDipole> cs_dipoles;

private:
   static int integrand(const int *ndim, const cubareal xx[],
                        const int *ncomp, cubareal ff[], void *userdata);
};

#endif /* XSECTION_HNONC_H_ */
