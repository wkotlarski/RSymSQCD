#ifndef SRC_XSECTION_SC_H_
#define SRC_XSECTION_SC_H_

#include "XSection.hpp"
#include "dilog.hpp"
#include "CSDipole.hpp"
// to get min(x,y) function
#include <algorithm>

//extern "C" {
//   double qq_s3Ls3R(double*, double*, double*, double*);
//}

class XSection_SC: public XSection {

public:
   std::array<double, 3> integrate();
   static double f(double*, double*, double*, double*);
   static std::vector<CSDipole> cs_dipoles;

  private:
   static int integrand_sc(const int*, const cubareal[],
              const int*, cubareal[], void*);
   static int integrand_c1(const int*, const cubareal[],
              const int*, cubareal[], void*);
   static int integrand_c2(const int*, const cubareal[],
              const int*, cubareal[], void*);
};

#endif /* SRC_XSECTION_SC_H_ */
