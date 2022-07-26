#ifndef SRC_XSECTION_SC_H_
#define SRC_XSECTION_SC_H_

#include "XSection.hpp"

#include "cuba.h"
#include "splitting_kernels.hpp"

class XSection_SC: public XSection {
public:
   XSection_SC(
      double m1, double m2,
      std::function<double(double, double, double, double, double)> f_soft,
      double dS, double dC,
      double muR, double muF,
      std::vector<std::array<int, 3>> flav,
      const LHAPDF::PDF* const pdf,
      std::array<std::pair<SplittingKernel, std::function<double(double, double)>>, 2> sp
   ) : XSection(m1, m2, muR, muF, flav, pdf), f_soft_(f_soft), dS_(dS), dC_(dC), sp_(sp) {};

   std::array<double, 3> integrate();
   int integrand_sc(const int*, const cubareal[],
              const int*, cubareal[], void*);
   int integrand_c1(const int*, const cubareal[],
              const int*, cubareal[], void*);
   int integrand_c2(const int*, const cubareal[],
              const int*, cubareal[], void*);

private:
   const double dS_;
   const double dC_;
   std::function<double(double, double, double, double, double)> f_soft_;
   std::array<std::pair<SplittingKernel, std::function<double(double, double)>>, 2> sp_;
};

#endif // SRC_XSECTION_SC_H_
