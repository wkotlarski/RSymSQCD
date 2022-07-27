#ifndef SRC_XSECTION_SC_H_
#define SRC_XSECTION_SC_H_

#include "XSection.hpp"

#include "cuba.h"
#include "splitting_kernels.hpp"

class XSection_SC: public XSection {
public:
   XSection_SC(
      XSectionParameters const& parameters,
      double m1, double m2,
      std::function<double(double, double, double, double, double)> f_soft,
      double dS, double dC,
      std::vector<std::array<int, 3>> flav,
      std::array<std::pair<SplittingKernel, std::function<double(double, double)>>, 2> sp,
      int integration_precision, int integration_verbosity
   ) : XSection(parameters, m1, m2, flav, integration_precision, integration_verbosity), f_soft_(f_soft), dS_(dS), dC_(dC), sp_(sp) {};

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
