#ifndef SRC_XSECTION_SC_H_
#define SRC_XSECTION_SC_H_

#include "cuba.h"

#include "XSection.hpp"

class XSection_SC: public XSection {
public:
   XSection_SC(double m1_, double m2_, std::function<double(double, double, double, double, double)> f_soft, double dS, double muR, std::vector<std::array<int, 3>> flav)
      : m1(m1_), m2(m2_), f_soft_(f_soft), flav_(flav), dS_(dS), muR_(muR) {};
   std::array<double, 3> integrate();
   int integrand_sc(const int*, const cubareal[],
              const int*, cubareal[], void*);
   int integrand_c1(const int*, const cubareal[],
              const int*, cubareal[], void*);
   int integrand_c2(const int*, const cubareal[],
              const int*, cubareal[], void*);

private:
   const double m1;
   const double m2;
   const double dS_;
   const double muR_;
   std::vector<std::array<int, 3>> flav_ {};
   std::function<double(double, double, double, double, double)> f_soft_;
};

#endif // SRC_XSECTION_SC_H_
