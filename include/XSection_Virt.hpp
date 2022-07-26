#ifndef XSECTION_VIRT_H_
#define XSECTION_VIRT_H_

#include "XSection.hpp"

class XSection_Virt : public virtual XSection {
public:
   XSection_Virt(
      double m1_, double m2_,
      std::function<double(double, double, double, double, double, int, double)> f_,
      std::vector<std::array<int, 3>> flav,
      double muR_, double muF_
   ) : m1(m1_), m2(m2_), f(f_), flav_(flav), muR(muR_), muF(muF_) {};
   std::array<double, 3> integrate();
   int integrand(const int *ndim, const double xx[],
      const int *ncomp, double ff[], void *userdata);

private:
   const double m1;
   const double m2;
   const double muR;
   const double muF;
   std::function<double(double, double, double, double, double, int, double)> f;
   std::vector<std::array<int, 3>> flav_ {};
};

#endif // XSECTION_VIRT_H_
