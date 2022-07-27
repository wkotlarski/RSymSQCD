#ifndef XSECTION_VIRT_H_
#define XSECTION_VIRT_H_

#include "XSection.hpp"

class XSection_Virt : public virtual XSection {
public:
   XSection_Virt(
      XSectionParameters const& parameters,
      double m1, double m2,
      std::function<double(double, double, double, double, double, int, double)> f_,
      std::vector<std::array<int, 3>> flav
   ) : XSection(parameters, m1, m2, flav), f(f_) {};

   std::array<double, 3> integrate();
   int integrand(const int *ndim, const double xx[],
      const int *ncomp, double ff[], void *userdata);

private:
   std::function<double(double, double, double, double, double, int, double)> f;
};

#endif // XSECTION_VIRT_H_
