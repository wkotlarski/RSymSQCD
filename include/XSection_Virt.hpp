#ifndef XSECTION_VIRT_H_
#define XSECTION_VIRT_H_

#include "XSection.hpp"

class XSection_Virt : public virtual XSection {
public:
   XSection_Virt(
      double m1, double m2,
      std::function<double(double, double, double, double, double, int, double)> f_,
      std::vector<std::array<int, 3>> flav,
      double muR, double muF,
      const LHAPDF::PDF* const pdf
   ) : XSection(m1, m2, muR, muF, flav, pdf), f(f_) {};
   std::array<double, 3> integrate();
   int integrand(const int *ndim, const double xx[],
      const int *ncomp, double ff[], void *userdata);

private:
   std::function<double(double, double, double, double, double, int, double)> f;
};

#endif // XSECTION_VIRT_H_
