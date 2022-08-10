#ifndef XSECTION_VIRT_H_
#define XSECTION_VIRT_H_

#include "XSection.hpp"

class XSection_Virt : public virtual XSection {
public:
   XSection_Virt(
      XSectionParameters const& parameters,
      double m1, double m2,
      std::function<double(double, double, double, double, double, int, double)> f_,
      std::vector<std::array<int, 3>> flav,
      int integration_precision, int integration_verbosity
   ) : XSection{parameters, m1, m2, flav, integration_precision, integration_verbosity}, f{f_}
   {};

   std::array<double, 3> integrate();
   double integrand(const double xx[]);

private:
   std::function<double(double, double, double, double, double, int, double)> f;
};

#endif // XSECTION_VIRT_H_
