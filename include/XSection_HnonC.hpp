#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include "XSection.hpp"

#include <array>

class XSection_HnonC : public XSection {
public:
   XSection_HnonC(
      XSectionParameters const& parameters,
      double m1, double m2,
      std::function<double(double, std::array<std::array<double, 4>, 5>)> f_,
      double dS, double dC,
      std::vector<std::array<int, 3>> const& flav,
      int integration_precision, int integration_verbosity
   ) : XSection{parameters, m1, m2, flav, integration_precision, integration_verbosity}, f{f_}, dS_{dS}, dC_{dC}
   {};

   std::array<double, 3> integrate();
   double integrand(std::array<double, 7> const&);

private:
   const double dS_;
   const double dC_;
   std::function<double(double, std::array<std::array<double, 4>, 5>)> f;
};

#endif // XSECTION_HNONC_H_
