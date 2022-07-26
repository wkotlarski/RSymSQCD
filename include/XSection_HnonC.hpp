#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include "XSection.hpp"

#include <array>

class XSection_HnonC : public XSection {
public:
   XSection_HnonC(
      double m1, double m2,
      std::function<double(double, std::array<std::array<double, 4>, 5>)> f_,
      double dS, double dC,
      double muR, double muF,
      std::vector<std::array<int, 3>> const& flav,
      const LHAPDF::PDF* const pdf
   ) : XSection(m1, m2, muR, muF, flav, pdf), f(f_), dS_(dS), dC_(dC) {};

   std::array<double, 3> integrate();
   int integrand(const int *ndim, const double xx[],
                 const int *ncomp, double ff[], void *userdata);

private:
   const double dS_;
   const double dC_;
   std::function<double(double, std::array<std::array<double, 4>, 5>)> f;
};

#endif // XSECTION_HNONC_H_
