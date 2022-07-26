#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include "XSection.hpp"

#include <array>

class XSection_HnonC : public XSection {
public:
   XSection_HnonC(double m1_, double m2_, std::function<double(double, std::array<std::array<double, 4>, 5>)> f_, std::vector<std::array<int, 3>> const& flav) : m1(m1_), m2(m2_), f(f_), flav_(flav) {};
   std::array<double, 3> integrate();
   int integrand(const int *ndim, const double xx[],
                 const int *ncomp, double ff[], void *userdata);

private:
   const double m1;
   const double m2;
   std::function<double(double, std::array<std::array<double, 4>, 5>)> f;
   std::vector<std::array<int, 3>> flav_ {};
};

#endif // XSECTION_HNONC_H_
