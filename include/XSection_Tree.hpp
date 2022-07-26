#ifndef XSECTION_TREE_MRSSM_H_
#define XSECTION_TREE_MRSSM_H_

#include "XSection.hpp"

#include <functional>

class XSection_Tree : public XSection {
public:
   XSection_Tree(
      double m1, double m2,
      std::function<double(double, double, double)> f_,
      std::vector<std::array<int, 3>> flav,
      double muR, double muF,
      const LHAPDF::PDF* const pdf
   ) : XSection(m1, m2, muR, muF, flav, pdf), f(f_) {};

   std::array<double, 3> integrate();
   void show_settings();
   int integrand(const int *ndim, const double xx[],
      const int *ncomp, double ff[], void *userdata);

private:
   std::function<double(double, double, double)> f;
};

#endif // XSECTION_TREE_MRSSM_H_
