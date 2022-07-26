#ifndef XSECTION_TREE_MRSSM_H_
#define XSECTION_TREE_MRSSM_H_

#include "XSection.hpp"

#include <functional>

class XSection_Tree : public XSection {
public:
   XSection_Tree(double m1_, double m2_, std::function<double(double, double, double)> f_, std::vector<std::array<int, 3>> flav, double muR_, double muF_)
      : m1(m1_), m2(m2_), f(f_), flav_(flav), muR(muR_), muF(muF_) {};
   std::array<double, 3> integrate();
   void show_settings();
   int integrand(const int *ndim, const double xx[],
      const int *ncomp, double ff[], void *userdata);

private:
   const double m1;
   const double m2;
   const double muR;
   const double muF;
   std::function<double(double, double, double)> f;
   std::vector<std::array<int, 3>> flav_ {};
};

#endif // XSECTION_TREE_MRSSM_H_
