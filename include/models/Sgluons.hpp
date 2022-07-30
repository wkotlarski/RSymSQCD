#ifndef SGLUONS_H_
#define SGLUONS_H_

#include <array>

struct SgluonParameters {
   double mt;
   double mO;
};

class Sgluons {
public:
   Sgluons(SgluonParameters const&);
      double matrixSgluonsTree_qqbar_OO(double, double, double) {return 0.;}
      double matrixSgluonsTree_gg_OO(double, double, double);
      double sigmaSgluonsTree_qqbar_OO(double, double);

      double sgluons_qqbar_OOg_soft(double, double, double, double, double) const;
      double sgluons_qqbar_OOg_hard(double, std::array<std::array<double, 4>, 5> const&) const;
      double sgluons_gg_OOg_soft(double, double, double, double, double) const;
      double sgluons_gg_OOg_hard(double, std::array<std::array<double, 4>, 5> const&) const;

private:
   double mO;
};

#endif
