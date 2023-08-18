#ifndef SGLUONS_H_
#define SGLUONS_H_

#include <array>

struct SgluonParameters {
   double mO;
   double mt;
};

class Sgluons {
public:
   Sgluons(SgluonParameters const&);
      double matrixSgluonsTree_qqbar_OO(double, double, double, int) {return 0.;}
      double matrixSgluonsTree_gg_OO(double, double, double, int) {return 0.;}
      double sigmaSgluonsTree_qqbar_OO(double, double, double = 0.);
      double sigmaSgluonsTree_gg_OO(double, double, double = 0.);

      double sgluons_qqbar_OOg_soft(double, double, double, double, double) const;
      double sgluons_qqbar_OOg_hard(double, std::array<std::array<double, 4>, 5> const&) const;
      double sgluons_gg_OOg_soft(double, double, double, double, double) const;
      double sgluons_gg_OOg_hard(double, std::array<std::array<double, 4>, 5> const&) const;

private:
   const double mO;
   const double mt;
};

#endif
