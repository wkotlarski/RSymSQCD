#ifndef XSECTION_H_
#define XSECTION_H_

#include <array>
#include <vector>
#include <functional>

namespace LHAPDF {
   class PDF;
}

struct XSectionParameters {
   double sqrtS;
   double muR;
   double muF;
   LHAPDF::PDF* pdf;
};

class XSection {
public:
   XSection(XSectionParameters const& parameters, double m1, double m2, std::vector<std::array<int, 3>> const& flav, int integration_precision, int integration_verbosity)
      : m1_(m1), m2_(m2),
        sqrtS_(parameters.sqrtS), muR_(parameters.muR), muF_(parameters.muF),
        flav_(flav), pdf_(parameters.pdf),
        integration_precision_(integration_precision), integration_verbosity_(integration_verbosity)
   {};

protected:

   const double sqrtS_;
   const double m1_;
   const double m2_;
   const double muR_;
   const double muF_;
   const std::vector<std::array<int, 3>> flav_ {};
   const LHAPDF::PDF* const pdf_;
   const int integration_verbosity_;
   const int integration_precision_;
};

#endif // XSECTION_H_
