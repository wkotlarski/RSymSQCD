#ifndef XSECTION_H_
#define XSECTION_H_

#include <array>
#include "LHAPDF/LHAPDF.h"
#include "include/cuba.h"

#include "mathematica_wrapper.hpp"

class XSection {
  protected:
    static constexpr double S_sqrt = 13000;
    static constexpr double S = pow( S_sqrt, 2);
    static constexpr double gluino_width = 0;
    static constexpr double gluino_mass = 1000;
    //
    static constexpr std::array< std::array<double, 2>, 6 > squark_mass {{
          {{1500, 1500}},
          {{1500, 1500}},
          {{1500, 1500}},
          {{1500, 1500}},
          {{1500, 1500}},
          {{1500, 1500}}
        }};
    static constexpr std::array< std::array<double, 2>, 6 > squark_width {{
      {{0., 0.}},
      {{0., 0.}},
      {{0., 0.}},
      {{0., 0.}},
      {{0., 0.}},
      {{0., 0.}}
    }};
    static constexpr double sgluon_mass = 1e+3;
    static constexpr std::array<double, 2> sgluons_width {{ 0., 0.}};
    // 5-flavor scheme (only massive top)
    static constexpr double top_quark_mass = 173.21;
    static const LHAPDF::PDF* pdf;
  public:
    virtual double integrate() = 0;
};

#endif /* XSECTION_H_ */
