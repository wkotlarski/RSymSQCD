#ifndef XSECTION_H_
#define XSECTION_H_

#include "LHAPDF/LHAPDF.h"
#include "include/cuba.h"

#include "constants.hpp"
#include "mathematica_wrapper.hpp"

class XSection {

  public:
    virtual std::array<double, 3> integrate() = 0;
    static void init (void);

  protected:
    static constexpr double S_sqrt { 13e+3 };
    static constexpr double S { pow( S_sqrt, 2) };
    static constexpr double gluino_width { 0 };
    static constexpr double gluino_mass { 1e+3 };
    //
    static std::array< std::array<double, 2>, 6 > squark_mass;
    static constexpr std::array< std::array<double, 2>, 6 > squark_width {{
      {{0., 0.}},
      {{0., 0.}},
      {{0., 0.}},
      {{0., 0.}},
      {{0., 0.}},
      {{0., 0.}}
    }};
    static constexpr double sgluon_mass = 1e+4;
    static constexpr std::array<double, 2> sgluons_width {{ 0., 0.}};
    // 5-flavor scheme (only massive top)
    static constexpr double top_quark_mass = 173.21;
    static const LHAPDF::PDF* pdf_lo;
    static const LHAPDF::PDF* pdf_nlo;
    /*
     *  TODO: add random generation of phase space point which can be used to check
     *  IR finiteness
     */
};

#endif /* XSECTION_H_ */
