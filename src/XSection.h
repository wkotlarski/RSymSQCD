/*
 * XSection.h
 *
 *  Created on: 27 kwi 2016
 */

#ifndef XSECTION_H_
#define XSECTION_H_

#include "LHAPDF/LHAPDF.h"

class XSection {
  protected:
    double S_sqrt = 13000;
    double gluino_width = 0;
    double gluino_mass = 1000;
    //
    std::array<double, 12> squark_mass {
      {1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500, 1500}};
    std::array< std::array<double, 2>, 6 > squark_width {{
      {{0., 0.}}, {{0., 0.}}, {{0., 0.}}, {{0., 0.}}, {{0., 0.}}, {{0., 0.}}
    }};
    double sgluon_mass;
    std::array<double, 2> sgluons_width {{ 0., 0.}};
    double top_quark_mass = 173.21;
    const LHAPDF::PDF* pdf;
  public:
    XSection();
    virtual ~XSection();
};

#endif /* XSECTION_H_ */
