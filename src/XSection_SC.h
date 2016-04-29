/*
 * XSectionSC.h
 *
 *  Created on: 28 kwi 2016
 *      Author: Navir
 */

#ifndef SRC_XSECTION_SC_H_
#define SRC_XSECTION_SC_H_

#include "XSection.h"

class XSection_SC: public XSection {
  public:
    XSection_SC();
    virtual ~XSection_SC();
    double integrate();
};

#endif /* SRC_XSECTION_SC_H_ */
