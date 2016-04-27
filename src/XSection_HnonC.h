/*
 * XSection_HnonC.h
 *
 *  Created on: 22 kwi 2016
 */

#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include "XSection.h"

class XSection_HnonC : public virtual XSection {
  private:
    //process

  public:
    XSection_HnonC();
    virtual ~XSection_HnonC();
    void show_settings();
};



#endif /* XSECTION_HNONC_H_ */
