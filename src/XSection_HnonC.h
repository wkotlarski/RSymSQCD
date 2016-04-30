#ifndef XSECTION_HNONC_H_
#define XSECTION_HNONC_H_

#include <iostream>
#include "XSection.h"

#include "LHAPDF/LHAPDF.h"

//  
#include "rk/rk.hh"
#include "rk/geom3.hh"

class XSection_HnonC : public virtual XSection {
  private:
   static int integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata);

  public:
    XSection_HnonC();
    virtual ~XSection_HnonC();
    double integrate();
    void show_settings();
};

#endif /* XSECTION_HNONC_H_ */
