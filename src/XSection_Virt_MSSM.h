#ifndef XSECTION_VIRT_MSSM_H_
#define XSECTION_VIRT_MSSM_H_

#include <iostream>
#include "XSection.h"

#include "LHAPDF/LHAPDF.h"


class XSection_Virt_MSSM : public virtual XSection {
  private:
  static  int integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata);

  public:
    XSection_Virt_MSSM();
    virtual ~XSection_Virt_MSSM();
    std::array<double, 3> integrate();
    void show_settings();
};



#endif /* XSECTION_VIRT_MSSM_H_ */
