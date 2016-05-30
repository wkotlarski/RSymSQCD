#ifndef XSECTION_TREE_MRSSM_H_
#define XSECTION_TREE_MRSSM_H_

#include <iostream>
#include "XSection.h"

#include "LHAPDF/LHAPDF.h"


class XSection_Tree_MRSSM : public virtual XSection {
  private:
  static  int integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata);

  public:
    XSection_Tree_MRSSM();
    virtual ~XSection_Tree_MRSSM();
    std::array<double, 3> integrate();
    void show_settings();
};



#endif /* XSECTION_TREE_MRSSM_H_ */
