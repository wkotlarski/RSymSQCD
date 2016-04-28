/*
 * XSection_HnonC.cpp
 *
 *  Created on: 22 kwi 2016
 */

#include "XSection_HnonC.h"

XSection_HnonC::XSection_HnonC() {

  // TODO Auto-generated constructor stub

}

XSection_HnonC::~XSection_HnonC() {
  // TODO Auto-generated destructor stub
}

//double XSection_HnonC::integrate() {
//  return 1.;
//}

void XSection_HnonC::show_settings() {
  std::cout << S_sqrt << std::endl;
        double x1 = 0.1;
        double m = 1000;
        std::cout << pdf->xfxQ(21, x1, m)/x1 << std::endl;
}


int XSection_HnonC::integrand(const int *ndim, const cubareal xx[],
              const int *ncomp, cubareal ff[], void *userdata) {
    
    return 1;
}
