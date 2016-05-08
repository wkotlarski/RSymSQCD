#ifndef XSECTION_REAL_HPP
#define XSECTION_REAL_HPP

#include "XSection.h"

class XSection_Real : public virtual XSection {
  /*
   *  This class keeps things common to hard and soft 
   *  real emission corrections
   */
  protected:
    static double dS;
    static double dC;

};
#endif /* XSECTION_REAL_HPP */

