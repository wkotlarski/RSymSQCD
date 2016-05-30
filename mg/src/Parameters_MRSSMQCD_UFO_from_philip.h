//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_MRSSMQCD_UFO_from_philip_H
#define Parameters_MRSSMQCD_UFO_from_philip_H

#include <complex> 

#include "read_slha.h"
using namespace std; 

class Parameters_MRSSMQCD_UFO_from_philip
{
  public:

    static Parameters_MRSSMQCD_UFO_from_philip * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double mdl_wOp, mdl_wOs, mdl_wdglu, mdl_WH, mdl_WW, mdl_WZ, mdl_WT,
        mdl_ymtau, mdl_ymt, aS, mdl_Gf, aEWM1, mdl_MD3, mdl_mOp, mdl_mOs,
        mdl_mstr, mdl_msbr, mdl_mscr, mdl_mssr, mdl_msur, mdl_msdr, mdl_mstl,
        mdl_msbl, mdl_mscl, mdl_mssl, mdl_msul, mdl_msdl, mdl_MH, mdl_MZ,
        mdl_Mnu, mdl_MTA, mdl_MT, mdl_yb, mdl_MZ__exp__2, mdl_MZ__exp__4,
        mdl_sqrt__2, mdl_MH__exp__2, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee,
        mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw,
        mdl_vev, mdl_vev__exp__2, mdl_lam, mdl_wein, mdl_yt, mdl_ytau, mdl_muH,
        mdl_ee__exp__2, mdl_sw__exp__2, mdl_cw__exp__2;
    std::complex<double> mdl_I1a33, mdl_I4a33, mdl_complexi, mdl_I2a33,
        mdl_I3a33;
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G; 
    std::complex<double> mdl_G__exp__2; 
    // Model couplings independent of aS

    // Model couplings dependent on aS
    std::complex<double> GC_17, GC_15, GC_14, GC_13, GC_12, GC_11, GC_10; 

    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_MRSSMQCD_UFO_from_philip * instance; 
}; 

#endif  // Parameters_MRSSMQCD_UFO_from_philip_H

