//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "Parameters_MRSSMQCD_UFO_from_philip.h"

// Initialize static instance
Parameters_MRSSMQCD_UFO_from_philip * Parameters_MRSSMQCD_UFO_from_philip::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_MRSSMQCD_UFO_from_philip * Parameters_MRSSMQCD_UFO_from_philip::getInstance()
{
  if (instance == 0)
    instance = new Parameters_MRSSMQCD_UFO_from_philip(); 

  return instance; 
}

void Parameters_MRSSMQCD_UFO_from_philip::setIndependentParameters(SLHAReader&
    slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  mdl_wOp = slha.get_block_entry("decay", 3000022, 1.000000e+01); 
  mdl_wOs = slha.get_block_entry("decay", 3000021, 1.000000e+01); 
  mdl_wdglu = slha.get_block_entry("decay", 1000021, 1.000000e+01); 
  mdl_WH = slha.get_block_entry("decay", 25, 4.070000e-03); 
  mdl_WW = slha.get_block_entry("decay", 24, 2.085000e+00); 
  mdl_WZ = slha.get_block_entry("decay", 23, 2.495200e+00); 
  mdl_WT = slha.get_block_entry("decay", 6, 1.508336e+00); 
  mdl_ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  mdl_ymt = slha.get_block_entry("yukawa", 6, 1.720000e+02); 
  aS = slha.get_block_entry("sminputs", 3, 1.184000e-01); //alphaS
  mdl_Gf = slha.get_block_entry("sminputs", 2, 1.166370e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.279000e+02); 
  mdl_MD3 = slha.get_block_entry("msoft", 1, 1.000000e+03);//gluino 
  mdl_mOp = slha.get_block_entry("mass", 3000022, 1.000000e+03); //sgluon
  mdl_mOs = slha.get_block_entry("mass", 3000021, 1.000000e+03); //sgluon
  mdl_mstr = slha.get_block_entry("mass", 2000006, 1.000000e+03); 
  mdl_msbr = slha.get_block_entry("mass", 2000005, 1.000000e+03); 
  mdl_mscr = slha.get_block_entry("mass", 2000004, 1.000000e+03); 
  mdl_mssr = slha.get_block_entry("mass", 2000003, 1.000000e+03); 
  mdl_msur = slha.get_block_entry("mass", 2000002, 1.000000e+03); 
  mdl_msdr = slha.get_block_entry("mass", 2000001, 1.000000e+03); 
  mdl_mstl = slha.get_block_entry("mass", 1000006, 1.000000e+03); 
  mdl_msbl = slha.get_block_entry("mass", 1000005, 1.000000e+03); 
  mdl_mscl = slha.get_block_entry("mass", 1000004, 1.000000e+03); 
  mdl_mssl = slha.get_block_entry("mass", 1000003, 1.000000e+03); 
  mdl_msul = slha.get_block_entry("mass", 1000002, 1.000000e+03); 
  mdl_msdl = slha.get_block_entry("mass", 1000001, 1.000000e+03); 
  mdl_MH = slha.get_block_entry("mass", 25, 1.250000e+02); 
  mdl_MZ = slha.get_block_entry("mass", 23, 9.118760e+01); 
  mdl_Mnu = slha.get_block_entry("mass", 16, 1.000000e-09); 
  mdl_MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  mdl_MT = slha.get_block_entry("mass", 6, 1.720000e+02); 
  mdl_yb = 0.; 
  mdl_MZ__exp__2 = pow(mdl_MZ, 2.); 
  mdl_MZ__exp__4 = pow(mdl_MZ, 4.); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_MH__exp__2 = pow(mdl_MH, 2.); 
  mdl_I1a33 = mdl_yb; 
  mdl_I4a33 = mdl_yb; 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_aEW = 1./aEWM1; 
  mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = sqrt(mdl_aEW); 
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
  mdl_MW__exp__2 = pow(mdl_MW, 2.); 
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
  mdl_cw = sqrt(1. - mdl_sw2); 
  mdl_sqrt__sw2 = sqrt(mdl_sw2); 
  mdl_sw = mdl_sqrt__sw2; 
  mdl_g1 = mdl_ee/mdl_cw; 
  mdl_gw = mdl_ee/mdl_sw; 
  mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee; 
  mdl_vev__exp__2 = pow(mdl_vev, 2.); 
  mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2); 
  mdl_wein = mdl_Mnu/mdl_vev__exp__2; 
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev; 
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev; 
  mdl_muH = sqrt(mdl_lam * mdl_vev__exp__2); 
  mdl_I2a33 = mdl_yt; 
  mdl_I3a33 = mdl_yt; 
  mdl_ee__exp__2 = pow(mdl_ee, 2.); 
  mdl_sw__exp__2 = pow(mdl_sw, 2.); 
  mdl_cw__exp__2 = pow(mdl_cw, 2.); 
}
void Parameters_MRSSMQCD_UFO_from_philip::setIndependentCouplings()
{

}
void Parameters_MRSSMQCD_UFO_from_philip::setDependentParameters()
{
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = pow(G, 2.); 
}
void Parameters_MRSSMQCD_UFO_from_philip::setDependentCouplings()
{
  GC_17 = mdl_complexi * mdl_G__exp__2; 
  GC_15 = mdl_complexi * G * mdl_sqrt__2; 
  GC_14 = -(mdl_complexi * G * mdl_sqrt__2); 
  GC_13 = G; 
  GC_12 = mdl_complexi * G; 
  GC_11 = -(mdl_complexi * G); 
  GC_10 = -G; 
}

// Routines for printing out parameters
void Parameters_MRSSMQCD_UFO_from_philip::printIndependentParameters()
{
  cout <<  "MRSSMQCD_UFO_from_philip model parameters independent of event kinematics:" <<
      endl;
  cout << setw(20) <<  "mdl_wOp " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_wOp << endl;
  cout << setw(20) <<  "mdl_wOs " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_wOs << endl;
  cout << setw(20) <<  "mdl_wdglu " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_wdglu << endl;
  cout << setw(20) <<  "mdl_WH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WH << endl;
  cout << setw(20) <<  "mdl_WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WW << endl;
  cout << setw(20) <<  "mdl_WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WZ << endl;
  cout << setw(20) <<  "mdl_WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WT << endl;
  cout << setw(20) <<  "mdl_ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymtau << endl;
  cout << setw(20) <<  "mdl_ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymt << endl;
  cout << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  cout << setw(20) <<  "mdl_Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gf << endl;
  cout << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  cout << setw(20) <<  "mdl_MD3 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MD3 << endl;
  cout << setw(20) <<  "mdl_mOp " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_mOp << endl;
  cout << setw(20) <<  "mdl_mOs " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_mOs << endl;
  cout << setw(20) <<  "mdl_mstr " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_mstr << endl;
  cout << setw(20) <<  "mdl_msbr " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_msbr << endl;
  cout << setw(20) <<  "mdl_mscr " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_mscr << endl;
  cout << setw(20) <<  "mdl_mssr " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_mssr << endl;
  cout << setw(20) <<  "mdl_msur " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_msur << endl;
  cout << setw(20) <<  "mdl_msdr " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_msdr << endl;
  cout << setw(20) <<  "mdl_mstl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_mstl << endl;
  cout << setw(20) <<  "mdl_msbl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_msbl << endl;
  cout << setw(20) <<  "mdl_mscl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_mscl << endl;
  cout << setw(20) <<  "mdl_mssl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_mssl << endl;
  cout << setw(20) <<  "mdl_msul " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_msul << endl;
  cout << setw(20) <<  "mdl_msdl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_msdl << endl;
  cout << setw(20) <<  "mdl_MH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MH << endl;
  cout << setw(20) <<  "mdl_MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MZ << endl;
  cout << setw(20) <<  "mdl_Mnu " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Mnu << endl;
  cout << setw(20) <<  "mdl_MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MTA << endl;
  cout << setw(20) <<  "mdl_MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MT << endl;
  cout << setw(20) <<  "mdl_yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yb << endl;
  cout << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2 << endl;
  cout << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4 << endl;
  cout << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2 << endl;
  cout << setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MH__exp__2 << endl;
  cout << setw(20) <<  "mdl_I1a33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I1a33 << endl;
  cout << setw(20) <<  "mdl_I4a33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I4a33 << endl;
  cout << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi << endl;
  cout << setw(20) <<  "mdl_aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_aEW << endl;
  cout << setw(20) <<  "mdl_MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MW << endl;
  cout << setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__aEW << endl;
  cout << setw(20) <<  "mdl_ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ee << endl;
  cout << setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw2 << endl;
  cout << setw(20) <<  "mdl_cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cw << endl;
  cout << setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__sw2 << endl;
  cout << setw(20) <<  "mdl_sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw << endl;
  cout << setw(20) <<  "mdl_g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_g1 << endl;
  cout << setw(20) <<  "mdl_gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gw << endl;
  cout << setw(20) <<  "mdl_vev " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_vev << endl;
  cout << setw(20) <<  "mdl_vev__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_vev__exp__2 << endl;
  cout << setw(20) <<  "mdl_lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_lam << endl;
  cout << setw(20) <<  "mdl_wein " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_wein << endl;
  cout << setw(20) <<  "mdl_yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yt << endl;
  cout << setw(20) <<  "mdl_ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ytau << endl;
  cout << setw(20) <<  "mdl_muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_muH << endl;
  cout << setw(20) <<  "mdl_I2a33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I2a33 << endl;
  cout << setw(20) <<  "mdl_I3a33 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_I3a33 << endl;
  cout << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2 << endl;
  cout << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2 << endl;
}
void Parameters_MRSSMQCD_UFO_from_philip::printIndependentCouplings()
{
  cout <<  "MRSSMQCD_UFO_from_philip model couplings independent of event kinematics:" <<
      endl;

}
void Parameters_MRSSMQCD_UFO_from_philip::printDependentParameters()
{
  cout <<  "MRSSMQCD_UFO_from_philip model parameters dependent on event kinematics:" <<
      endl;
  cout << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS << endl;
  cout << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  cout << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2 << endl;
}
void Parameters_MRSSMQCD_UFO_from_philip::printDependentCouplings()
{
  cout <<  "MRSSMQCD_UFO_from_philip model couplings dependent on event kinematics:" <<
      endl;
  cout << setw(20) <<  "GC_17 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_17 << endl;
  cout << setw(20) <<  "GC_15 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_15 << endl;
  cout << setw(20) <<  "GC_14 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_14 << endl;
  cout << setw(20) <<  "GC_13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_13 << endl;
  cout << setw(20) <<  "GC_12 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_12 << endl;
  cout << setw(20) <<  "GC_11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_11 << endl;
  cout << setw(20) <<  "GC_10 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_10 << endl;
}


