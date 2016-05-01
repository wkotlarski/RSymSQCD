#include <iostream>

#include "XSection_SC.h"
#include "XSection_HnonC.h"

using namespace std;

/*
 *  This is really shitty solution!
 *  Those variables have to be static.
 *  TODO: How to initialize dS and dC in main at run-time?
 */
const LHAPDF::PDF* XSection::pdf = LHAPDF::mkPDF("MMHT2014nlo68cl", 0);
double XSection_HnonC::dS = 1e-3;
double XSection_HnonC::dC = 1e-4;

int main(int argc, char* argv[]) {

  XSection_HnonC hc;
  array<double, 3> xsection_HnonC = hc.integrate();
  cout << xsection_HnonC.at(0) << " +/- " << xsection_HnonC.at(1)
      << " fb ( p-value = " << xsection_HnonC.at(2) << " )\n";
  return 1;
}
