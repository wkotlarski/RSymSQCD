#include <chrono>

#include "XSection_SC.h"
#include "XSection_HnonC.h"

using namespace std;

/*
 *  This is really shitty solution!
 *  Those variables have to be static.
 *  TODO: How to initialize dS and dC in main at run-time?
 */
const LHAPDF::PDF* XSection::pdf_lo = LHAPDF::mkPDF("MMHT2014lo68cl", 0);
const LHAPDF::PDF* XSection::pdf_nlo = LHAPDF::mkPDF("MMHT2014nlo68cl", 0);
double XSection_HnonC::dS = 1e-3;
double XSection_HnonC::dC = 1e-4;

int main(int argc, char* argv[]) {

  // some preparing might go here

  // accual calculation
  auto t1 = chrono::steady_clock::now();

  XSection_HnonC hc;
  array<double, 3> xsection_HnonC = hc.integrate();
  auto t2 = chrono::steady_clock::now();
  cout << "\nHard - non-collinear part took " <<
      chrono::duration_cast<chrono::seconds>(t1-t2).count() << endl;
  cout << "Result: " << xsection_HnonC.at(0) << " +/- " << xsection_HnonC.at(1)
      << " fb ( p-value = " << xsection_HnonC.at(2) << " )\n";

  XSection_SC sc;
  array<double, 3> xsection_SC = sc.integrate();
  auto t3 = chrono::steady_clock::now();
    cout << "\nSoft and/or collinear part took " <<
        chrono::duration_cast<chrono::seconds>(t3-t2).count() << endl;
    cout << "Result: " << xsection_SC.at(0) << " +/- " << xsection_SC.at(1)
        << " fb ( p-value = " << xsection_SC.at(2) << " )\n";

  cout << "\n";
  return 1;
}
