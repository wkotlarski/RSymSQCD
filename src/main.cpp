#include <chrono>

#include "XSection_Virt.h"
#include "XSection_SC.hpp"
#include "XSection_HnonC.h"

using namespace std;

/*  This is really shitty solution!
 *  Those variables have to be static.
 *  TODO: How to initialize dS and dC in main at run-time?
 */
//std::cout << "initializing LO PDF\n";
const LHAPDF::PDF* XSection::pdf_lo = LHAPDF::mkPDF("MMHT2014lo68cl", 0);
//std::cout << "initializing NLO PDF\n";
const LHAPDF::PDF* XSection::pdf_nlo = LHAPDF::mkPDF("MMHT2014nlo68cl", 0);

double XSection_Real::dS = 2e-3;
double XSection_Real::dC = dS / 50.;

std::array< std::array<double, 2>, 6 > XSection::squark_mass {{
          {1500, 1500},
          {1500, 1500},
          {1500, 1500},
          {1500, 1500},
          {1500, 1500},
          {1500, 1500}
    }};

int main(int argc, char* argv[]) {
    
  // some preparing might go here
    std::complex<double> z1 (1000.7, 0);
    std::cout << setprecision(16) << dilogarithm::dilog(z1) << endl;
  // accual calculation
  auto t1 = chrono::steady_clock::now();

  XSection_Virt v;
  //array<double, 3> xsection_v = v.integrate();
  auto t2 = chrono::steady_clock::now();
//  cout << "\nBorn part took " <<
//        chrono::duration_cast<chrono::seconds>(t2-t1).count() << endl;
//    cout << "Result: " << xsection_v.at(0) << " +/- " << xsection_v.at(1)
//        << " fb ( p-value = " << xsection_v.at(2) << " )\n";

    XSection_HnonC hc;
    array<double, 3> xsection_HnonC = hc.integrate();
    auto t3 = chrono::steady_clock::now();
    cout << "\nHard - non-collinear part took " <<
      chrono::duration_cast<chrono::seconds>(t3-t2).count() << endl;
    cout << "Result: " << xsection_HnonC.at(0) << " +/- " << xsection_HnonC.at(1)
      << " fb ( p-value = " << xsection_HnonC.at(2) << " )\n";

    XSection_SC sc;
    array<double, 3> xsection_SC = sc.integrate();
    auto t4 = chrono::steady_clock::now();
    cout << "\nSoft and/or collinear part took " <<
        chrono::duration_cast<chrono::seconds>(t4-t3).count() << endl;
    cout << "Result: " << xsection_SC.at(0) << " +/- " << xsection_SC.at(1)
        << " fb ( p-value = " << xsection_SC.at(2) << " )\n";

    cout << "Total real emission:\n";
    cout    << xsection_HnonC.at(0) + xsection_SC.at(0) << " "
            << sqrt( pow(xsection_HnonC.at(1),2) + pow(xsection_SC.at(1),2) )
            <<  '\n';

  
  return 1;
}
