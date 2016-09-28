#include <chrono>
#include <stdio.h>

#include "XSection_Tree.h"
#include "XSection_Virt.h"
#include "XSection_SC.hpp"
#include "XSection_HnonC.h"

using namespace std;

double XSection_Real::dS = 1e-4;
double XSection_Real::dC = 1e-5;

// why do I have to write this?
// why isn't init() enough?
std::array< std::array<double, 2>, 6 > XSection::squark_mass; 
double XSection::gluino_mass; 
double XSection::sgluon_mass; 
Process *XSection::processID; 
double XSection::muF;
double XSection::muR;
const LHAPDF::PDF* XSection::pdf_nlo;
const LHAPDF::PDF* XSection::pdf_lo;

int main(int argc, char* argv[]) {

   auto start = chrono::steady_clock::now();

    double m_squark = 1500.;//M_min + i/(num-1.)*(M_maxSq - M_min);
    double m_gluino = 1000.;//M_min + j/(num-1.)*(M_maxGlu - M_min);
    Process process("MRSSM,uu_suLsuR");

    XSection::init(m_squark, m_gluino, 5000., &process);      //squark mass, gluino mass, sgluon mass, string which specifies the matrix element
    
   auto t0 = chrono::steady_clock::now();
   array<double,3> xsection_tree, xsection_virt;
   XSection_Tree tree;
   xsection_tree = tree.integrate();
   auto t1 = chrono::steady_clock::now();
    
   cout  << "\nBorn part took " 
         << chrono::duration_cast<chrono::seconds>(t1-t0).count() 
         << " s" << endl;
   
   cout  << "Result: " << xsection_tree.at(0) 
         << " +/- " << xsection_tree.at(1)
         << " fb ( p-value = " << xsection_tree.at(2) << " )\n";

   XSection_Virt virt;
   xsection_virt = virt.integrate();
   auto t2 = chrono::steady_clock::now();

   cout  << "\nVirtual part took " 
         << chrono::duration_cast<chrono::seconds>(t2-t1).count() 
         << " s" << endl;
   cout  << "Result: " << xsection_virt.at(0) << " +/- " 
         << xsection_virt.at(1) << " fb ( p-value = " 
         << xsection_virt.at(2) << " )\n";

   XSection_SC sc;
   array<double, 3> xsection_SC = sc.integrate();
   auto t3 = chrono::steady_clock::now();  
    
   cout << "\nSoft and/or collinear part took " 
         << chrono::duration_cast<chrono::seconds>(t3-t2).count() << " s" << endl;
    cout << "Result: " << xsection_SC.at(0) << " +/- " << xsection_SC.at(1)
         << " fb ( p-value = " << xsection_SC.at(2) << " )\n";  
    
    XSection_HnonC hc;
    array<double, 3> xsection_HnonC = hc.integrate();
    auto t4 = chrono::steady_clock::now();
    
    cout << "\nHard - non-collinear part took " 
         << chrono::duration_cast<chrono::seconds>(t4-t3).count() << " s" << endl;
    cout << "Result: " << xsection_HnonC.at(0) << " +/- " << xsection_HnonC.at(1)
         << " fb ( p-value = " << xsection_HnonC.at(2) << " )\n";

    cout << endl;
    cout << "Total real emission:\n";
    cout << xsection_HnonC.at(0) + xsection_SC.at(0) << " "
         << sqrt( pow(xsection_HnonC.at(1),2) + pow(xsection_SC.at(1),2) )
         <<  '\n' << endl;
   
    cout << "NLO Result: " 
         << xsection_tree.at(0) + xsection_virt.at(0) 
            + xsection_HnonC.at(0) + xsection_SC.at(0) 
         << " fb +/- " << sqrt(pow(xsection_tree.at(1),2) + 
            pow(xsection_virt.at(1),2) + pow(xsection_HnonC.at(1),2) +
            pow(xsection_SC.at(1),2)) << " fb" << endl;
    
    cout << "Total time needed: "
         << chrono::duration_cast<chrono::seconds>(t4-t0).count()
         << " s\n" << endl;

    auto end = chrono::steady_clock::now();

   double total_time = chrono::duration_cast<chrono::seconds>(end-start).count();

   cout << "Whole calcluation took ";
   if ( total_time < 60 ) {
      cout << total_time << " seconds\n";
   }
   else if (total_time >= 60 && total_time < 3600 ) {
      cout << total_time/60.0 << " minutes\n";
   }
   else if (total_time >= 3600 ) {
      cout << total_time/3600.0 << " hours\n";
   }

   return 0;
}
