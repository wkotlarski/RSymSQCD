#include <chrono>

#include "XSection_Tree.hpp"
#include "XSection_Virt.hpp"
#include "XSection_SC.hpp"
#include "XSection_HnonC.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

using namespace std;

array<double,3> add(array<double,3> x, array<double,3> y) {
   return array<double,3> { x.at(0) + y.at(0), sqrt( pow(x.at(1), 2) + pow(y.at(1), 2) ),
     x.at(2) + y.at(2) };
}

// why do I have to write this?
// why isn't init() enough?
std::array< std::array<double, 2>, 6 > XSection::squark_mass; 
Process *XSection::processID; 
double XSection::muF;
double XSection::muR;
double XSection::m1;
double XSection::m2;
double XSection::prec_virt;
double XSection::prec_sc;
double XSection::prec_hnc;
double XSection::S;
double XSection::S_sqrt;
double XSection::dC;
double XSection::dS;
boost::property_tree::ptree XSection::pt;
const LHAPDF::PDF* XSection::pdf;

int main(int argc, char* argv[]) {
   
   boost::property_tree::ptree pt;
   boost::property_tree::ini_parser::read_ini("config.ini", pt);
   
   array<double,3> temp, xsection_tree, xsection_virt, xsection_SC, xsection_HnonC;
   
   auto start = chrono::steady_clock::now();
   auto t0 = chrono::steady_clock::now();
   
   if( string( argv[2] ) == "LO" ) {
      if ( string(argv[1]) == "pp_OsOs" ) {         
         Process process1("sgluons-gg_OO", pt);
         XSection::init( &process1, pt, 1, 1, 1 );
         XSection_Tree tree;
         temp = tree.integrate();
      
         Process process2("sgluons-qqbar_OO", pt);
         XSection::init( &process2, pt, 1, 1, 1 );
         xsection_tree = add(tree.integrate(), temp);
         cout << xsection_tree.at(0) << endl;
         
      } else if( string(argv[1]) == "pp_suLsuR" ) {
         Process process1("MRSSM-uu_suLsuR", pt);
         XSection::init( &process1, pt, 1, 1, 1 );
         XSection_Tree tree;
         xsection_tree = tree.integrate();
      } else {
         cout << "LO process not implemented\n";
      }
   } else if ( string( argv[2] ) == "NLO" ) {
      if( string(argv[1]) == "pp_suLsuR" ) {
         Process process1("MRSSM-uu_suLsuR", pt);
         XSection::init( &process1, pt, pow(10, -atoi(argv[3])), pow(10, -atoi(argv[4])), pow(10, -atoi(argv[5])) );
         XSection_Tree tree;
         xsection_tree = tree.integrate();
      
         XSection_Virt virt;
         xsection_virt = virt.integrate();
      
         XSection_SC sc;
         xsection_SC = sc.integrate();
      
         XSection_HnonC hc;
         xsection_HnonC = hc.integrate();
      } else {
         cout << "NLO process not implemented\n";
      }
   } else {
      cout << "Second command line argument must be 'LO' or 'NLO'." << endl;
      return 0;
   }
   
   auto t1 = chrono::steady_clock::now();
    
   cout  << "\nBorn part took " 
         << chrono::duration_cast<chrono::seconds>(t1-t0).count() 
         << " s" << endl;
/*
   auto t2 = chrono::steady_clock::now();

   cout  << "\nVirtual part took " 
         << chrono::duration_cast<chrono::seconds>(t2-t1).count() 
         << " s" << endl;


   auto t3 = chrono::steady_clock::now();  
    
   cout << "\nSoft and/or collinear part took " 
         << chrono::duration_cast<chrono::seconds>(t3-t2).count() << " s" << endl;
    

    
   
   cout << "\nHard - non-collinear part took " 
        << chrono::duration_cast<chrono::seconds>(end-t3).count() << " s" << endl;
   cout << endl;

   double total_time = chrono::duration_cast<chrono::seconds>(end-start).count();
 */  
   
   auto end = chrono::steady_clock::now();
   
   cout << "\nRun summary\n";
   cout << "Time: " << chrono::duration_cast<chrono::minutes>(end-start).count()
        << " minutes\n";
   cout << "---------------------------------------------------------------" << endl;
   cout << setprecision(5);
   cout << setw(12) << "tree:" << setw(13) << xsection_tree.at(0) 
         << " +/- " << setprecision(1) << xsection_tree.at(1)
         << " fb ( p-value = " << setw(8) << xsection_tree.at(2) << " )\n";
   if( string( argv[2] ) == "NLO" ) { 
   cout << setprecision(5);
   cout << setw(12) << "virtual:" << setw(13) << xsection_virt.at(0) << " +/- " 
         << setprecision(1) << xsection_virt.at(1) << " fb ( p-value = " 
         << setw(8) << xsection_virt.at(2) << " )\n";

   cout << setprecision(5);
   cout << setw(12) << "real (soft):" << setw(13) << xsection_SC.at(0) << " +/- " << setprecision(1) << xsection_SC.at(1)
         << " fb ( p-value = " << setw(8) << xsection_SC.at(2) << " )\n"; 
   cout << setprecision(5);
   cout << setw(12) << "real (hard):" << setw(13) << xsection_HnonC.at(0) << " +/- " << setprecision(1) << xsection_HnonC.at(1)
         << " fb ( p-value = " << setw(8) << xsection_HnonC.at(2) << " )\n";
   cout << "---------------------------------------------------------------" << endl;
   cout << setprecision(5);
   cout << setw(12) << "sum:" << setw(13)
        << xsection_tree.at(0) + xsection_virt.at(0) + xsection_HnonC.at(0) + xsection_SC.at(0) 
         << " +/- " << setprecision(1) << sqrt(pow(xsection_tree.at(1),2) + 
            pow(xsection_virt.at(1),2) + pow(xsection_HnonC.at(1),2) +
            pow(xsection_SC.at(1),2)) << " fb" << endl;
   }
   cout << '\n';
   
   return 0;
}
