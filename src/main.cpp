#include <chrono>
#include <stdio.h>

#include "XSection_Tree_MRSSM.h"
#include "XSection_Tree_MSSM.h"
#include "XSection_Virt_MRSSM.h"
#include "XSection_Virt_MSSM.h"
#include "XSection_SC.hpp"
#include "XSection_HnonC.h"

using namespace std;

double XSection_Real::dS = 1e-4;
double XSection_Real::dC = dS / 50.;

// why do I have to write this?
// why isn't init() enough?
std::array< std::array<double, 2>, 6 > XSection::squark_mass; 
double XSection::muF;
double XSection::muR;
const LHAPDF::PDF* XSection::pdf_nlo;
const LHAPDF::PDF* XSection::pdf_lo;

int main(int argc, char* argv[]) {
    
    ofstream myfile;
    char file_MRSSM_uu_susu[] =
         "MRSSM_1L_uu_susu_msg=2000GeV.txt";
    myfile.open(file_MRSSM_uu_susu);   

    /* array of squark masses */
    int num = 100;
    double m_sq[num];
    m_sq[0] = 100.;
    for(int i=0;i<=num-2;i++)
    {
        m_sq[i+1] = m_sq[0] + 1900.*(i+2)/num;
    }
    myfile << "m_squark,Tree_XSection,Virt_XSection,Real_XSection,1L_XSection";

auto start = chrono::steady_clock::now();
for(int j=0; j<=num-1;j++)
{
    myfile << m_sq[j] << "  ";    

    // this function initiales parameters in XSection class 
    // at runtime reading values from text file
    XSection::init(m_sq[j]);
    
    // format terminal output
    cout << scientific << setprecision(4);

    // actual calculation
    auto t0 = chrono::steady_clock::now();
    bool MRSSM = true;
    auto t1 = chrono::steady_clock::now(), t2 = chrono::steady_clock::now();
    array<double,3> xsection_tree, xsection_virt;
    if (MRSSM)
    {
        XSection_Tree_MRSSM tree_MRSSM;
        xsection_tree = tree_MRSSM.integrate();
        t1 = chrono::steady_clock::now();
      
        XSection_Virt_MRSSM virt_MRSSM;
        xsection_virt = virt_MRSSM.integrate();
        t2 = chrono::steady_clock::now();
    }
    else // MSSM
    {
        XSection_Tree_MSSM tree_MSSM;
        xsection_tree = tree_MSSM.integrate();
        t1 = chrono::steady_clock::now();
      
        XSection_Virt_MSSM virt_MSSM;
        xsection_virt = virt_MSSM.integrate();
        t2 = chrono::steady_clock::now();
    }

    XSection_SC sc;
    array<double, 3> xsection_SC = sc.integrate();
    auto t3 = chrono::steady_clock::now();  
  
    XSection_HnonC hc;
    array<double, 3> xsection_HnonC = hc.integrate();
    auto t4 = chrono::steady_clock::now();
    
    cout << "\nBorn part took " 
         << chrono::duration_cast<chrono::seconds>(t1-t0).count() << " s" << endl;
    cout << "Result: " << xsection_tree.at(0) << " +/- " << xsection_tree.at(1)
         << " fb ( p-value = " << xsection_tree.at(2) << " )\n";
    
    cout << "\nVirtual part took " 
         << chrono::duration_cast<chrono::seconds>(t2-t1).count() << " s" << endl;
    cout << "Result: " << xsection_virt.at(0) << " +/- " << xsection_virt.at(1)
         << " fb ( p-value = " << xsection_virt.at(2) << " )\n";
    
    cout << "\nSoft and/or collinear part took " 
         << chrono::duration_cast<chrono::seconds>(t3-t2).count() << " s" << endl;
    cout << "Result: " << xsection_SC.at(0) << " +/- " << xsection_SC.at(1)
         << " fb ( p-value = " << xsection_SC.at(2) << " )\n";  
    
    cout << "\nHard - non-collinear part took " 
         << chrono::duration_cast<chrono::seconds>(t4-t3).count() << " s" << endl;
    cout << "Result: " << xsection_HnonC.at(0) << " +/- " << xsection_HnonC.at(1)
         << " fb ( p-value = " << xsection_HnonC.at(2) << " )\n";

    cout << "Total real emission:\n";
    cout << xsection_HnonC.at(0) + xsection_SC.at(0) << " "
         << sqrt( pow(xsection_HnonC.at(1),2) + pow(xsection_SC.at(1),2) )
         <<  '\n' << endl;
    cout << "NLO Result: " << xsection_tree.at(0) + xsection_virt.at(0) 
                            + xsection_HnonC.at(0) + xsection_SC.at(0) << 
                            " fb +/-" << sqrt(pow(xsection_tree.at(1),2) + 
                                              pow(xsection_virt.at(1),2) +
                                              pow(xsection_HnonC.at(1),2) +
                                              pow(xsection_SC.at(1),2)) <<
                            " fb" << endl;
    cout << "Total time needed: "
         << chrono::duration_cast<chrono::seconds>(t4-t0).count()
         << " s\n" << endl;
    myfile << xsection_tree.at(0) << "  " << xsection_virt.at(0) << "   " << xsection_HnonC.at(0)+xsection_SC.at(0) << "    ";
    myfile << xsection_tree.at(0) + xsection_virt.at(0) + xsection_HnonC.at(0) + xsection_SC.at(0) << "\n";
}  
auto end = chrono::steady_clock::now();
myfile.close();                           
cout << "Whole calcluation took " << chrono::duration_cast<chrono::seconds>(end-start).count()/3600. << " hours" << endl;
  return 0;
}
