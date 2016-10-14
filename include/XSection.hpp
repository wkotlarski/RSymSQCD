#ifndef XSECTION_H_
#define XSECTION_H_

#include <boost/property_tree/ptree.hpp>
#include "LHAPDF/LHAPDF.h"
#include "include/cuba.h"
#include <string>

#include "constants.hpp"
#include "mathematica_wrapper.hpp"
#include "Process.hpp"

class XSection {

  public:
    //Xsection(Process);
    virtual std::array<double, 3> integrate() = 0;
    static void init (Process*, boost::property_tree::ptree, double, double, double);
    static Process* processID;

  protected:
    
   static double dS;
   static double dC;
    static boost::property_tree::ptree ptr;
    static double S_sqrt;
    static double S;
    static double prec_virt;
    static double prec_sc;
    static double prec_hnc;
    static std::array< std::array<double, 2>, 6 > squark_mass;
    static double muR;
    static double muF;
    static double gluino_mass;

    static const LHAPDF::PDF* pdf_lo;
    static const LHAPDF::PDF* pdf_nlo;
    /*
     *  TODO: add random generation of phase space point which can be used to check
     *  IR finiteness
     */
};

#endif /* XSECTION_H_ */
