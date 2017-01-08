#ifndef XSECTION_H_
#define XSECTION_H_

#include <array>

#include <boost/property_tree/ptree.hpp>

#include "LHAPDF/LHAPDF.h"
#include "include/cuba.h"

#include "Process.hpp"

// every class with at least one pure virtual function 
// is an abstract base class
class XSection {

   public:
      // pure virtual function (abstract function)
      virtual std::array<double, 3> integrate() = 0;
      
      static void init (Process*, boost::property_tree::ptree, double, double, double);
      static Process* processID;

   protected:
      static double dS;
      static double dC;
      static boost::property_tree::ptree pt;
      static double S_sqrt;
      static double S;
      static double prec_virt;
      static double prec_sc;
      static double prec_hnc;
      static double mu_r;
      static double mu_f;
      static double m1;
      static double m2;
    
      static const LHAPDF::PDF* pdf;
};

#endif /* XSECTION_H_ */
