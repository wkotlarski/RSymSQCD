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
      static std::array< std::array<double, 2>, 6 > squark_mass;
      //static double muR;
      //static double muF;
      static double mu_r;
      static double mu_f;
      static double m1;
      static double m2;
      
      // SU(3) group factors
      static constexpr double CF = 4/3.;
      static constexpr double CA = 3.;
    
      static const LHAPDF::PDF* pdf;
};

#endif /* XSECTION_H_ */
