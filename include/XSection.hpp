#ifndef XSECTION_H_
#define XSECTION_H_

#include "Process.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/program_options/variables_map.hpp>
#include "LHAPDF/LHAPDF.h"

#include <array>

// every class with at least one pure virtual function
// is an abstract base class
class XSection {

   public:
      // pure virtual function (abstract function)
      virtual std::array<double, 3> integrate() = 0;

      static void init (Process*, boost::property_tree::ptree const&, boost::program_options::variables_map const&);
      static Process* processID;

   protected:
      static double dS;
      static double dC;
      static boost::property_tree::ptree pt;
      static boost::program_options::variables_map vm;
      static double S_sqrt;
      static double S;
      static double mu_r;
      static double mu_f;
      static double m1;
      static double m2;

      static const LHAPDF::PDF* pdf;
};

#endif /* XSECTION_H_ */
