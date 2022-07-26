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
   void init (Process const&, boost::property_tree::ptree const&, boost::program_options::variables_map const&);
   Process processID;

protected:
   double dS;
   double dC;
   boost::property_tree::ptree pt;
   boost::program_options::variables_map vm;
   double S_sqrt;
   double S;
   double mu_r;
   double mu_f;

   const LHAPDF::PDF* pdf;
};

#endif /* XSECTION_H_ */
