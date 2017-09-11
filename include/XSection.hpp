#ifndef XSECTION_H_
#define XSECTION_H_

#include <array>
#include <vector>

#include <boost/property_tree/ptree.hpp>
#include <boost/program_options/variables_map.hpp>

#include "LHAPDF/LHAPDF.h"
#include "include/cuba.h"

#include "IMatrixElements.h"

class XSection {

public:
   static void init (boost::property_tree::ptree, boost::program_options::variables_map);
   virtual std::array<double, 3> integrate() = 0;

protected:
   static std::vector<std::vector<Particle>> particles;
   static IMatrixElements* model;

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
