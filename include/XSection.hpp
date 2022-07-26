#ifndef XSECTION_H_
#define XSECTION_H_

#include "Process.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/program_options/variables_map.hpp>

#include <array>
#include <vector>

namespace LHAPDF {
   class PDF;
}

// every class with at least one pure virtual function
// is an abstract base class
class XSection {
public:
   XSection(double m1, double m2, double muR, double muF, std::vector<std::array<int, 3>> const& flav, const LHAPDF::PDF* const pdf)
      : m1_(m1), m2_(m2), muR_(muR), muF_(muF), flav_(flav), pdf_(pdf) {};

   void init (boost::property_tree::ptree const&, boost::program_options::variables_map const&);

protected:
   boost::property_tree::ptree pt;
   boost::program_options::variables_map vm;

   double sqrtS_;
   const double m1_;
   const double m2_;
   const double muR_;
   const double muF_;
   const std::vector<std::array<int, 3>> flav_ {};
   const LHAPDF::PDF* const pdf_;
};

#endif // XSECTION_H_
