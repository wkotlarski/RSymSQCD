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

struct XSectionParameters {
   double sqrtS;
   double muR;
   double muF;
   LHAPDF::PDF* pdf;
};

class XSection {
public:
   XSection(XSectionParameters const& parameters, double m1, double m2, std::vector<std::array<int, 3>> const& flav)
      : m1_(m1), m2_(m2), sqrtS_(parameters.sqrtS), muR_(parameters.muR), muF_(parameters.muF), flav_(flav), pdf_(parameters.pdf) {};

   void init (boost::program_options::variables_map const&);

protected:
   boost::program_options::variables_map vm;

   const double sqrtS_;
   const double m1_;
   const double m2_;
   const double muR_;
   const double muF_;
   const std::vector<std::array<int, 3>> flav_ {};
   const LHAPDF::PDF* const pdf_;
};

#endif // XSECTION_H_
