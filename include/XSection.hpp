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
static int particle_to_int (Particle p) {
   switch (p) {
      case Particle::d:
         return 1;
      case Particle::dbar:
         return -1;
      case Particle::u:
         return 2;
      case Particle::ubar:
         return -2;
      case Particle::s:
         return 3;
      case Particle::sbar:
         return -3;
      case Particle::c:
         return 4;
      case Particle::cbar:
         return -4;
      case Particle::b:
         return 5;
      case Particle::bbar:
         return -5;
   }
}
};
#endif /* XSECTION_H_ */
