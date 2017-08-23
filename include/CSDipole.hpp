// the idea behind this class
// ask the code to create you a dipole for given process, of given type (II,FF, ...) with given legs.


#ifndef RSYMSQCD_COLOR_CONNECTED_BORN_ME_H
#define RSYMSQCD_COLOR_CONNECTED_BORN_ME_H

#include <vector>
#include <string>
#include <iostream>

#include <boost/property_tree/ptree.hpp>

#include "Vec4D.hpp"
#include "Color_Connected_MEs.h"

#include "dilog.hpp"

constexpr double K_q = (7/2. - pi*pi/6.)*CF;
constexpr double K_sq = (4 - pi*pi/6.)*CF;
constexpr double gamma_sq = 2*CF;

enum class DipoleType {
       FF, II, FI, IF
};

enum class DipoleV {
   gQ,k
};

class CSDipole {

public:
   CSDipole(boost::property_tree::ptree pt, std::string s, DipoleType type, unsigned int i, unsigned int j)
            : Born_(Color_Connected_MEs(i, j, s, pt)), type_(type), emit_(i), spec_(j) {
      std::cout << "Initialized Catani-Seymour dipole of type ";
      if (type_ == DipoleType::FF) {
         std::cout << "FF";
      } else if (type_ == DipoleType::IF) {
         std::cout << "IF";
      } else if (type_ == DipoleType::II) {
         std::cout << "II";
      } else if (type_ == DipoleType::FI) {
         std::cout << "FI";
      }
      std::cout << " with leg " << i << " as the emitter and leg " << j << " as the spectator " << '\n';
   };
   double eval_unintegrated_dipole(std::vector<Vec4D<double>> const&) const;
   double eval_integrated_dipole(std::vector<Vec4D<double>> const&) const;
   double eval_P(std::vector<Vec4D<double>> const&, double) const;
   double eval_K(std::vector<Vec4D<double>> const&, double) const;
   unsigned int get_emitter() const { return emit_; }
   unsigned int get_spectator() const { return spec_; }
   const Color_Connected_MEs Born_;

private:
   const DipoleType type_;
   const unsigned int emit_;
   const unsigned int spec_;
   double V_S(double) const;
   double V_NS(double) const;
};
#endif //RSYMSQCD_COLOR_CONNECTED_BORN_ME_H
