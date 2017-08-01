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

enum class DipoleType {
       FF, II, FI, IF
};

class CSDipole {

public:
    CSDipole(boost::property_tree::ptree pt, std::string s, DipoleType type, unsigned int i, unsigned int j)
            : Born_(Color_Connected_MEs(i, j, s, pt)), type_(type), i_(i), j_(j) {
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
       std::cout << " with leg " << i << " as the emitter and leg " << j << " as the spectator " << j << '\n';
    };
    double eval(std::vector<Vec4D<double>> const&) const;

private:
    const DipoleType type_;
    const unsigned int i_;
    const unsigned int j_;
    const Color_Connected_MEs Born_;
};
#endif //RSYMSQCD_COLOR_CONNECTED_BORN_ME_H
