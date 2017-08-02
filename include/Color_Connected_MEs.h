#ifndef RSYMSQCD_COLOR_CONNECTED_MES_H
#define RSYMSQCD_COLOR_CONNECTED_MES_H

#include <vector>
#include <string>

#include <boost/property_tree/ptree.hpp>

#include "MRSSM.h"
#include "Vec4D.hpp"
#include "constants.hpp"

class Color_Connected_MEs : public MRSSM {

public:
    Color_Connected_MEs(unsigned int emitter, unsigned int spectator, std::string s,  boost::property_tree::ptree pt) :
            MRSSM(pt),
            mu_r(pt.get<double>("collider setup.mu_r")),
            emitter_(emitter),
            spectator_(spectator)
    {};
   double get_ME2_value(const std::vector<Vec4D<double>>&) const;

private:
   const unsigned int emitter_;
   const unsigned int spectator_;
   const double mu_r;
   double ColorMatrix(const std::string&, const std::string&) const;
   double uu_suLsuR(double, double, double, double, double, double) const;
};


#endif //RSYMSQCD_COLOR_CONNECTED_MES_H
