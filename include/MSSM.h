#ifndef RSYMSQCD_MSSM_H
#define RSYMSQCD_MSSM_H

#include "IMatrixElements.h"

class MSSM : public IMatrixElements {

public:
   MSSM(boost::property_tree::ptree const&);

   double BornME(std::vector<Particle>, double, double) const noexcept {return 0.;}
   double BornCCME(std::vector<Particle>, int, int, EpsOrd, std::vector<Vec4D<double>> const&) const noexcept {return 0.;}
   double VirtualME(std::vector<Particle>, EpsOrd, double, double) const noexcept {return 0.;}
   double RealME(std::vector<Particle>, std::vector<Vec4D<double>> const&) const noexcept {return 0.;};

private:
};

#endif //RSYMSQCD_MSSM_H
