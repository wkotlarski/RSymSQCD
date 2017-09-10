#ifndef RSYMSQCD_SM_H
#define RSYMSQCD_SM_H

#include "IMatrixElements.h"

class SM : public IMatrixElements {

public:
   SM(boost::property_tree::ptree const& ptree);

   double BornME(std::vector<Particle>, double, double) const noexcept;
   double BornCCME(std::vector<Particle>, int, int, EpsOrd, std::vector<Vec4D<double>> const&) const noexcept;
   double VirtualME(std::vector<Particle>, EpsOrd, double, double) const noexcept;
   double RealME(std::vector<Particle>, std::vector<Vec4D<double>> const&) const noexcept;

private:
   const double MB_;
   const double Alfa_;

   double eebar_bbbar_B(double, double) const noexcept;
   double eebar_bbbar_CCB(int, int, EpsOrd, double, double) const noexcept;
   double eebar_bbbar_V(EpsOrd, double, double) const noexcept;
   double eebar_bbbar_R (std::vector<Vec4D<double>> const&) const noexcept;
};

#endif //RSYMSQCD_SM_H
