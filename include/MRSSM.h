#ifndef RSYMSQCD_MRSSM_H
#define RSYMSQCD_MRSSM_H

#include <string>
#include <cmath>
#include <iostream>


#include "IMatrixElements.h"

class MRSSM : public IMatrixElements {

public:
   MRSSM(boost::property_tree::ptree const& pt) :
//      pdf (LHAPDF::mkPDF( pt.get<std::string>("collider setup.pdf") , 0)),
      MassTop (pt.get<double>("masses.top")),
      MassGlu (pt.get<double>("masses.gluino")),
      MasssigmaO (pt.get<double>("masses.pseudoscalar_sgluon")),
      MassphiO (sqrt( pow(MasssigmaO,2) + 4.0 * pow(MassGlu, 2))),
      MassSuL (pt.get<double>("masses.suL")),
      MassSuR (pt.get<double>("masses.suR")),
      MassSdL (pt.get<double>("masses.sdL")),
      MassSdR (pt.get<double>("masses.sdR")),
      MassSsL (pt.get<double>("masses.ssL")),
      MassSsR (pt.get<double>("masses.ssR")),
      MassScL (pt.get<double>("masses.scL")),
      MassScR (pt.get<double>("masses.scR")),
      MassSbL (pt.get<double>("masses.sbL")),
      MassSbR (pt.get<double>("masses.sbR")),
      MassStL (pt.get<double>("masses.stL")),
      MassStR (pt.get<double>("masses.stR")) {};

   double BornME(std::vector<Particle>, EpsOrd, double, double) const noexcept {return 0.;}
   double BornCCME(std::vector<Particle>, int, int, EpsOrd, std::vector<Vec4D<double>> const&) const noexcept {return 0.;}
   double VirtualME(std::vector<Particle>, EpsOrd, double, double) const noexcept {return 0.;}
   double RealME(std::vector<Particle>, std::vector<Vec4D<double>> const&) const noexcept {return 0.;}

protected:
   const double MasssigmaO;
   const double MassphiO;
   const double MassGlu;
   const double MassTop;
   const double MassSuL;
   const double MassSuR;
   const double MassSdL;
   const double MassSdR;
   const double MassSsL;
   const double MassSsR;
   const double MassScL;
   const double MassScR;
   const double MassSbL;
   const double MassSbR;
   const double MassStL;
   const double MassStR;
};


#endif //RSYMSQCD_MRSSM_H
