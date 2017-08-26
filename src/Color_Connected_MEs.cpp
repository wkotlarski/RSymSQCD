#include <cmath>

#include "../include/Color_Connected_MEs.h"
#include "ColorFull/Core/Tree_level_gluon_basis.h"

double Color_Connected_MEs::get_ME2_value(
      unsigned int emitter,
      unsigned int spectator,
      std::vector<Vec4D<double>> const& p
   ) const {

   auto [k12, k23, k24, k13, k14, k34] = std::tuple(
           p[0]*p[1], p[1]*p[2], p[1]*p[3], p[0]*p[2], p[0]*p[3], p[2]*p[3]
   );
   return eebar_ttbar(emitter, spectator, k12, k23, k24, k13, k14, k34);
}

double Color_Connected_MEs::ColorMatrix( int emitter,  int spectator, std::string const& col_str1, std::string const& col_str2 ) const {
   ColorFull::Col_amp Ca1(col_str1);
   ColorFull::Col_amp Ca2(col_str2);

   ColorFull::Col_functions Col_fun;
   // FormCalc numbers particles from 1, hence i_+1
   if(emitter < 0 && spectator < 0) {
      return Col_fun.double_num(
              Col_fun.scalar_product(Ca1, Ca2)
      );
   } else {
      return Col_fun.double_num(
              Col_fun.scalar_product(
                      // emit gluon from emitter and attach it to the spectator
                      Col_fun.emit_gluon(Ca1, emitter+1, 77), Col_fun.emit_gluon(Ca2, spectator+1, 77)
              )
      );
   }
}

double Color_Connected_MEs::uu_suLsuR(
      unsigned int emitter, unsigned int spectator,
      double k12, double k23, double k24, double k13, double k14, double k34
   ) const {

   const double Alfas2 = pow (pdf->alphasQ( mu_r ), 2);
   double m1 = MassSuL;
   double m2 = MassSuR;
   double temp = 1/9.*(8*Alfas2*(2*k13*k23 - k12*pow(m1,2))*pow(pi,2)*(-3*(ColorMatrix(emitter, spectator, "[{1,3}{2,4}]","[{1,4}{2,3}]") + ColorMatrix(emitter, spectator, "[{1,4}{2,3}]","[{1,3}{2,4}]"))*
      (2*(2*k13 + 2*k14 - pow(m1,2) - pow(m2,2))*pow(MassGlu,2) + 2*pow(MassGlu,4) + pow(-2*k13 + pow(m1,2),2) +
      pow(-2*k14 + pow(m2,2),2)) + ColorMatrix(emitter, spectator, "[{1,3}{2,4}]","[{1,3}{2,4}]")*
      (2*(18*k13 + 2*k14 - 9*pow(m1,2) - pow(m2,2))*pow(MassGlu,2) + 10*pow(MassGlu,4) + 9*pow(-2*k13 + pow(m1,2),2) +
      pow(-2*k14 + pow(m2,2),2)) + ColorMatrix(emitter, spectator, "[{1,4}{2,3}]","[{1,4}{2,3}]")*
      (2*(2*k13 + 18*k14 - pow(m1,2) - 9*pow(m2,2))*pow(MassGlu,2) + 10*pow(MassGlu,4) + pow(-2*k13 + pow(m1,2),2) +
      9*pow(-2*k14 + pow(m2,2),2)))*pow(2*k13 - pow(m1,2) + pow(MassGlu,2),-2)*pow(2*k14 - pow(m2,2) + pow(MassGlu,2),-2))/
      9.;
   return temp;
}

double Color_Connected_MEs::eebar_ttbar(
        unsigned int emitter, unsigned int spectator,
        double k12, double k23, double k24, double k13, double k14, double k34
) const {
   const double MB2 = MassSuL*MassSuL;
   const double Alfa2 = pow (137., -2);
   const double T = MB2 - 2*k13;
   const double U = MB2 - 2*k23;
   const double S = 2*k12;
   double temp = 4*(8*Alfa2*ColorMatrix(emitter, spectator, "[{3,4}]","[{3,4}]")*pow(pi,2)*pow(S,-2)*(2*MB2*(S - T - U) + 2*pow(MB2,2) + pow(T,2) + pow(U,2)))/9.;
   return temp;
}


