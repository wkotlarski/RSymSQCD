#include <cmath>

#include "../include/Color_Connected_MEs.h"
#include "ColorFull/Core/Tree_level_gluon_basis.h"

double Color_Connected_MEs::get_ME2_value(const std::vector<Vec4D<double>> & v) const {
   const double Alfas2 = pow(pdf->alphasQ( mu_r ), 2);
   double k12 = v[0]*v[1];
   double k23 = v[1]*v[2];
   double k24 = v[1]*v[3];
   double k13 = v[0]*v[2];
   double k14 = v[0]*v[3];
   double k34 = v[2]*v[3];
   double m1 = MassSuL;
   double m2 = MassSuR;
   double temp = 1/9.*(8*Alfas2*(2*k13*k23 - k12*pow(m1,2))*pow(pi,2)*(-3*(ColorMatrix("[{1,3}{2,4}]","[{1,4}{2,3}]") + ColorMatrix("[{1,4}{2,3}]","[{1,3}{2,4}]"))*
      (2*(2*k13 + 2*k14 - pow(m1,2) - pow(m2,2))*pow(MassGlu,2) + 2*pow(MassGlu,4) + pow(-2*k13 + pow(m1,2),2) +
      pow(-2*k14 + pow(m2,2),2)) + ColorMatrix("[{1,3}{2,4}]","[{1,3}{2,4}]")*
      (2*(18*k13 + 2*k14 - 9*pow(m1,2) - pow(m2,2))*pow(MassGlu,2) + 10*pow(MassGlu,4) + 9*pow(-2*k13 + pow(m1,2),2) +
      pow(-2*k14 + pow(m2,2),2)) + ColorMatrix("[{1,4}{2,3}]","[{1,4}{2,3}]")*
      (2*(2*k13 + 18*k14 - pow(m1,2) - 9*pow(m2,2))*pow(MassGlu,2) + 10*pow(MassGlu,4) + pow(-2*k13 + pow(m1,2),2) +
      9*pow(-2*k14 + pow(m2,2),2)))*pow(2*k13 - pow(m1,2) + pow(MassGlu,2),-2)*pow(2*k14 - pow(m2,2) + pow(MassGlu,2),-2))/
      9.;
   return temp * sqrt(Alfas2);
}
double Color_Connected_MEs::ColorMatrix(std::string const& col_str1, std::string const& col_str2 ) const {
   ColorFull::Col_amp Ca1(col_str1);
   ColorFull::Col_amp Ca2(col_str2);

   ColorFull::Col_functions Col_fun;
   // FormCalc numbers particles from 1, hence i_+1
   double temp = Col_fun.double_num(
           Col_fun.scalar_product(
                   // emit gluon from emitter and attach it to the spectator
                   Col_fun.emit_gluon(Ca1, emitter_+1, 77), Col_fun.emit_gluon(Ca2, spectator_+1, 77)
           )
   );
   return temp/CF;
}
