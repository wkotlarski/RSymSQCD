#include <iostream>

#include "CS_FF_Dipoles.hpp"
#include "constants.hpp"

//using namespace ColorFull;
using namespace std;

CS_FF_Dipole::CS_FF_Dipole(
      double (*me)(valarray<valarray<double>> const&),
      int j, int i, int k,
      valarray<valarray<double>> const& p) : CS_Dipole(nullptr, i, j, k, p) {

   assert( i != j && i != k && j !=k);
   q_ = p_[i_] + p_[j_] + p_[k_];
   ptilde_ = mom_shuffle();
   born_me2_ = me;
}

inline double CS_FF_Dipole::yijk() {
   return dot(p_[i_], p_[j_])/(dot(p_[i_], p_[j_]) + dot(p_[i_], p_[k_]) + dot(p_[j_], p_[k_]));
}

inline double CS_FF_Dipole::z( int l ) {
   return dot(p_[l], p_[k_])/(dot(p_[i_], p_[k_]) + dot(p_[j_], p_[k_]));
}

// @todo the only place not symmetric in ij
inline double CS_FF_Dipole::Vijk() {
   return 8 * pi * Alfas * CF * (
      2/(1 - z(i_)*(1-yijk()))
      - beta(ptilde_[i_], ptilde_[k_])/beta(p_[i_]+p_[j_], p_[k_]) * (1 + z(i_) + mi*mi/dot(p_[i_], p_[j_]))
   );
}

inline double CS_FF_Dipole::Born(valarray<valarray<double>> const& p) {
   double k12 = dot(ptilde_[0], ptilde_[1]);
   double k23 = dot(ptilde_[1], ptilde_[2]);
   double k24 = dot(ptilde_[1], ptilde_[3]);
   double k13 = dot(ptilde_[0], ptilde_[2]);
   double k14 = dot(ptilde_[0], ptilde_[3]);
   double k34 = dot(ptilde_[2], ptilde_[3]);
   const string cos = "[{3, 4}]";
   return ColorMatrix(cos,cos) * 8 * Alfa2 * pi_sqr * (mi*mi*k12 + k14*k23 + k13*k24)/(9.*pow(k12,2));
}

// eq. 5.9 from CS'97: combine particles i & j -> ij using k as the spectator
valarray<valarray<double>> CS_FF_Dipole::mom_shuffle() {
   const double q2 = dot(q_, q_);
   auto mij {mi};
   valarray<valarray<double>> temp (4);
   temp[0] = p_[0];
   temp[1] = p_[1];
   temp[k_] =
      0.5*(q2 + mk*mk - mij*mij)/q2 * q_
      + sqrt(lambda(q2, mij*mij, mk*mk)/lambda(q2, dot(p_[i_] + p_[j_], p_[i_] + p_[j_]), mk*mk))
         * (p_[k_] - dot(q_, p_[k_])/q2 * q_);
   temp[i_] = q_ - temp[k_];
   return temp;
}

double CS_FF_Dipole::unsu(const valarray<valarray<double>>& p) {
   double k12 = dot(p_[0], p_[1]);
   double k23 = dot(p_[1], p_[2]);
   double k25 = dot(p_[1], p_[4]);
   double k24 = dot(p_[1], p_[3]);
   double k35 = dot(p_[2], p_[4]);
   double k34 = dot(p_[2], p_[3]);
   double k45 = dot(p_[3], p_[4]);
   double k13 = dot(p_[0], p_[2]);
   double k14 = dot(p_[0], p_[3]);
   double k15 = dot(p_[0], p_[4]);
   double MB2 = pow(mi,2);
   double unsu2 = (-512*Alfa2*Alfas*((-2*k14*k24*MB2 - k12*k35*MB2 + k12*MB2*(k12 + MB2) + (-2*k14*k24 + k12*(k14 + k24))*(-k13 - k23 + k34 + MB2))/(4.*(k35*k35)) -
       (k15*k15*k23 - k15*(k23*k23) - 2*(k12*k12)*(k13 + k23) - k13*k13*k25 - k15*k23*k25 - 2*k15*k23*k35 - 2*k14*k23*MB2 + k13*(k15*(k23 - k25) + k25*(k23 + k25) + 2*(-2*k23 - k25)*k35 - 2*k24*MB2) +
          k12*(4*(k13*k13) + 4*(k23*k23) + k23*(2*k15 + 3*k25 - 2*k35) - (k15 + k25 - 2*k35)*k35 + k13*(3*k15 - 2*(-2*k23 - k25 + k35)) + (-k15 - k25 + 2*(k34 + k35))*MB2))/(4.*k35*k45) +
       (k15*(-(k13*k23) + k23*k23) + k13*k13*k25 - k13*k23*k25 + 2*k13*k23*k35 + k12*k12*MB2 - 2*k13*k23*MB2 + k12*(MB2*MB2 + k15*(k13 + MB2) + k25*(k23 + MB2) - k35*(k13 + k23 + MB2)))/(4.*(k45*k45)))*pow(pi,3))/
   (9.*(k12*k12));
   assert(unsu2 > 0);
   return unsu2;
}

// eq. 5.2 in CS'2002; color structure is inside of Born
double CS_FF_Dipole::get_dipoles_value(const valarray<valarray<double>>& p) {
   return - Born(ptilde_) * Vijk() / (dot(p_[i_] + p_[j_], p_[i_] + p_[j_]) - mi*mi);
}

double CS_FF_Dipole::ColorMatrix(const string& col_str1, const string& col_str2 ) {
   ColorFull::Col_amp Ca1(col_str1);
   ColorFull::Col_amp Ca2(col_str2);

   ColorFull::Col_functions Col_fun;
   // FormCalc numbers particles from 1, hence i_+1
   double temp = Col_fun.double_num(
      Col_fun.scalar_product(
         // emit gluon from emitter i_ and attach it to the spectator k_
         Col_fun.emit_gluon(Ca1, i_+1, 77), Col_fun.emit_gluon(Ca2, k_+1, 77)
      )
   );
   return temp/CF;
}

