#include <iostream>

#include "../include/CS_FF_Dipoles.hpp"
#include "../include/constants.hpp"

//using namespace ColorFull;
using namespace std;

CS_FF_Dipole::CS_FF_Dipole(
      double (*me)(vector<Vec4D<double>> const&),
      int i, int j, int k,
      vector<Vec4D<double>> const& p) : CS_Dipole {nullptr, i, j, k, p} {

   assert( i != j && i != k && j !=k);
   q_ = p_[i_] + p_[j_] + p_[k_];
   ptilde_ = mom_shuffle();
   born_me2_ = me;
}

inline double CS_FF_Dipole::yijk() {
   return p_[i_]*p_[j_]/( p_[i_]*p_[j_] + p_[i_]*p_[k_] + p_[j_]*p_[k_]);
}

inline double CS_FF_Dipole::z( int l ) {
   return p_[l]*p_[k_]/(p_[i_]*p_[k_] + p_[j_]*p_[k_]);
}

// @todo the only place not symmetric in ij
inline double CS_FF_Dipole::Vijk() {
   return 8 * pi * Alfas * CF * (
      2/(1 - z(j_)*(1-yijk()))
      - beta(ptilde_[j_], ptilde_[k_])/beta(p_[i_]+p_[j_], p_[k_]) * (1 + z(j_) + mj*mj/(p_[i_]*p_[j_]))
   );
}

inline double CS_FF_Dipole::Born(vector<Vec4D<double>> const& v) {
   double k12 = v[0]*v[1];
   double k23 = v[1]*v[2];
   double k24 = v[1]*v[3];
   double k13 = v[0]*v[2];
   double k14 = v[0]*v[3];
   double k34 = v[2]*v[3];
   double m1 = 1.5e+3;
   double m2 = 1.5e+3;
   double Alfas2 = pow(0.1184, 2);
   double MassGlu = 2e+3;
    return 1/9.*(8*Alfas2*(2*k13*k23 - k12*pow(m1,2))*pow(pi,2)*(-3*(ColorMatrix("[{1,3}{2,4}]","[{1,4}{2,3}]") + ColorMatrix("[{1,4}{2,3}]","[{1,3}{2,4}]"))*
                                                            (2*(2*k13 + 2*k14 - pow(m1,2) - pow(m2,2))*pow(MassGlu,2) + 2*pow(MassGlu,4) + pow(-2*k13 + pow(m1,2),2) +
                                                             pow(-2*k14 + pow(m2,2),2)) + ColorMatrix("[{1,3}{2,4}]","[{1,3}{2,4}]")*
                                                                                          (2*(18*k13 + 2*k14 - 9*pow(m1,2) - pow(m2,2))*pow(MassGlu,2) + 10*pow(MassGlu,4) + 9*pow(-2*k13 + pow(m1,2),2) +
                                                                                           pow(-2*k14 + pow(m2,2),2)) + ColorMatrix("[{1,4}{2,3}]","[{1,4}{2,3}]")*
                                                                                                                        (2*(2*k13 + 18*k14 - pow(m1,2) - 9*pow(m2,2))*pow(MassGlu,2) + 10*pow(MassGlu,4) + pow(-2*k13 + pow(m1,2),2) +
                                                                                                                         9*pow(-2*k14 + pow(m2,2),2)))*pow(2*k13 - pow(m1,2) + pow(MassGlu,2),-2)*pow(2*k14 - pow(m2,2) + pow(MassGlu,2),-2))/
           9.;
}

// eq. 5.9 from CS'97: combine particles i & j -> ij using k as the spectator
vector<Vec4D<double>> CS_FF_Dipole::mom_shuffle() {
   const double q2 = q_*q_;
   auto mij {mj};
   vector<Vec4D<double>> temp (4);
   temp[0] = p_[0];
   temp[1] = p_[1];
   temp[k_] =
      0.5*(q2 + mk*mk - mij*mij)/q2 * q_
      + sqrt(lambda(q2, mij*mij, mk*mk)/lambda(q2, (p_[i_] + p_[j_])*(p_[i_] + p_[j_]), mk*mk))
         * (p_[k_] - (q_*p_[k_])/q2 * q_);
   temp[j_] = q_ - temp[k_];
   return temp;
}

double CS_FF_Dipole::unsu(const vector<Vec4D<double>>& p) {
   double k12 = p_[0]*p_[1];
   double k23 = p_[1]*p_[2];
   double k25 = p_[1]*p_[4];
   double k24 = p_[1]*p_[3];
   double k35 = p_[2]*p_[4];
   double k34 = p_[2]*p_[3];
   double k45 = p_[3]*p_[4];
   double k13 = p_[0]*p_[2];
   double k14 = p_[0]*p_[3];
   double k15 = p_[0]*p_[4];
   double MB2 = 0.;//pow(mj,2);
   double m1 = 0.;//mj;
   // e+e- -> bbg
   /*
   double unsu2 = (-512*Alfa2*Alfas*((-2*k14*k24*MB2 - k12*k35*MB2 + k12*MB2*(k12 + MB2) + (-2*k14*k24 + k12*(k14 + k24))*(-k13 - k23 + k34 + MB2))/(4.*(k35*k35)) -
       (k15*k15*k23 - k15*(k23*k23) - 2*(k12*k12)*(k13 + k23) - k13*k13*k25 - k15*k23*k25 - 2*k15*k23*k35 - 2*k14*k23*MB2 + k13*(k15*(k23 - k25) + k25*(k23 + k25) + 2*(-2*k23 - k25)*k35 - 2*k24*MB2) +
          k12*(4*(k13*k13) + 4*(k23*k23) + k23*(2*k15 + 3*k25 - 2*k35) - (k15 + k25 - 2*k35)*k35 + k13*(3*k15 - 2*(-2*k23 - k25 + k35)) + (-k15 - k25 + 2*(k34 + k35))*MB2))/(4.*k35*k45) +
       (k15*(-(k13*k23) + k23*k23) + k13*k13*k25 - k13*k23*k25 + 2*k13*k23*k35 + k12*k12*MB2 - 2*k13*k23*MB2 + k12*(MB2*MB2 + k15*(k13 + MB2) + k25*(k23 + MB2) - k35*(k13 + k23 + MB2)))/(4.*(k45*k45)))*pow(pi,3))/
   (9.*(k12*k12));
    */
   // qq -> e+e-
   /*
   double unsu2 = (-256*Alfa2*Alfas*((2*(k13*k13)*k23 + 3*k13*k14*k23 + k14*k14*k23 - k14*(k23*k23) - k13*k13*k24 - k13*k14*k24 + 3*k13*k23*k24 + 3*k14*k23*k24 - k13*(k24*k24) - 2*k13*k23*k34 - k14*k23*k34 - k13*k24*k34 -
                       2*k14*k24*k34 - 2*k13*k23*(m1*m1) - k14*k23*(m1*m1) + k13*k24*(m1*m1) + 2*k15*(-k23 - 2*k25)*(m1*m1) + k12*k12*(2*k23 + m1*m1) +
                       k12*(-3*k14*k23 - 4*k23*k24 + 2*k23*k34 + 2*k24*k34 + 3*k15*(m1*m1) + k23*(m1*m1) - k25*(m1*m1) + k35*(m1*m1) + k13*(-4*k23 + k24 - m1*m1)))/(4.*(k15*k15)) +
                      (-(k14*k14*k23) - k14*(k23*k23) - k13*k13*k24 - k14*k23*k24 - k14*k23*k34 - 2*k14*k24*k34 + k14*k23*(m1*m1) - 4*k15*k25*(m1*m1) + k12*k12*(2*k13 + m1*m1) -
                       k13*(-2*(k23*k23) + k14*(-3*k23 - 3*k24) - 3*k23*k24 - k24*k24 + 2*k23*k34 + k24*k34 + 2*k23*(m1*m1) + k24*(m1*m1) + 2*k25*(m1*m1)) +
                       k12*(k14*k23 + 2*k14*k34 - k15*(m1*m1) - k23*(m1*m1) + 3*k25*(m1*m1) + k35*(m1*m1) + k13*(-4*k14 - 4*k23 - 3*k24 + 2*k34 + m1*m1)))/(4.*(k25*k25)) -
                      (k14*k14*k23 + 2*k14*(k23*k23) + k13*k13*(k23 + 2*k24) - k14*k23*k34 - k14*k23*(m1*m1) + k12*k12*(k13 + k23 + 2*k34 + 3*(m1*m1)) + k13*(k23*k23 + k24*k24 - k24*(k34 + m1*m1)) +
                       k12*(-(k13*k13) + 2*(k34*k34) + 3*pow(m1,4) + k13*(-k14 - 2*k23 - 2*k24 - k34 - m1*m1) - (k23 + k24)*(k23 + m1*m1) - k14*(2*k23 + k34 + m1*m1) + k34*(-k23 - k24 + 5*(m1*m1))))/(2.*k15*k25))*pow(pi,3))/
   (9.*pow(k34 + m1*m1,2));
    */
   double unsu2 = (16*Alfa2*Alfas*(-(6*(k14*k14)*k23 - 10*k14*(k23*k23) + 2*(k13*k13)*k24 + 6*k14*k23*k24 + 9*(k12*k12)*k34 - 22*k14*k23*k34 - 16*k14*k24*k34 +
                      k13*(-16*(k23*k23) + 2*(k24*k24) + 2*k14*(3*k23 + k24) + k23*(2*k24 - 16*k34) - 18*k24*k34) +
                      k12*(16*(k23*k23) + k13*(-2*k24 - 9*k34) + 7*k23*k34 + 7*k24*k34 + 9*(k34*k34) - k14*(6*k23 + 9*k34)))/(2.*(k15*k15)) -
                    (8*(k12*(k12 + k14) - k14*k24)*(-k13 - k23 + k34))/(k35*k35) +
                    (-16*(k14*k14)*k23 - 32*k14*(k23*k23) + 32*(k13*k13)*k24 - 16*k14*k23*k24 - 4*(k12*k12)*(5*k13 + k23 - 8*k34) +
                     k13*(32*k23*k24 + 32*(k24*k24) + k14*(-16*k23 + 16*k24) - 16*k24*k34) +
                     k12*(20*(k13*k13) - 12*k14*k23 + 4*(k23*k23) - 12*k23*k24 + k13*(20*k14 + 8*k23 - 12*k24 - 36*k34) + 12*k23*k34 - 16*k24*k34 + 16*(k34*k34)))/(2.*k15*k35))*
    pow(pi,3))/(9.*(k24-1e+4)*(k24-1e+4));
   return unsu2;
}

// eq. 5.2 in CS'2002; color structure is inside of Born
double CS_FF_Dipole::get_dipoles_value(const vector<Vec4D<double>>& p) {
   return - Born(ptilde_) * Vijk() / (2. * p_[i_]*p_[j_]);
}

double CS_FF_Dipole::ColorMatrix(const string& col_str1, const string& col_str2 ) {
   ColorFull::Col_amp Ca1(col_str1);
   ColorFull::Col_amp Ca2(col_str2);

   ColorFull::Col_functions Col_fun;
   // FormCalc numbers particles from 1, hence i_+1
   double temp = Col_fun.double_num(
      Col_fun.scalar_product(
         // emit gluon from emitter i_ and attach it to the spectator k_
         Col_fun.emit_gluon(Ca1, j_+1, 77), Col_fun.emit_gluon(Ca2, k_+1, 77)
      )
   );
   return temp/CF;
}

