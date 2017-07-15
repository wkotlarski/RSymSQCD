#include <iostream>
#include "CSDipoles.hpp"
#include "constants.hpp"

//using namespace ColorFull;
using namespace std;

CSDipole::CSDipole(double (*me)(valarray<valarray<double>> const&), 
                   int j, int i, int k,  
                   valarray<valarray<double>> const& p) {
   assert( i != j && i != k && j !=k);
   i_ = i;
   j_ = j;
   k_ = k;
   p_ = p;
   q_ = p_[i_] + p_[j_] + p_[k_];
   ptilde_ = ff_mom_shuffle();
   born_me2_ = me;
}

// eq. 5.7 of CS'1997
inline double CSDipole::lambda(double x, double y, double z) {
   return x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
}

// eq. 5.6: 
inline double CSDipole::beta(valarray<double> p, valarray<double> q) {
   double res { sqrt(lambda(dot(p+q,p+q), dot(p,p), dot(q,q)))/(dot(p+q,p+q) - dot(p,p) - dot(q,q))};
   constexpr double prec = 1e-3;
   if(res >= 0. && res <= 1.) {
      return res;
   } else if ( abs(res - 1) < prec ) {
      return 1;
   } else if ( abs(res) < prec ) {
      return 0.;
   } else {
      cout << "Velocity " << 1-res << endl;
   }
}

inline double CSDipole::dot(const valarray<double>& v1, const valarray<double>& v2) {
   assert (v1.size() == 4 && v2.size() == 4);
   double sum = v1[0] * v2[0];
   for(int i = 1; i < 4; ++i) sum -= v1[i] * v2[i];
   return sum;
}

inline double CSDipole::yijk() {
   return dot(p_[i_], p_[j_])/(dot(p_[i_], p_[j_]) + dot(p_[i_], p_[k_]) + dot(p_[j_], p_[k_]));
}

inline double CSDipole::z( int l ) {
   return dot(p_[l], p_[k_])/(dot(p_[i_], p_[k_]) + dot(p_[j_], p_[k_]));
}

inline double CSDipole::Vijk() {
   return 8 * pi * Alfas * CF * (
      2/(1 - z(j_)*(1-yijk())) 
      - beta(ptilde_[i_], ptilde_[k_])/beta(p_[i_]+p_[j_], p_[k_]) * (1 + z(j_) + mi*mi/dot(p_[i_], p_[j_]))
   );
}

inline double CSDipole::x(int i) {
   return 2 * dot(p_[i], q_)/dot(q_, q_);
}

inline double CSDipole::Born(valarray<valarray<double>> const& p) {
   double k12 = dot(ptilde_[0], ptilde_[1]);
   double k23 = dot(ptilde_[1], ptilde_[2]);
   double k24 = dot(ptilde_[1], ptilde_[3]);
   double k13 = dot(ptilde_[0], ptilde_[2]);
   double k14 = dot(ptilde_[0], ptilde_[3]);
   double k34 = dot(ptilde_[2], ptilde_[3]);
   const string cos = "[{3, 4}]";
   return Mat(cos,cos) * 8 * Alfa2 * pi_sqr * (mi*mi*k12 + k14*k23 + k13*k24)/(9.*pow(k12,2));
}

// eq. 5.9: combine particles i & j -> ij using k as the spectator
valarray<valarray<double>> CSDipole::ff_mom_shuffle() {
   const double q2 = dot(q_, q_);
   auto mij {mi};
   valarray<valarray<double>> ptilde__ (4);
   ptilde__[0] = p_[0];
   ptilde__[1] = p_[1];
   ptilde__[k_] = 
      0.5*(q2 + mk*mk - mij*mij)/q2 * q_
      + sqrt(lambda(q2, mij*mij, mk*mk)/lambda(q2, dot(p_[i_] + p_[j_], p_[i_] + p_[j_]), mk*mk)) 
         * (p_[k_] - dot(q_, p_[k_])/q2 * q_);
   ptilde__[i_] = q_ - ptilde__[k_];
   return ptilde__;
}

double CSDipole::unsu(const valarray<valarray<double>>& p) {
   
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
   // double x1 = x(p,2);
   // double x2 = x(p,3);
   //double unsu = 1/(1-x2) * (2/(2-x1-x2) - (1+x1)) 
   //   + 1/(1-x1) * (2/(2-x1-x2) - (1+x2));  
   double MB2 = pow(0.,2);
   double unsu2 = (-512*Alfa2*Alfas*((-2*k14*k24*MB2 - k12*k35*MB2 + k12*MB2*(k12 + MB2) + (-2*k14*k24 + k12*(k14 + k24))*(-k13 - k23 + k34 + MB2))/(4.*(k35*k35)) - 
       (k15*k15*k23 - k15*(k23*k23) - 2*(k12*k12)*(k13 + k23) - k13*k13*k25 - k15*k23*k25 - 2*k15*k23*k35 - 2*k14*k23*MB2 + k13*(k15*(k23 - k25) + k25*(k23 + k25) + 2*(-2*k23 - k25)*k35 - 2*k24*MB2) + 
          k12*(4*(k13*k13) + 4*(k23*k23) + k23*(2*k15 + 3*k25 - 2*k35) - (k15 + k25 - 2*k35)*k35 + k13*(3*k15 - 2*(-2*k23 - k25 + k35)) + (-k15 - k25 + 2*(k34 + k35))*MB2))/(4.*k35*k45) + 
       (k15*(-(k13*k23) + k23*k23) + k13*k13*k25 - k13*k23*k25 + 2*k13*k23*k35 + k12*k12*MB2 - 2*k13*k23*MB2 + k12*(MB2*MB2 + k15*(k13 + MB2) + k25*(k23 + MB2) - k35*(k13 + k23 + MB2)))/(4.*(k45*k45)))*pow(pi,3))/
   (9.*(k12*k12));
   assert(unsu2 > 0);
   return unsu2;
}
double CSDipole::eval_dipole(const valarray<valarray<double>>& p) {
   double s12 = 2 * dot(p_[0], p_[1]);
   /*
   double x1, x2;
   if (i_ == 2) {
      x1 = x(p, i_); x2 = x(p, k_); 
   }
   else {
      x1 = x(p, k_); x2 = x(p, i_); 
   }
   */
   double mu = mi/sqrt(s12);
   double mu2 = pow(mu, 2);
   //return 1/(1-x2) * (2/(2-x1-x2)-(1+x1))+(1-x1)/x2 
   //   + 1/(1-x1) * (2/(2-x1-x2) - (1+x2))+(1-x2)/x1;
   //cout << i_ << ' ' << j_ << ' ' << k_ << ' ' << dot(temp[i_], temp[j_]) << ' ' << dot(p_in[i_], p_in[j_]) << '\n';
   return Born(ptilde_) * Vijk() / (2 * dot(p_[i_], p_[j_]));
}

double CSDipole::Mat(const string& col_str1, const string& col_str2 ) {
   //ColorFull::Col_amp Ca1(col_str1);
   //ColorFull::Col_amp Ca2(col_str2);
   //ColorFull::Col_functions Col_fun;
   // FormCalc numbers particles from 1, hence i_+1
   /*
   double temp = Col_fun.double_num(
      Col_fun.scalar_product(
         // emit gluon from emitter i_ and attach it to the spectator k_
         Col_fun.emit_gluon(Ca1, i_+1, 77), Col_fun.emit_gluon(Ca2, k_+1, 77)
      )
   );
   */
   return 3.; //temp/CF;
}

