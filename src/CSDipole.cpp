#include "CSDipole.hpp"
#include "../include/constants.hpp"
#include "ColorFull/Core/Tree_level_gluon_basis.h"

std::vector<Vec4D<double>> d_to_vec_2(double* d1) {
   return std::vector<Vec4D<double>> {
           Vec4D<double> {d1[1], d1[2], d1[3], d1[4]},
           Vec4D<double> {d1[5], d1[6], d1[7], d1[8]},
           Vec4D<double> {d1[9], d1[10], d1[11], d1[12]},
           Vec4D<double> {d1[13], d1[14], d1[15], d1[16]}
   };
}
extern "C" double* c_interface_to_cs_unintegrated_dipoles(
        int,
        double, double, double, double,
        double, double, double, double,
        double, double, double, double,
        double, double, double, double,
        double, double, double, double,
        int, int, int
);

extern "C" void free_mom(double*);

double CSDipole::eval_unintegrated_dipole(std::vector<Vec4D<double>> const& p) const {
   int type_to_int;
   if (type_ == DipoleType::FF) {
      type_to_int = 1;
   } else if (type_ == DipoleType::II) {
      type_to_int = 2;
   } else if (type_ == DipoleType::IF) {
      type_to_int = 3;
   } else if (type_ == DipoleType::FI) {
      type_to_int = 4;
   }

   double* d1 = c_interface_to_cs_unintegrated_dipoles(
           type_to_int,
           p[0].t_, p[0].x_, p[0].y_, p[0].z_,
           p[1].t_, p[1].x_, p[1].y_, p[1].z_,
           p[2].t_, p[2].x_, p[2].y_, p[2].z_,
           p[3].t_, p[3].x_, p[3].y_, p[3].z_,
           p[4].t_, p[4].x_, p[4].y_, p[4].z_,
           4, i_, j_
   );
   auto res = Born_.get_ME2_value(d_to_vec_2(d1))*d1[0];
   free_mom(d1);
   return res;
}
double lambda(double x, double y, double z) {
   x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
}
double vtijk(double s12) {
   double mu2 = 1500*1500/s12;
   return sqrt(lambda(1., mu2, mu2))/(1. - 2*mu2);
}

double rho(int n, double mu_j, double mu_k) {
   return sqrt(1 - )
}
double CSDipole::eval_integrated_dipole(std::vector<Vec4D<double>> const& p) const {
   return 0.;//Born_.get_ME2_value(p);
}
