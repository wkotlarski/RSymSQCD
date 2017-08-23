#include "CSDipole.hpp"
#include "../include/constants.hpp"
#include "ColorFull/Core/Tree_level_gluon_basis.h"
#include "dilog.hpp"

using namespace dilogarithm;

std::vector<Vec4D<double>> d_to_vec_2(double* d1) {
   return std::vector<Vec4D<double>> {
           Vec4D<double> {d1[1], d1[2], d1[3], d1[4]},
           Vec4D<double> {d1[5], d1[6], d1[7], d1[8]},
           Vec4D<double> {d1[9], d1[10], d1[11], d1[12]},
           Vec4D<double> {d1[13], d1[14], d1[15], d1[16]}
   };
}
extern "C" double* c_unintegrated_dipoles(
        int,
        double,
        double, double, double, double,
        double, double, double, double,
        double, double, double, double,
        double, double, double, double,
        double, double, double, double,
        int, int, int
);

extern "C" double c_integrated_dipoles(
        int,
        double,
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

   double* d1 = c_unintegrated_dipoles(
           type_to_int,
           Born_.pdf->alphasQ(Born_.mu_r),
           p[0].t_, p[0].x_, p[0].y_, p[0].z_,
           p[1].t_, p[1].x_, p[1].y_, p[1].z_,
           p[2].t_, p[2].x_, p[2].y_, p[2].z_,
           p[3].t_, p[3].x_, p[3].y_, p[3].z_,
           p[4].t_, p[4].x_, p[4].y_, p[4].z_,
           4, emit_, spec_
   );
   auto res = Born_.get_ME2_value(emit_, spec_, d_to_vec_2(d1))*d1[0]/CF;
   free_mom(d1);
   return res;
}

double CSDipole::eval_integrated_dipole(std::vector<Vec4D<double>> const& p) const {
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

   double d1 = c_integrated_dipoles(
           type_to_int,
           Born_.pdf->alphasQ(Born_.mu_r),
           p[0].t_, p[0].x_, p[0].y_, p[0].z_,
           p[1].t_, p[1].x_, p[1].y_, p[1].z_,
           p[2].t_, p[2].x_, p[2].y_, p[2].z_,
           p[3].t_, p[3].x_, p[3].y_, p[3].z_,
           4, emit_, spec_
   );
   auto res = Born_.get_ME2_value(emit_, spec_, p)*d1;
   std::cout << "I: " << res << std::endl;
   return res;
}

double CSDipole::eval_P(std::vector<Vec4D<double>> const& p, double x) const {
   // P-term appears only for the initial state emitters
   if (type_ == DipoleType::FF || type_ == DipoleType::FI) return 0.;
   double sja = 2.*p[emit_]*p[spec_];
   double mu_f = Born_.mu_f;
   double alpha_s = Born_.pdf->alphasQ(Born_.mu_r);
   // eq. 6.54 of CS'02
   return alpha_s/(2.*pi) * 1./CF * Born_.get_ME2_value(spec_, emit_, p)*log(mu_f/sja);
}

double CSDipole::eval_K(std::vector<Vec4D<double>> const& p, double x) const {
   // K-term appears only for the initial state emitters
   if (type_ == DipoleType::FF || type_ == DipoleType::FI) return 0.;
   double sja = 2.*p[emit_]*p[spec_];
   double alpha_s = Born_.pdf->alphasQ(Born_.mu_r);
   // eq. 6.55 of CS'02
   return alpha_s/(2.*pi);
}
