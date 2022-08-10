#include "XSection_HnonC.hpp"
#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include "cuba.h"
#include "LHAPDF/LHAPDF.h"
#include "rk/rk.hh"
#include "rk/geom3.hh"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <execution>
#include <iostream>

// @todo this is absolutely necessary but I don't know why
// without it the uu -> uLuRg gives for BMP1 57.1.. +/- 0.005
// with it                                   57.8.. +/- 0.005

namespace {

// x[*nvec][*ndim] and f[*nvec][*ncomp]
int forwarder(const int *ndim, const double xx[][7],
   const int *ncomp, double ff[], void *userdata,
   int const* nvec, int const* cores, double const* weight, int const* iter) {
   std::vector<std::array<double, 7>> input;
   input.resize(*nvec);
   for(int i=0; i<*nvec; ++i) {
      for(int j = 0; j<7; ++j) {
         input.at(i).at(j) = xx[i][j];
      }
   }
   std::vector<double> output;
   output.resize(*nvec);
   std::transform(
      std::execution::par_unseq,
      input.cbegin(), input.cend(), output.begin(),
      [&userdata](std::array<double, 7> const& x) -> double {
         return static_cast<XSection_HnonC*>(userdata)->integrand(x);
      }
   );
   for (int i=0; i<*nvec; ++i) {
      ff[i] = output.at(i);
   }
   return 0;
}

} // anonymous

std::array<double, 3> XSection_HnonC::integrate() {

   //  integral dimension, number of integrands
   static constexpr int ndim = 7;
   static constexpr int ncomp = 1;
   //  accuraccy
   const double accuracy_rel {std::pow(10., -integration_precision_)};
   static constexpr double accuracy_abs {1e-12};

   static constexpr int neval_min = 10'000;
   long long int neval;
   static constexpr long long int neval_max {1'000'000'000'000};

   // technical (Vegas specific) stuff
   static constexpr int nstart = 200'000;
   static constexpr int nincrease = 10000;
   static constexpr int nbatch = 1000;
   static constexpr int gridno = 0;
   static constexpr int seed = 0;
   static constexpr int nvec = 1000;
   const char* state_file = "";
   int nregions, fail;

   double integral[ncomp], error[ncomp], prob[ncomp];

   const char* env_cubacores = std::getenv("CUBACORES");
   int nn = 0;
   const int pn = 10'000; // this is Cuba's default, see arXiv:1408.6373
   cubacores(&nn, &pn);

   llVegas(ndim, ncomp, (integrand_t)forwarder, this, nvec,
      accuracy_rel, accuracy_abs, integration_verbosity_, seed,
      neval_min, neval_max, nstart, nincrease, nbatch,
      gridno, state_file, NULL,
      &neval, &fail, integral, error, prob );

   /* @todo:
    * if CUBACORES was undefined, Cuba automatically set number of cores
    * based on the current system load. I don't know how to get this value
    * so I cannot restore it... */
   if (env_cubacores) {
      const int ncores = std::atoi(env_cubacores);
      cubacores(&ncores, &pn);
   }

   std::array <double, 3> result_finite
      {integral[0], error[0], prob[0]};

   return result_finite;
}

double XSection_HnonC::integrand(std::array<double, 7> const& xx) {

   const double m_sqr = Sqr(m1_);

   /*
   *  3-body phase space parametrization based on
   *  http://www.t39.ph.tum.de/T39_files/T39_people_files/duell_files/Dipl-MultiPion.pdf
   */

   // failsafe (this should never happen)
   // but sometimes does for suave
   assert(xx[0] >= 0 && xx[0] <= 1 // outgoing parton energy
       && xx[1] >= 0 && xx[1] <= 1 // particle 1 energy
       && xx[2] >= 0 && xx[2] <= 1 // angle
       && xx[3] >= 0 && xx[3] <= 1 // angle
       && xx[4] >= 0 && xx[4] <= 1 // angle
       && xx[5] >= 0 && xx[5] <= 1 // Bjorken x
       && xx[6] >= 0 && xx[6] <= 1 // Bjorken x
   );

   const double S = Sqr(sqrtS_);
   const double x1 = 4.*m_sqr/S      + (1.-4.*m_sqr/S)      * xx[5];
   const double x2 = 4.*m_sqr/(S*x1) + (1.-4.*m_sqr/(S*x1)) * xx[6];
   const double shat = x1*x2*S;
   const double shat_sqrt = std::sqrt(shat);

   const double Ej_max = 0.5*shat_sqrt - 2.*m_sqr/shat_sqrt;

   if (Ej_max < 0.5*dS_*shat_sqrt) {
      return 0;
   }

   const double Ej = 0.5*dS_*shat_sqrt + (Ej_max - 0.5*dS_*shat_sqrt)*xx[0];

   const double c = shat - 2.*shat_sqrt*Ej;
   const double sqrtc = std::sqrt(Sqr(c) - 4.*c*m_sqr);
   // Eq. 4.5
   const double E1_max = ((shat_sqrt - Ej)*c + Ej*sqrtc)/(2.*c);
   const double E1_min = ((shat_sqrt - Ej)*c - Ej*sqrtc)/(2.*c);
   const double E1 = E1_min + (E1_max - E1_min)*xx[1];
   assert(E1 >= m1_);

   // Eq. 4.2 with E2 = Ej
   const double p1 = std::sqrt((E1-m1_)*(E1+m1_));
   const double cosx = (shat - 2*shat_sqrt*(E1+Ej) + 2.*Ej*E1)/(2.*Ej*p1);

   // check if due to numerics |cos(x)| is not > 1
   // if yes, return 0 and continue
   if ( cosx > 1 || cosx < -1)  {
      std::cout << "Warning! 1 - |cos(x)| = " << 1 - std::abs(cosx) << "  - Skipping the phase space point.\n";
      return 0.;
   }

   // initialize matrix of particles momenta
   std::array<std::array<double, 4>, 5> p;

   // incoming partons momenta
   p[0] = {0.5*shat_sqrt, 0., 0.,  0.5*shat_sqrt};
   p[1] = {0.5*shat_sqrt, 0., 0., -0.5*shat_sqrt};

   // 1st final state momenta
   const double sinpix2 = std::sin(pi*xx[2]);
   const double cospix2 = std::cos(pi*xx[2]);
   const double sin2pix3 = std::sin(two_pi*xx[3]);
   const double cos2pix3 = std::cos(two_pi*xx[3]);
   p[2][0] = E1;
   p[2][1] = p1 * sinpix2 * cos2pix3;
   p[2][2] = p1 * sinpix2 * sin2pix3;
   p[2][3] = p1 * cospix2;

   const geom3::UnitVector3 rotation_axis(
      sinpix2 * cos2pix3,
      sinpix2 * sin2pix3,
      cospix2
   );

   // construct rotation by angle xx[4] * two_pi
   // around final state momentum
   const geom3::Rotation3 rot(rotation_axis, xx[4]*two_pi);

   /*
    *  kinematics was solved for the cosx of angle
    *  between parton and SUSY particle
    *
    *  probably it would be also ok to skip it hoping that
    *  periodicity of sin and cos would solve the thing
    */
   double parton_theta, parton_phi;
   if( pi*xx[2] + std::acos(cosx) < pi ) {
      parton_theta = pi*xx[2] + std::acos(cosx);
      parton_phi = two_pi*xx[3];
   }
   else {
      parton_theta = two_pi - pi*xx[2] - std::acos(cosx);
      parton_phi = two_pi*xx[3] + pi;
   }

   const double sinparton_theta = std::sin(parton_theta);
   const geom3::Vector3 p_parton(
      Ej * sinparton_theta * std::cos(parton_phi),
      Ej * sinparton_theta * std::sin(parton_phi),
      Ej * std::cos(parton_theta)
   );

   // rotate parton momentum
   const geom3::Vector3 p_temp = rot.rotate(p_parton);

   // set parton momenta
   const std::array<double, 4> p_temp_2 = {Ej, p_temp.x(), p_temp.y(), p_temp.z()};

   // 2nd sgluon momenta
   p[3][0] = shat_sqrt - E1 - Ej;
   for (int i = 1; i < 4; ++i) p[3][i] = - p[2][i] - p_temp_2[i];

   assert( std::abs(
      (pow(p[2][0], 2) - pow(p[2][1], 2) - pow(p[2][2], 2) - pow(p[2][3], 2))/(m1_ * m1_) - 1) < 1e-10
         && p[2][0] >= m1_
   );

   // write parton momentum to momentum matrix p
   for(int i = 0; i < 4; ++i) p[4][i] = p_temp_2[i];

   const double t15 = -2.*(p[0][0]*p[4][0] - p[0][3]*p[4][3]);
   const double t25 = -2.*(p[1][0]*p[4][0] - p[1][3]*p[4][3]);

   // check if we are not in the collinear region
   // if yes, return
   if (-t15 < dC_*shat_sqrt*Ej || -t25 < dC_*shat_sqrt*Ej) {
      return 0.;
   }

   double ME2 = f(pdf_->alphasQ(muR_), p);
   assert(!std::isnan(ME2) && ME2 >= 0);
   /*
   std::cout << setprecision(17);
   for (const auto& row : p) {
      for (const double& col : row) {
         std::cout << col << ' ';
      }
      std::cout << '\n';
   }
   std::cout << ME2 << '\n';
   std::cout << "================================================\n";
   */

   // some final factors
   ME2 *= to_fb/shat*sinpix2/(128.*pi_sqr);

   double pdf_flux = 0.0;
   for (const auto& f : flav_) {
      pdf_flux += f.at(2) * pdf_->xfxQ(f.at(0), x1, muF_) * pdf_->xfxQ(f.at(1), x2, muF_);
   }
   ME2 *=  pdf_flux/(x1*x2);

   const double jacobian =  (xx[5]*(-4*dS_*Sqr(m1_) + (-1 + dS_)*xx[5]*xx[6]*(-S + 4*Sqr(m1_)))*
     (xx[5]*xx[6]*xx[0]*(S - 4*Sqr(m1_)) +
       dS_*(-1 + xx[0])*(-(S*xx[5]*xx[6]) + 4*(-1 + xx[5]*xx[6])*Sqr(m1_)))/S*
     Sqr(S - 4*Sqr(m1_))*std::pow(S*xx[5] - 4*(-1 + xx[5])*Sqr(m1_),-1)*
     std::pow(S*xx[5]*xx[6] + (4 - 4*xx[5]*xx[6])*Sqr(m1_),-1)*
     std::pow((-1 + dS_)*S*xx[5]*xx[6]*(-1 + xx[0]) +
       4*(1 + xx[5]*xx[6]*(-1 + xx[0]) - dS_*(-1 + xx[5]*xx[6])*(-1 + xx[0]))*Sqr(m1_),-1)*
     std::sqrt((-1 + xx[0])*(S*xx[5]*xx[6]*(-1 + dS_ + xx[0] - dS_*xx[0]) +
         4*(-1 + xx[5]*xx[6] + dS_*(-1 + xx[5]*xx[6])*(-1 + xx[0]) - xx[5]*xx[6]*xx[0])*Sqr(m1_))*
       (-4*dS_*Sqr(m1_) + (-1 + dS_)*xx[5]*xx[6]*(-S + 4*Sqr(m1_)))))*0.25;

   return ME2 * std::abs(jacobian);
}
