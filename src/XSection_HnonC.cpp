#include <cmath>
#include <cassert>

#include "XSection_HnonC.hpp"
#include "constants.hpp"
// needed to do Euler rotation
#include "rk/rk.hh"
#include "rk/geom3.hh"

std::array<double, 3> XSection_HnonC::integrate() {

   //  integral dimension, number of integrands
   static constexpr int ndim = 7;
   static constexpr int ncomp = 1;
   //  accuraccy
   const double accuracy_rel {std::pow(10., -vm["precision-hard"].as<int>())};
   static constexpr double accuracy_abs {1e-12};

   static constexpr int neval_min = 10'000;
   long long int neval;
   static constexpr long long int neval_max {1'000'000'000'000};

   // technical (Vegas specific) stuff
   static constexpr int nstart = 200'000;
   static constexpr int nincrease = 10000;
   static constexpr int nbatch = 1000;
   static constexpr int gridno = 0;
   const int flags = vm["verbosity-hard"].as<int>();
   static constexpr int seed = 0;
   const char* state_file = "";
   int nregions, fail;

   cubareal integral[ncomp], error[ncomp], prob[ncomp];
   llVegas( ndim, ncomp, integrand, NULL, 1,
      accuracy_rel, accuracy_abs, flags, seed,
      neval_min, neval_max, nstart, nincrease, nbatch,
      gridno, state_file, NULL,
      &neval, &fail, integral, error, prob );

   std::array <double, 3> result_finite
      {integral[0], error[0], prob[0]};

   return result_finite;
}

int XSection_HnonC::integrand(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

   /*
   *  3-body phase space parametrization based on
   *  http://www.t39.ph.tum.de/T39_files/T39_people_files/duell_files/Dipl-MultiPion.pdf
   */

   const double m_sqr = m1*m1;
   double x1 = 1;
   double x2 = 1;
   x1 = 4.*m_sqr/S      + (1. - 4.*m_sqr/S)      * xx[*ndim-2];
   x2 = 4.*m_sqr/(S*x1) + (1. - 4.*m_sqr/(S*x1)) * xx[*ndim-1];

   const double shat = x1*x2*S;
   const double shat_sqrt = std::sqrt(shat);

   /*  s35   = (p3 + p5)^2 = (p1 + p2 - p4)^2
    *        = S - 2 (p1 + p2) * p4 + m2^2
    *  since p1 + p2 doesn't have 3-momentum component, the result depends only on energy
    *        = S - 2 sqrtS * E2 + m2^2
    */
   const double s35_min = m1*m1;
   const double s35_max = Sqr(shat_sqrt - m2);
   const double s35 = s35_min + (s35_max - s35_min)*xx[0];
   const double E2 = -0.5*(s35 - shat - m1*m1)/shat_sqrt;
   assert(E2 >= m2);
   const double p2 = std::sqrt((E2-m2)*(E2+m2));

   const double c = shat - 2.*shat_sqrt*E2 + m1*m1 + m2*m2;
   // Eq. 4.5
   const double E1_max = ((shat_sqrt-E2)*c + p2 * std::abs(c-2.*m_sqr))/(2.*(c-m_sqr));
   const double E1_min = ((shat_sqrt-E2)*c - p2 * std::abs(c-2.*m_sqr))/(2.*(c-m_sqr));

   // analogusly to s35, s45 = S - 2 sqrtS * E1 + m1^2
   const double s45_min = shat + m2*m2 - 2 * shat_sqrt * E1_max;
   const double s45_max = shat + m2*m2 - 2 * shat_sqrt * E1_min;
   const double s45 = s45_min + (s45_max - s45_min)*xx[1];
   const double E1 = -0.5*(s45 - shat - m2*m2)/shat_sqrt;
   assert(E1 >= m1);

   ff[0] = 0.;
   if (shat_sqrt - E1 - E2 < 0.5*dS*shat_sqrt) {
      return 0;
   }

   const double p1 = std::sqrt((E1-m1)*(E1+m1));

   // eq. 4.2
   const double cosx = 0.5*(shat - 2.*shat_sqrt*(E1+E2) + 2.*E2*E1 + m1*m1 + m2*m2)/(p1*p2);

   // check if due to numerics |cos(x)| is not > 1
   // if yes but reasonable, return 0 and continue
   // assert(std::abs(cosx) - 1. < 1e-7);
   if (cosx > 1. || cosx < -1.)  {
      // std::cout << "Warning! 1 - |cos(x)| = " << 1. - std::abs(cosx) << "  - Skipping the phase space point.\n";
      return 0;
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

   geom3::UnitVector3 rotation_axis(
      sinpix2 * cos2pix3,
      sinpix2 * sin2pix3,
      cospix2
   );

   // construct rotation by angle xx[4] * two_pi
   // around final state momentum
   geom3::Rotation3 rot(rotation_axis, xx[4]*two_pi);

   /*
    *  kinematics was solved for the cosx of angle
    *  between parton and SUSY particle
    *
    *  probably it would be also ok to skip it hoping that
    *  periodicity of sin and cos would solve the thing
    */
   double parton_theta, parton_phi;
   if (pi*xx[2] + std::acos(cosx) < pi) {
      parton_theta = pi*xx[2] + std::acos(cosx);
      parton_phi = two_pi*xx[3];
   }
   else {
      parton_theta = two_pi - pi*xx[2] - std::acos(cosx);
      parton_phi = two_pi*xx[3] + pi;
   }

   const double sinparton_theta = std::sin(parton_theta);
   const geom3::Vector3 p_parton(
      p2 * sinparton_theta * std::cos(parton_phi),
      p2 * sinparton_theta * std::sin(parton_phi),
      p2 * std::cos(parton_theta)
   );

   // rotate parton momentum
   const geom3::Vector3 p_temp = rot.rotate(p_parton);

   // 2nd final state momenta
   p[3] = {E2, p_temp.x(), p_temp.y(), p_temp.z()};

   // parton
   p[4][0] = shat_sqrt - E1 - E2;
   assert (p[4][0] >= 0);
   for(int i = 1; i < 4; ++i) {
      p[4][i] = - p[2][i] - p[3][i];
   }

   const double t15 = -2.*(p[0][0]*p[4][0] - p[0][3]*p[4][3]);
   const double t25 = -2.*(p[1][0]*p[4][0] - p[1][3]*p[4][3]);

   // check if we are not in the collinear region
   // if yes, return
   if (-t15 < dC*shat_sqrt*p[4][0] || -t25 < dC*shat_sqrt*p[4][0]) {
      return 0;
   }

   double ME2 = (processID->*processID->matrixelementReal_HnonC)(p);
   assert(!std::isnan(ME2) && ME2 >= 0);

   // some final factors
   const double fac = to_fb/(128.*shat*pi_sqr) * sinpix2;

   double pdf_flux = 0.0;
   for (const auto& f : processID->flav) {
      pdf_flux += f.at(2) * pdf->xfxQ(f.at(0), x1, mu_f) * pdf->xfxQ(f.at(1), x2, mu_f);
   }
   pdf_flux /= x1 * x2;

   ME2 *=  pdf_flux;

   const double m {m1};
   const double J = ((-Sqr(m) + s35)*(-2*m + shat_sqrt)*shat_sqrt*std::sqrt(std::pow(m,4) + Sqr(s35 - shat) - 2*Sqr(m)*(s35 + shat)))/s35;
   ME2 *= std::abs(J);

   // Jakobian of (E2, E1) -> (s35, s45) change
   ME2 *= 0.25/shat;

   // Jakobian of Bjorken vars. mapping xx0, xx1 -> x1, x2
   const double xx0 = xx[*ndim-2];
   ME2 *= Sqr(-4.*std::pow(m,2) + S)*xx0/(S*(-4.*Sqr(m)*(-1. + xx0) + S*xx0));

   ME2 *= fac;

   ff[0] = ME2;

   return 0;
}
