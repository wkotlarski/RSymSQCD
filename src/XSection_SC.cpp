#include <numeric>
#include "XSection_SC.hpp"

std::vector<CSDipole> XSection_SC::cs_dipoles;

std::vector<Vec4D<double>> s12_and_costh_to_p (double s, double costh) {
   double p = sqrt(s/4. - pow(1500,2));
   std::vector<Vec4D<double>> temp = {
           Vec4D<double> { sqrt(s)/2., 0, 0, sqrt(s)/2.},
           Vec4D<double> { sqrt(s)/2., 0, 0, -sqrt(s)/2.},
           Vec4D<double> { sqrt(s)/2., 0., p*sqrt(1-costh*costh), p*costh},
           Vec4D<double> { sqrt(s)/2., 0., -p*sqrt(1-costh*costh), -p*costh}
   };
   return temp;
}
std::array<double, 3> XSection_SC::integrate() {
  
   // integrals dimensions, number of integrands
   constexpr int ndim { 4 }, ncomp { 1 };
   // accuraccy
   double accuracy_rel_sc { pow( 10., -vm["precision-sc"].as<int>() ) }, 
            accuracy_rel_c { pow( 10., -vm["precision-sc"].as<int>() ) };
   constexpr double accuracy_abs { 1e-12 };

   constexpr int neval_min = 10'000;
   long long int neval;
   constexpr long long int neval_max { 1'000'000'000 }; 
   const int verbose = vm["verbosity-sc"].as<int>();

   int nregions, fail;

   cubareal integral_c1[ncomp], error_c1[ncomp], prob_c1[ncomp];
   llCuhre(ndim, ncomp, integrand_c1, NULL, 1,
      accuracy_rel_c, accuracy_abs, verbose,
      neval_min, neval_max, 1, NULL, NULL,
      &nregions, &neval, &fail, integral_c1, error_c1, prob_c1);

   cubareal integral_c2[ncomp], error_c2[ncomp], prob_c2[ncomp];
   llCuhre(ndim, ncomp, integrand_c2, NULL, 1,
      accuracy_rel_c, accuracy_abs, verbose,
      neval_min, neval_max, 1, NULL, NULL,
      &nregions, &neval, &fail, integral_c2, error_c2, prob_c2);

   std::array <double, 3> result_finite { 
      integral_c1[0] + integral_c2[0],
      sqrt( pow ( error_c1[0], 2) + pow ( error_c2[0], 2) ),
      + prob_c2[0]
   };
  
  return result_finite;
}

int XSection_SC::integrand_sc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
    
   // integration variables
   double x1 = 4. * pow(m1, 2)/S + (1 - 4. * pow(m1, 2)/S ) * xx[0];
	double x2 = 4. * pow(m1, 2) /(S * x1) + (1 - 4. * pow(m1, 2)/(S * x1)) * xx[1];
   double th = xx[2] * pi;
    
   double s12 = x1 * x2 * S;
   double Alfas = pdf->alphasQ( mu_r );
    
   double pdf_flux = 0.0;
//   for (const auto& f : processID->flav) {
//      pdf_flux += f.at(2) * pdf->xfxQ( f.at(0), x1, mu_f ) * pdf->xfxQ( f.at(1), x2, mu_f );
//   }
   pdf_flux /= x1 * x2;
   
//   ff[0] = pdf_flux * (processID->*processID->matrixelementReal_SC)(s12, th);
   ff[0] *= to_fb;
    
   // jakobian
   ff[0] *= pi*pow(-4*pow(m1, 2) + S,2)*xx[0] / 
           (S*(-4*pow(m1, 2)*(-1 + xx[0]) + S*xx[0]));
   return 0;
}

// z = 1 piece of K & P terms
int XSection_SC::integrand_c1(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {
    
   // integration variables
   double x1 = 4. * pow(m1, 2)/S + (1 - 4. * pow(m1, 2)/S ) * xx[0];
   double x2 = 4. * pow(m1, 2) /(S * x1) + (1 - 4. * pow(m1, 2)/(S * x1)) * xx[1];
   double costh = -1. + 2.*xx[2];

   double s12 = x1 * x2 * S;
   double sja;
   auto p = s12_and_costh_to_p(s12, costh);
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow( Alfas, 2);
    
   double pdf_flux = 0.0;
//   for (const auto& inner : processID->flav) {
//      pdf_flux += inner.at(2) * pdf->xfxQ( inner.at(0), x1, mu_f ) * pdf->xfxQ( inner.at(1), x2, mu_f );
//   }
   pdf_flux /= x1 * x2;
   
   ff[0] = 0.;

   double GAMMA_Q = 3/2. * CF;
   double mj2 = pow(1500, 2);
   double mj = 1500;
   double Qja2;

   // x=1 part of K for IF
   for (const auto& dipole : cs_dipoles) {
      // eq. 6.55 of CS'02
      if (dipole.get_emitter() < 2 && dipole.get_spectator() > 1) {
         sja = 2.*p[dipole.get_emitter()]*p[dipole.get_spectator()];
         Qja2 = sja + mj2;
         double sig = 1.; //dipole.Born_.get_ME2_value(-1, -1, s12_and_costh_to_p(s12, costh));
         double sigij = 1; //dipole.Born_.get_ME2_value(dipole.get_spectator() , dipole.get_emitter() , s12_and_costh_to_p(s12, costh));
         ff[0] -= Alfas/(2.*pi) * (
            // z=1 part of K
            - (3/2.*CF + (7./2. - pi*pi/6.)*CF - 5/6.*pi*pi*CF) * sig
            // z=1 part of italic K
            // including z=1 part of equation C.15
            - (-3/2. + mj2/sja*log(mj2/(sja+mj2)) + 1/2.*mj2/(sja+mj2) - (mj2/sja*log(mj2/sja)+mj2/(2*Qja2)) + (3/2. - 2))*sigij
            // non K-terms (including remnants of P)
            - 1./CF * GAMMA_Q
            * (log((sja - 2*mj*sqrt(sja + mj2) + 2*mj2)/sja) + 2*mj/(sqrt(sja+mj2)+mj))*sigij
         );
      }
   }

   // x=1 part of K for II
   for (const auto& dipole : cs_dipoles) {
      // eq. 6.68 of CS'02
      if (dipole.get_emitter() < 2 && dipole.get_spectator() < 2) {
         //ff[0] += Alfas/(2.*pi) * (pi*pi/3.)
         //         * dipole.Born_.get_ME2_value(dipole.get_spectator() , dipole.get_emitter() , s12_and_costh_to_p(s12, costh));
      }
   }

   double beta = sqrt(1. - 4*1500*1500/s12);
   ff[0] *= pdf_flux * to_fb/(2.*s12) * beta/(16*pi);
   ff[0] *= 2.*(pow(-4*pow(m1, 2) + S,2)*xx[0])/(S*(-4*pow(m1, 2)*(-1 + xx[0]) + S*xx[0]));

    
   return 0;
}  

int XSection_SC::integrand_c2(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

   // scale integration variables as Cuba works in a unit hipercube
   double x1 = 4. * pow(m1, 2)/S + (1 - 4. * pow(m1, 2)/S ) * xx[0];
	double x2 = 4. * pow(m1, 2) /(S * x1) + (1 - 4. * pow(m1, 2)/(S * x1)) * xx[1];
   double z = x1 + ( 1 - dS - x1 ) * xx[2];
   double costh = -1. + 2.*xx[3];
   double s12 = x1 * x2 * S;

   if( x1/z > 1. || z*s12 < 4*m1*m1) {
      ff[0] = 0;
      return 0;
   }
     
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow( Alfas, 2);
   double mj2 = pow(1500, 2);
   auto p = s12_and_costh_to_p(s12, costh);
   double sja;
   double muQ = 1500/sqrt(s12);
   ff[0] = 0.0;

   // ------------------------------- P - terms ------------------------------------------------------------------------
   double dipole_sum1 = std::accumulate(
           cs_dipoles.begin(), cs_dipoles.end(), 0.,
           [s12,costh,z](double current,  CSDipole& el) {
              return current + el.eval_P(s12_and_costh_to_p(s12, costh), z);
           }
   );
   double dipole_sum2 = std::accumulate(
           cs_dipoles.begin(), cs_dipoles.end(), 0.,
           [s12,costh,z](double current,  CSDipole& el) {
              return current + el.eval_P(s12_and_costh_to_p(s12, costh), 1.);
           }
   );
//   for (const auto& f : processID->flav) {
//      for (const auto& dipole : cs_dipoles) {
//         if (dipole.get_emitter() == 0) {
//            ff[0] += f.at(2) * pdf->xfxQ(f.at(0), x1 / z, mu_f) / (x1 / z) * pdf->xfxQ(f.at(1), x2, mu_f) / x2
//                     * (-CF * (1. + z) * dipole_sum1/z + 2. * CF / (1. - z) * (dipole_sum1 - dipole_sum2));
//         }
//         else if (dipole.get_emitter() == 1) {
//            ff[0] += f.at(2) * pdf->xfxQ(f.at(0), x2, mu_f) / x2 * pdf->xfxQ(f.at(1), x1 / z, mu_f) / (x1 / z)
//                     * (-CF * (1. + z) * dipole_sum1/z + 2. * CF / (1. - z) * (dipole_sum1 - dipole_sum2));
//         }
//      }
//   }
   // ------------------------------------------------------------------------------------------------------------------

   // ------------------------------ z-dependent part of K for II ------------------------------------------------------
   for (const auto& dipole : cs_dipoles) {
      // eq. 6.68 of CS'02
      if (dipole.get_emitter() < 2 && dipole.get_spectator() < 2) {
         ff[0] += Alfas/(2.*pi) * (
         //dipole.Born_.get_ME2_value(dipole.get_spectator() , dipole.get_emitter() , s12_and_costh_to_p(z*s12, costh))
                   1 * (-(1.+z)*log(1-z) + 2.*log(1-z)/(1-z)) -
            // dipole.Born_.get_ME2_value(dipole.get_spectator() , dipole.get_emitter() , s12_and_costh_to_p(s12, costh))
           1 * 2.*log(1-z)/(1-z)
         );
      }
   }
   // ------------------------------------------------------------------------------------------------------------------

   // z-dependent part of K insertion operator for IF
//   for (const auto& f : processID->flav) {
//      for (const auto& dipole : cs_dipoles) {
//         // eq. 6.68 of CS'02
//         if (dipole.get_emitter() < 2 && dipole.get_spectator() > 1) {
//            sja = 2.*p[dipole.get_emitter()]*p[dipole.get_spectator()];
//            double diff = 1; //dipole.Born_.get_ME2_value(-1, -1 , s12_and_costh_to_p(z*s12, costh)) - dipole.Born_.get_ME2_value(-1, -1 , s12_and_costh_to_p(s12, costh));
//            double sz =  1; //dipole.Born_.get_ME2_value(-1, -1 , s12_and_costh_to_p(z*s12, costh));
//            double sz1 = 1; //dipole.Born_.get_ME2_value(-1, -1 , s12_and_costh_to_p(s12, costh));
//            double szij =  1; //dipole.Born_.get_ME2_value(dipole.get_spectator(), dipole.get_emitter() , s12_and_costh_to_p(z*s12, costh));
//            double sz1ij = 1; //dipole.Born_.get_ME2_value(dipole.get_spectator(), dipole.get_emitter(), s12_and_costh_to_p(s12, costh));
//            double diffij = szij - sz1ij;
//            // eq. 6.55
//            ff[0] += Alfas/(2.*pi) * (
//                     // K-bar part
//                        sz * (-CF*(1+z)*log((1-z)/z) + CF*(1-z))
//                        + CF*(2/(1-z)*log((1-z)/z))*diff
//                     // K part
//                        - (
//                           2*log(1-z)/(1-z)*diff
//                           - 2*log(2-z)/(1-z)*szij
//                           // eq. 5.58
//                     + (
//                                                                                       ((1-z)/(2*pow(1-z+muQ*muQ,2))-2/(1-z)*(1+log(1-z+muQ*muQ))) * (szij * sz1ij)
//                                   + 2/(1-z)*(log(2+muQ*muQ-z)*szij - log(1+muQ*muQ)*sz1ij)
//
//                                                                               )
//                           + 2/(1-z)*(log((2-z)*sja/((2-z)*sja+mj2))*szij - log(sja/(sja+mj2))*sz1ij)
//                           - (sja*sja*(1-z)/pow(2*(sja*(1-z)+mj2),2)*sz1ij - 0.)
//                  )
//                  -1./CF*(-CF*(1.-z)*log((1-z)*sja/((1-z)*sja+mj2)))
//            );
//         }
//      }
//   }

   double beta = sqrt(1. - 4.*1500*1500/s12);
   ff[0] *= to_fb * 1./(2.*s12) * beta/(16.*pi);
   
   // check the Jakobian
   ff[0] *= (pow(-4*pow(m1, 2) + S,2)*xx[0]*(4*pow(m1, 2)*(-1 + xx[0]) - S*(-1 + dS + xx[0])))/
            (pow(S,2)*(-4*pow(m1, 2)*(-1 + xx[0]) + S*xx[0]));
   
   return 0;
}
