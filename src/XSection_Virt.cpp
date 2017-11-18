#include <vector>
#include <numeric>
#include "XSection_Virt.hpp"
#include "../include/Vec4D.hpp"
#include "../include/IMatrixElements.h"

#include "clooptools.h"
std::vector<CSDipole> XSection_Virt::cs_dipoles;

std::vector<Vec4D<double>> mandelstam_to_p (double s, double t) {
   double mass = 5;
   double tp = t - pow(mass,2);
   double p = sqrt(s/4. - pow(mass,2));
   double costh = (s/2. + tp)/(p*sqrt(s));
   std::vector<Vec4D<double>> temp = {
           Vec4D<double> { sqrt(s)/2., 0, 0, sqrt(s)/2.},
           Vec4D<double> { sqrt(s)/2., 0, 0, -sqrt(s)/2.},
           Vec4D<double> { sqrt(s)/2., 0., p*sqrt(1-costh*costh), p*costh},
           Vec4D<double> { sqrt(s)/2., 0., -p*sqrt(1-costh*costh), -p*costh}
   };
   return temp;
}

int XSection_Virt::integrand(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

   double m1 = 5;
   double x1 = 1;
   double x2 = 1;
   double x1min = 0;
   double x2min = 0;
   if (pt.get<std::string>("collider setup.collider") == "pp"
       || pt.get<std::string>("collider setup.collider") == "ppbar") {
      double x1min = 4. * pow(m1, 2) / S;
      x1 = x1min + (1. - x1min) * xx[*ndim-2];
      double x2min = 4. * pow(m1, 2) / (S * x1);
      x2 = x2min + (1. - x2min) * xx[*ndim-1];
   }
   double s = x1 * x2 * S;     //partonic
   double Tmin = pow( m1, 2 ) - s/2. - sqrt( pow(s, 2)/4 -
                  pow( m1, 2 )*s);
   double Tmax = pow( m1, 2 ) - s/2. + sqrt( pow(s, 2)/4. -
                  pow( m1, 2 )*s);
   double T = xx[0]*(Tmax-Tmin) + Tmin;
   double jacobian = (Tmax-Tmin) * (1.-x1min) * (1.-x2min);


   // in debug mode check cancelation of double poles
   //assert(
//      std::cout << "double " <<
//         abs(
//               std::accumulate(
//                     cs_dipoles.begin(), cs_dipoles.end(), 0.,
//                     [s,T](double current,  CSDipole& el) {
//                        return current + el.eval_integrated_dipole(-2, mandelstam_to_p(s, T));
//                     }
//               ) + (model->VirtualME)(particles[0], EpsOrd::DoublePole, s, T)
//         )
//   //   < 1e-16
//      << std::endl;
   //);
   // and single poles
   assert(
//      std::cout << "single " <<
         abs(
               std::accumulate(
                     cs_dipoles.begin(), cs_dipoles.end(), 0.,
                     [s,T](double current,  CSDipole& el) {
                        return current + el.eval_integrated_dipole(-1, mandelstam_to_p(s, T));
                     }
               ) + (model->VirtualME)(particles[0], EpsOrd::SinglePole, s, T)
         )
                       // < 1e-16
//      << std::endl;
   );


   int FiniteGs = 1;
   double Dminus4 = 0;
   int Divergence = 0;     // O(eps)
double squaredMReal;
   // contraction with O(eps) from Dminus4
   Divergence = -1;           // O(eps)
   FiniteGs = 0;

   Dminus4 = -2.;
//   double squaredMRealMinus2 = (processID->*processID->matrixelementVirtual)(
//                         s, T, FiniteGs, Dminus4, Divergence);
//
//   double dSigmaPart3 = 2.*(squaredMRealMinus2 - squaredMReal)*
//                         (processID->h)*pi/(pow(4.*pi,2))/
//                         (processID->k)/(pow(s,2));

   // contraction with O(eps^2) prefactor of loop integral
   // and with product of O(eps) prefactors of phase space and loop integral
   Divergence = -2;
   Dminus4 = 0;

   // -------------------------

   FiniteGs = 1;
   Dminus4 = 0;
   Divergence = 0;     // O(eps)

   squaredMReal = (model->VirtualME)(particles[0], EpsOrd::Eps0, s, T);

   double pdf_flux = 0.0;
//   for (const auto& flav : processID->flav) {
      //pdf_flux += flav.at(2) * pdf->xfxQ( flav.at(0), x1, mu_f ) * pdf->xfxQ( flav.at(1), x2, mu_f );
//   }
   //pdf_flux /= (x1 * x2);

   double dipole_sum = std::accumulate(
           cs_dipoles.begin(), cs_dipoles.end(), 0.,
           [s,T](double current,  CSDipole& el) {
              return current + el.eval_integrated_dipole(0, mandelstam_to_p(s, T));
           }
   );
   ff[0] = (squaredMReal + 1*dipole_sum) * pi/pow(4.*pi,2)/pow(s,2) * jacobian * to_fb; //* pdf_flux;   // in femto barn
   return 0;

}


 std::array<double, 3> XSection_Virt::integrate() {
   constexpr int ndim = 2;
   constexpr int ncomp = 1;
   constexpr int nvec = 1;
   constexpr double accuracy_abs = 1e-12;
   constexpr int eval_min = 1000;
   constexpr int eval_max = 1000000;
   const int verbose = vm["verbosity-virt"].as<int>();        // adjust output 0 ... 3
   int nregions, neval, fail;
   cubareal integral[ncomp], error[ncomp], prob[ncomp];
   constexpr int last = 4;
   constexpr int key = 0;
   // Divonne specific
   constexpr int seed = 1;
   constexpr double border = 1e-9; // don't go lower than 1e-6 for 1e-4 relative accuraccy
   constexpr int key1 = 47;
   constexpr int key2 = 1;
   constexpr int key3 = 1;
   constexpr int maxpass = 5;
   constexpr double maxchisq = 10.;
   constexpr double mindeviation = .25;
   constexpr int ngiven = 0;
   constexpr int ldxgiven = ndim;
   constexpr int nextra = 0;
   const double prec_virt = pow( 10., -vm["precision-virt"].as<int>());

   static bool looptools_initialized = false;
   if(!looptools_initialized) {
      std::cout << std::endl;
      ltini();
      looptools_initialized = true;
   }
   /*
    *    LoopTools fails for phase space points near the border.
    *    That's why we have to use Divonne with option border.
    *    It's not clear if this is just a sign of numerical instability
    *    of LoopTools or something more serious.
    */
   Divonne(ndim, ncomp, integrand, NULL, nvec,
        prec_virt, accuracy_abs, verbose, seed,
        eval_min, eval_max, key1, key2, key3, maxpass,
        border, maxchisq, mindeviation,
        ngiven, ldxgiven, NULL, nextra, NULL,
        NULL, NULL,
        &nregions, &neval, &fail, integral, error, prob);

   std::array <double, 3> result{ integral[0], error[0], prob[0] };

   return result;
}
