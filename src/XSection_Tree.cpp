#include "XSection_Tree.hpp"
#include "constants.hpp"

/*
 * E1 = sqrtS/2 (1 + (m1^2 - m2^2)/s)
 * E2 = sqrtS/2 (1 + (m2^2 - m1^2)/s)
 * b  = sqrt(1 - 2(m1^2 + m2^2)/s + ((m1 + m2)(m1-m2)/s)^2)
 * p  = sqrtS/2 b
 *
 */
int XSection_Tree::integrand(const int *ndim, const cubareal xx[],
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

   double result = 0.0;
   for (const auto& v : particles) {
      double temp = (model->BornME)(v, EpsOrd::Eps0, s, T);
      if (
            v[0] == Particle::u ||
            v[0] == Particle::ubar || 
            v[0] == Particle::d ||
            v[0] == Particle::dbar ||
            v[0] == Particle::s ||
            v[0] == Particle::sbar ||
            v[0] == Particle::c ||
            v[0] == Particle::cbar ||
            v[0] == Particle::b ||
            v[0] == Particle::bbar) temp *= pdf->xfxQ( particle_to_int(v[0]), x1, mu_f ); 
      if (
            v[1] == Particle::u ||
            v[1] == Particle::ubar || 
            v[1] == Particle::d ||
            v[1] == Particle::dbar ||
            v[1] == Particle::s ||
            v[1] == Particle::sbar ||
            v[1] == Particle::c ||
            v[1] == Particle::cbar ||
            v[1] == Particle::b ||
            v[1] == Particle::bbar) temp *= pdf->xfxQ( particle_to_int(v[1]), x2, mu_f ); 
      result += temp;
   }
   result /= (x1 * x2);
    
   result *= M_PI/(pow(4.*pi,2))/(pow(s,2));
           
   ff[0] = result * jacobian * to_fb;

   return 0;
}

 std::array<double, 3> XSection_Tree::integrate() {
   constexpr int ndim = 2;
   constexpr int ncomp = 1;
   constexpr int nvec = 1;
   constexpr double accuracy_rel = 1e-6;
   constexpr double accuracy_abs = 1e-12;
   constexpr int eval_min = 1000;
   constexpr int eval_max = 1000000000;
   const int verbose = vm["verbosity-born"].as<int>();        // adjust shown output 0 ... 3
   constexpr int last = 4;
   constexpr int key = 0;
   int nregions, neval, fail;
   cubareal integral[ncomp], error[ncomp], prob[ncomp];

   Cuhre( ndim, ncomp, integrand, NULL, nvec,
      accuracy_rel, accuracy_abs, verbose | last,
      eval_min, eval_max, key, NULL, NULL,
      &nregions, &neval, &fail, integral, error, prob );

   std::array <double, 3> result{ integral[0], error[0], prob[0] };

   return result;
}
