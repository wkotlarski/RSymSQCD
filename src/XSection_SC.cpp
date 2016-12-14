#include "XSection_SC.hpp"

std::array<double, 3> XSection_SC::integrate() {
  
   // integrals dimensions, number of integrands
   constexpr int ndim { 3 }, ncomp { 1 };
   // accuraccy
   constexpr double accuracy_rel_sc { 1e-6 }, 
            accuracy_rel_c { 1e-6 };
   constexpr double accuracy_abs { 1e-12 };

   constexpr int neval_min = 10000;
   long long int neval;
   constexpr long long int neval_max { 100'000'000 }; 

   int nregions, fail;

   cubareal integral_sc[ncomp], error_sc[ncomp], prob_sc[ncomp];
   llCuhre(ndim, ncomp, integrand_sc, NULL, 1,
      prec_sc, accuracy_abs, 0,
      neval_min, neval_max, 1, NULL, NULL,
      &nregions, &neval, &fail, integral_sc, error_sc, prob_sc);
   
   cubareal integral_c1[ncomp], error_c1[ncomp], prob_c1[ncomp];
   llCuhre(ndim, ncomp, integrand_c1, NULL, 1,
      accuracy_rel_c, accuracy_abs, 0,
      neval_min, neval_max, 1, NULL, NULL,
      &nregions, &neval, &fail, integral_c1, error_c1, prob_c1);
  
   cubareal integral_c2[ncomp], error_c2[ncomp], prob_c2[ncomp];
   llCuhre(ndim, ncomp, integrand_c2, NULL, 1,
      accuracy_rel_c, accuracy_abs, 0,
      neval_min, neval_max, 1, NULL, NULL,
      &nregions, &neval, &fail, integral_c2, error_c2, prob_c2);

   std::array <double, 3> result_finite { 
      integral_sc[0]  + integral_c1[0] + integral_c2[0], 
      sqrt( pow( error_sc[0], 2)  + pow ( error_c1[0], 2) + pow ( error_c2[0], 2) ),
      prob_sc[0]  + prob_c2[0] 
   };
  
  return result_finite;
}

int XSection_SC::integrand_sc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
   
   // double m_sqr = pow( pt.get<double>("masses.suL"), 2);
    
   // integration variables
   double x1 = 4. * pow(m1, 2)/S + (1 - 4. * pow(m1, 2)/S ) * xx[0];
	double x2 = 4. * pow(m1, 2) /(S * x1) + (1 - 4. * pow(m1, 2)/(S * x1)) * xx[1];
   double th = xx[2] * pi;
    
   double s12 = x1 * x2 * S;
   double Alfas = pdf->alphasQ( mu_r );
    
   double pdf_flux = 0.0;
   for (const auto& inner : processID->flav) {
      pdf_flux += inner.at(2) * pdf->xfxQ( inner.at(0), x1, mu_f ) * pdf->xfxQ( inner.at(1), x2, mu_f );
   }
   
   ff[0] = to_fb * (processID->*processID->matrixelementReal_SC)(s12, th) 
            * pdf_flux/( x1 * x2 );
    
   // jakobian
   ff[0] *= pi*Power(-4*pow(m1, 2) + S,2)*xx[0] / 
           (S*(-4*pow(m1, 2)*(-1 + xx[0]) + S*xx[0]));
   return 0;
}

// @todo if muR != muF one needs one more term here

int XSection_SC::integrand_c1(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {
    
   // integration variables
   double x1 = 4. * pow(m1, 2)/S + (1 - 4. * pow(m1, 2)/S ) * xx[0];
   double x2 = 4. * pow(m1, 2) /(S * x1) + (1 - 4. * pow(m1, 2)/(S * x1)) * xx[1];
   double th = pi * xx[2];
    
   double s12 = x1 * x2 * S;
   double beta = sqrt( 1 - 4 * pow(m1, 2)/s12 );
    
   double Alfas = pdf->alphasQ( muR );
   double Alfas2 = pow( Alfas, 2);
    
   double pdf_flux = 0.0;
   for (const auto& inner : processID->flav) {
      pdf_flux += inner.at(2) * pdf->xfxQ( inner.at(0), x1, mu_f ) * pdf->xfxQ( inner.at(1), x2, mu_f );
   }
   
   ff[0] = 0*to_fb
      * pdf->xfxQ(2, x1, muF)/x1 * pdf->xfxQ(2, x2, muF)/x2
      * 4./3. * (2 * log(dS) + 3./2.);
    
   ff[0] *= (pi*Power(-4*pow(m1, 2) + S,2)*xx[0])/(S*(-4*pow(m1, 2)*(-1 + xx[0]) + S*xx[0]));
    
   return 0;
}  

int XSection_SC::integrand_c2(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {
   
   // integration variables
   double x1 = 4. * pow(m1, 2)/S + (1 - 4. * pow(m1, 2)/S ) * xx[0];
	double x2 = 4. * pow(m1, 2) /(S * x1) + (1 - 4. * pow(m1, 2)/(S * x1)) * xx[1];
   double y = x1 + ( 1 - dS - x1 ) * xx[2];
     
   double s12 = x1 * x2 * S;
   double beta = sqrt( 1 - 4 * pow(m1, 2)/s12 );
    
   double Alfas = pdf->alphasQ( mu_r );
   double Alfas2 = pow( Alfas, 2);
       
   double pdf_flux = 0.0;
   for (const auto& inner : processID->flav) {
      if( abs(inner.at(0)) == 1 || abs(inner.at(0)) == 2 || abs(inner.at(0)) == 3 || abs(inner.at(0)) == 4
              || abs(inner.at(0)) == 5 ) {
         pdf_flux += inner.at(2) * pdf->xfxQ( inner.at(0), std::min(x1/y, 1.), mu_f )/std::min(x1/y, 1.) 
                  * pdf->xfxQ( inner.at(1), x2, mu_f )/x2 * CF *
                 ( (1 - y) + (1 + y*y)/(1 - y) * log( dC * s12 * pow(1 - y, 2) / (2 * mu_f * mu_f * y) ) );
      }
      else if( inner.at(0) == 0 ) {
         pdf_flux += inner.at(2) * pdf->xfxQ( inner.at(0), std::min(x1/y, 1.), mu_f )/std::min(x1/y, 1.) 
                  * pdf->xfxQ( inner.at(1), x2, mu_f )/x2 * 2 * CA *
                 (y/(1-y) + (1-y)/y + y*(1-y)) * log( dC * s12 * pow(1 - y, 2) / (2 * mu_f * mu_f * y) );
      }
      if( abs(inner.at(1)) == 1 || abs(inner.at(1)) == 2 || abs(inner.at(1)) == 3 || abs(inner.at(1)) == 4
              || abs(inner.at(1)) == 5 ) {
         pdf_flux += inner.at(2) * pdf->xfxQ( inner.at(0), x1, mu_f )/x1
                  * pdf->xfxQ( inner.at(1), std::min(x2/y, 1.), mu_f )/std::min(x2/y, 1.) * CF *
                 ( (1 - y) + (1 + y*y)/(1 - y) * log( dC * s12 * pow(1 - y, 2) / (2 * mu_f * mu_f * y) ) );
      }
            else if( inner.at(1) == 0 ) {
         pdf_flux += inner.at(2) * pdf->xfxQ( inner.at(0), x1, mu_f )/x1
                  * pdf->xfxQ( inner.at(1), std::min(x2/y, 1.), mu_f )/std::min(x2/y, 1.) * 2 * CA *
                 (y/(1-y) + (1-y)/y + y*(1-y)) * log( dC * s12 * pow(1 - y, 2) / (2 * mu_f * mu_f * y) );
      }
   }
      
   ff[0] = to_fb * 1./y * pdf_flux 
            * Alfas/two_pi * (processID->*processID->sigmaPartTree)(s12);
   
    // multiply by jakobian of integration variable transformation
    ff[0] *= (Power(-4*pow(m1, 2) + S,2)*xx[0]*(4*pow(m1, 2)*(-1 + xx[0]) - S*(-1 + dS + xx[0])))/
        (Power(S,2)*(-4*pow(m1, 2)*(-1 + xx[0]) + S*xx[0]));
    
  return 0;
}
