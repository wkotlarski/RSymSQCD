#include "XSection_Tree.hpp"

/*
 * E1 = sqrtS/2 (1 + (m1^2 - m2^2)/s)
 * E2 = sqrtS/2 (1 + (m2^2 - m1^2)/s)
 * b  = sqrt(1 - 2(m1^2 + m2^2)/s + ((m1 + m2)(m1-m2)/s)^2)
 * p  = sqrtS/2 b
 *
 */
int XSection_Tree::integrand(const int *ndim, const cubareal xx[],
    const int *ncomp, cubareal ff[], void *userdata) {

    double mu_f = pt.get<double>("collider setup.mu_f");
    double x1min = 4. * pow(processID->m1,2)/S;
    double xmax = 1.;
    double x1 = x1min + (xmax - x1min ) * xx[0];
    double x2min = 4. * pow(processID->m1,2)/(S*x1);
    double x2 = x2min + (xmax - x2min) * xx[1];
    double s = S * x1 * x2;     //partonic 

    /* integration of |M^B|^2 */
    if(processID->partonic == false) {		
        double Tmin = pow( processID->m1, 2 ) - s/2. - sqrt( pow(s, 2)/4 -
                      pow( processID->m1, 2 )*s);
        double Tmax = pow( processID->m1, 2 ) - s/2. + sqrt( pow(s, 2)/4. -
                      pow( processID->m1, 2 )*s);
        double T = xx[2]*(Tmax-Tmin) + Tmin;
        double jacobian = (Tmax-Tmin)*(xmax-x1min)*(xmax-x2min);
        double squaredM = (processID->*processID->matrixelementTree)(
                           s, T);
        double dSigmaPart = squaredM*(processID->h)*M_PI/(pow(4.*M_PI,2))/
                         (processID->k)/(pow(s,2));
        double pdf_flux = 0.0;
        for (const auto& inner : processID->flav) {
               //std::cout << inner.at(0) << " " << inner.at(1) << '\n';
         pdf_flux += inner.at(2) * pdf->xfxQ( inner.at(0), x1, mu_f ) * pdf->xfxQ( inner.at(1), x2, mu_f );
         }            
        ff[0] = dSigmaPart * pdf_flux * jacobian*to_fb / (x1 * x2); 
    }
    
    /* integration of partonic cross section */
    else {    
    double pdf_flux = 0.;
    if( processID->f1 == 69 && processID->f2 == 69 ) {
        for ( int i = 1; i < 6; ++i) pdf_flux += 2. * pdf->xfxQ( i, x1, mu_f ) *
             pdf->xfxQ( -i, x2, mu_f );
    } else if( processID->f1 == 0 && processID->f2 == 0 ) {
        pdf_flux = pdf->xfxQ( processID->f1, x1, mu_f ) *
        pdf->xfxQ( processID->f2, x2, mu_f );
    } else {
        pdf_flux = pdf->xfxQ( processID->f1, x1, mu_f ) *
        pdf->xfxQ( processID->f2, x2, mu_f );
    }
    pdf_flux /= x1 * x2;
   
    ff[0] = (processID->*processID->sigmaPartTree)(s) * to_fb * pdf_flux *
       pow(-4.*pow(m1, 2) + S, 2)*xx[0] / 
       (S*(-4*pow(m1, 2)*(-1 + xx[0]) + S*xx[0]));
    }   
     
    return 0;
}


 std::array<double, 3> XSection_Tree::integrate() {
    constexpr int ndim = 3;
    constexpr int ncomp = 1;
    constexpr int nvec = 1;
    constexpr double accuracy_rel = 1e-4;
    constexpr double accuracy_abs = 1e-12;
    constexpr int eval_min = 1000;
    constexpr int eval_max = 1000000;
    constexpr int verbose = 0;        // adjust shown output 0 ... 3
    constexpr int last = 4;
    constexpr int key = 0;
    int nregions, neval, fail;
    cubareal integral[ncomp], error[ncomp], prob[ncomp];

      Cuhre(ndim, ncomp, integrand, NULL, nvec,
      accuracy_rel, accuracy_abs, verbose | last,
      eval_min, eval_max, key, NULL, NULL,
      &nregions, &neval, &fail, integral, error, prob);

    std::array <double, 3> result{integral[0], error[0], prob[0]}; 

    return result;
}
