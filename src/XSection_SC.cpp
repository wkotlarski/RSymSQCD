#include "XSection_SC.hpp"

std::array<double, 3> XSection_SC::integrate() {
  
  //  integral dimension, number of integrands
  constexpr int ndim { 3 }, ncomp { 1 };
  //  accuraccy
    constexpr double accuracy_rel_sc { 1e-3 }, 
            accuracy_rel_c { 1e-3 };
    constexpr double accuracy_abs { 1e-12 };

  constexpr int neval_min = 10000;
  long long int neval;
  constexpr long long int neval_max { 100'000'000 }; 
    // @TODO: read from external source strtoll( "1e+3", NULL, 10 );

  // technical (Vegas specific) stuff
  constexpr int nstart      = 200000;
  constexpr int nincrease   = 100;
  constexpr int nbatch      = 1000;
  constexpr int gridno      = 0;
  const char* state_file    = "";
  int nregions, fail;

  cubareal integral_sc[ncomp], error_sc[ncomp], prob_sc[ncomp];
  llVegas( ndim, ncomp, integrand_sc, NULL, 1,
           accuracy_rel_sc, accuracy_abs, 8 | 1, 0,
           neval_min, neval_max, nstart, nincrease, nbatch,
           gridno, state_file, NULL,
           &neval, &fail, integral_sc, error_sc, prob_sc );

  cubareal integral_c1[ncomp], error_c1[ncomp], prob_c1[ncomp];
  llVegas( ndim, ncomp, integrand_c1, NULL, 1,
           accuracy_rel_c, accuracy_abs, 8 | 1, 0,
           neval_min, neval_max, nstart, nincrease, nbatch,
           gridno, state_file, NULL,
           &neval, &fail, integral_c1, error_c1, prob_c1 );  
  
  cubareal integral_c2[ncomp], error_c2[ncomp], prob_c2[ncomp];
  llVegas( ndim, ncomp, integrand_c2, NULL, 1,
           accuracy_rel_c, accuracy_abs, 8 | 1, 0,
           neval_min, neval_max, nstart, nincrease, nbatch,
           gridno, state_file, NULL,
           &neval, &fail, integral_c2, error_c2, prob_c2 );
  std::cout << integral_sc[0] << ' ' << integral_c1[0] << ' ' <<  integral_c2[0] << '\n';
  std::array <double, 3> result_finite { 
      integral_sc[0] + integral_c1[0] + integral_c2[0], 
      sqrt( pow( error_sc[0], 2) + pow ( error_c1[0], 2) + pow ( error_c2[0], 2) ),
      prob_sc[0] + prob_c1[0] + prob_c2[0] 
  };
  
  return result_finite;
}

int XSection_SC::integrand_sc(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
    
    double m_sqr = pow( 1500., 2);
    double MGl2 = pow ( gluino_mass, 2 );
    double muR = 1500.;
    
    // integration variables
    double x1 = 4. * m_sqr/S + (1 - 4. * m_sqr/S ) * xx[0];
	double x2 = 4. * m_sqr /(S * x1) + (1 - 4. * m_sqr/(S * x1)) * xx[1];
    double th = xx[2] * pi;
    
    double s12 = x1 * x2 * S;
    double beta = sqrt( 1 - 4. * m_sqr/s12 );
    
    double Alfas = pdf_nlo->alphasQ(1500.);
    double Alfas2 = pow( Alfas, 2);
    
    std::complex<double> sigma = ((-42.666666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (42.666666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*MGl2*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
   (2.3703703703703704*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) - 
   (0.59259259259259259*Alfas*Alfas2*beta*(-1. + beta*beta)*(s12*s12)*(-1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
   (1.1851851851851852*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
   (0.59259259259259259*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (1.1851851851851852*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (42.666666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (42.666666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (21.333333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) - 
   (10.666666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.14814814814814815*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.59259259259259259*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.14814814814814815*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + 
        beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.59259259259259259*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + 
        beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.59259259259259259*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.14814814814814815*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + 
        beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.074074074074074074*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.14814814814814815*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.074074074074074074*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.14814814814814815*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (2.3703703703703704*Alfas*Alfas2*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (2.3703703703703704*Alfas*Alfas2*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.14814814814814815*Alfas*Alfas2*Power(-1. + beta*beta,4)*Power(s12,3)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.2962962962962963*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.2962962962962963*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.2962962962962963*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.2962962962962963*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*Power(1. + beta*Cos(th),2)*
      (-2.0794415416798359*beta*Sin(th) + beta*euler*Sin(th) + beta*Log(beta)*Sin(th) + beta*Log(dS)*Sin(th) - 1.*beta*Log(pi)*Sin(th) - 
        1.*beta*Log((muR*muR)/s12)*Sin(th) + beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 
        1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 
        2.*beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 
        1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 
        2.*beta*Log(Sin(th))*Sin(th)))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (9.4814814814814815*Alfas*Alfas2*(-1. + beta*beta)*(MGl2*MGl2)*s12*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (9.4814814814814815*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (4.7407407407407407*Alfas*Alfas2*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.59259259259259259*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (9.4814814814814815*Alfas*Alfas2*(MGl2*MGl2)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (9.4814814814814815*Alfas*Alfas2*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (4.7407407407407407*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.59259259259259259*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(s12,3)*Power(-1. + beta*Cos(th),3)*(1. + beta*Cos(th))*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),2)*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.1851851851851852*Alfas*Alfas2*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),3)*
      (-4.1588830833596719*beta*Sin(th) + 2.*beta*euler*Sin(th) + 2.*beta*Log(beta)*Sin(th) - 1.*Log((1. + beta)/(1. - 1.*beta))*Sin(th) + 
        2.*beta*Log(dS)*Sin(th) - 2.*beta*Log(pi)*Sin(th) - 2.*beta*Log((muR*muR)/s12)*Sin(th) + 2.*beta*Log(Sin(th))*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*(0.25*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.5*beta*Power(Log(beta),2) + 0.5*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(1.039720770839918 + 0.5*Log(pi)) + beta*(2.1620385626319064 + 2.0794415416798359*Log(pi) + 0.5*Power(Log(pi),2)) + 
           2.*(-0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)))*Log((muR*muR)/s12) + 0.5*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.5*beta*euler - 0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12)))*Sin(th) + 
        2.*(-0.5*beta*euler - 0.5*beta*Log(beta) - 0.5*beta*Log(dS) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) - 0.5*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(-1.*
         (0.25*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.5*beta*Power(Log(beta),2) + 0.5*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(1.039720770839918 + 0.5*Log(pi)) + beta*(2.1620385626319064 + 2.0794415416798359*Log(pi) + 0.5*Power(Log(pi),2)) + 
           2.*(-0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)))*Log((muR*muR)/s12) + 0.5*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.5*beta*euler - 0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12)))*Sin(th) + 
        2.*(-0.5*beta*euler - 0.5*beta*Log(beta) - 0.5*beta*Log(dS) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) - 0.5*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. - 1.*beta*Cos(th),3)*
      (-1.*(0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.25*beta*Power(Log(beta),2) + 0.25*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.51986038541995898 + 0.25*Log(pi)) + beta*(1.0810192813159532 + 1.039720770839918*Log(pi) + 0.25*Power(Log(pi),2)) + 
           2.*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)))*Log((muR*muR)/s12) + 0.25*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.25*beta*euler - 0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)))*Sin(th) + 
        2.*(-0.25*beta*euler - 0.25*beta*Log(beta) - 0.25*beta*Log(dS) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) - 0.25*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(1. - 1.*beta*Cos(th))*
      (-1.*(0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.25*beta*Power(Log(beta),2) + 0.25*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.51986038541995898 + 0.25*Log(pi)) + beta*(1.0810192813159532 + 1.039720770839918*Log(pi) + 0.25*Power(Log(pi),2)) + 
           2.*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)))*Log((muR*muR)/s12) + 0.25*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.25*beta*euler - 0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)))*Sin(th) + 
        2.*(-0.25*beta*euler - 0.25*beta*Log(beta) - 0.25*beta*Log(dS) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) - 0.25*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*(0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.125*beta*Power(Log(beta),2) + 0.125*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.25993019270997949 + 0.125*Log(pi)) + beta*(0.5405096406579766 + 0.51986038541995898*Log(pi) + 0.125*Power(Log(pi),2)) + 
           2.*(-0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)))*Log((muR*muR)/s12) + 0.125*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.125*beta*euler - 0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12)))*Sin(th) + 
        2.*(-0.125*beta*euler - 0.125*beta*Log(beta) - 0.125*beta*Log(dS) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12))*
         Log(Sin(th))*Sin(th) - 0.125*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)\
    + (0.037037037037037037*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-1.*(0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.125*beta*Power(Log(beta),2) + 0.125*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.25993019270997949 + 0.125*Log(pi)) + beta*(0.5405096406579766 + 0.51986038541995898*Log(pi) + 0.125*Power(Log(pi),2)) + 
           2.*(-0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)))*Log((muR*muR)/s12) + 0.125*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.125*beta*euler - 0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12)))*Sin(th) + 
        2.*(-0.125*beta*euler - 0.125*beta*Log(beta) - 0.125*beta*Log(dS) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12))*
         Log(Sin(th))*Sin(th) - 0.125*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)\
    + (0.037037037037037037*Alfas*Alfas2*s12*Power(1. - 1.*beta*Cos(th),2)*(1. + beta*Cos(th))*
      ((0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.125*beta*Power(Log(beta),2) + 0.125*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.25993019270997949 + 0.125*Log(pi)) + beta*(0.5405096406579766 + 0.51986038541995898*Log(pi) + 0.125*Power(Log(pi),2)) + 
           2.*(-0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)))*Log((muR*muR)/s12) + 0.125*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.125*beta*euler - 0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.125*beta*euler - 0.125*beta*Log(beta) - 0.125*beta*Log(dS) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12))*
         Log(Sin(th))*Sin(th) + 0.125*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)\
    + (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*Power(1. + beta*Cos(th),2)*
      ((0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.125*beta*Power(Log(beta),2) + 0.125*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.25993019270997949 + 0.125*Log(pi)) + beta*(0.5405096406579766 + 0.51986038541995898*Log(pi) + 0.125*Power(Log(pi),2)) + 
           2.*(-0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)))*Log((muR*muR)/s12) + 0.125*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.125*beta*euler - 0.125*beta*Log(beta) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.125*beta*euler - 0.125*beta*Log(beta) - 0.125*beta*Log(dS) + beta*(0.25993019270997949 + 0.125*Log(pi)) + 0.125*beta*Log((muR*muR)/s12))*
         Log(Sin(th))*Sin(th) + 0.125*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)\
    + (0.037037037037037037*Alfas*Alfas2*Power(1. - 1.*(beta*beta),2)*s12*
      ((0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.25*beta*Power(Log(beta),2) + 0.25*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.51986038541995898 + 0.25*Log(pi)) + beta*(1.0810192813159532 + 1.039720770839918*Log(pi) + 0.25*Power(Log(pi),2)) + 
           2.*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)))*Log((muR*muR)/s12) + 0.25*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.25*beta*euler - 0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.25*beta*euler - 0.25*beta*Log(beta) - 0.25*beta*Log(dS) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) + 0.25*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      ((0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.25*beta*Power(Log(beta),2) + 0.25*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.51986038541995898 + 0.25*Log(pi)) + beta*(1.0810192813159532 + 1.039720770839918*Log(pi) + 0.25*Power(Log(pi),2)) + 
           2.*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)))*Log((muR*muR)/s12) + 0.25*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.25*beta*euler - 0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.25*beta*euler - 0.25*beta*Log(beta) - 0.25*beta*Log(dS) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) + 0.25*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. - 1.*beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      ((0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.25*beta*Power(Log(beta),2) + 0.25*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.51986038541995898 + 0.25*Log(pi)) + beta*(1.0810192813159532 + 1.039720770839918*Log(pi) + 0.25*Power(Log(pi),2)) + 
           2.*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)))*Log((muR*muR)/s12) + 0.25*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.25*beta*euler - 0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.25*beta*euler - 0.25*beta*Log(beta) - 0.25*beta*Log(dS) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) + 0.25*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*Power(1. - 1.*(beta*beta),2)*s12*
      ((0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.25*beta*Power(Log(beta),2) + 0.25*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.51986038541995898 + 0.25*Log(pi)) + beta*(1.0810192813159532 + 1.039720770839918*Log(pi) + 0.25*Power(Log(pi),2)) + 
           2.*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)))*Log((muR*muR)/s12) + 0.25*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.25*beta*euler - 0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.25*beta*euler - 0.25*beta*Log(beta) - 0.25*beta*Log(dS) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) + 0.25*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. - 1.*beta*Cos(th),2)*(1. + beta*Cos(th))*
      ((0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.25*beta*Power(Log(beta),2) + 0.25*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.51986038541995898 + 0.25*Log(pi)) + beta*(1.0810192813159532 + 1.039720770839918*Log(pi) + 0.25*Power(Log(pi),2)) + 
           2.*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)))*Log((muR*muR)/s12) + 0.25*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.25*beta*euler - 0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.25*beta*euler - 0.25*beta*Log(beta) - 0.25*beta*Log(dS) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) + 0.25*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. - 1.*beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      ((0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.25*beta*Power(Log(beta),2) + 0.25*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(0.51986038541995898 + 0.25*Log(pi)) + beta*(1.0810192813159532 + 1.039720770839918*Log(pi) + 0.25*Power(Log(pi),2)) + 
           2.*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)))*Log((muR*muR)/s12) + 0.25*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.25*beta*euler - 0.25*beta*Log(beta) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.25*beta*euler - 0.25*beta*Log(beta) - 0.25*beta*Log(dS) + beta*(0.51986038541995898 + 0.25*Log(pi)) + 0.25*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) + 0.25*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. - 1.*beta*Cos(th),2)*
      ((0.25*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.5*beta*Power(Log(beta),2) + 0.5*beta*Power(Log(dS),2) - 
           2.*beta*Log(beta)*(1.039720770839918 + 0.5*Log(pi)) + beta*(2.1620385626319064 + 2.0794415416798359*Log(pi) + 0.5*Power(Log(pi),2)) + 
           2.*(-0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)))*Log((muR*muR)/s12) + 0.5*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.5*beta*euler - 0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.5*beta*euler - 0.5*beta*Log(beta) - 0.5*beta*Log(dS) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) + 0.5*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. + beta*Cos(th))*((0.25*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi)) + 0.5*beta*Power(Log(beta),2) + 
           0.5*beta*Power(Log(dS),2) - 2.*beta*Log(beta)*(1.039720770839918 + 0.5*Log(pi)) + 
           beta*(2.1620385626319064 + 2.0794415416798359*Log(pi) + 0.5*Power(Log(pi),2)) + 
           2.*(-0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)))*Log((muR*muR)/s12) + 0.5*beta*Power(Log((muR*muR)/s12),2) - 
           2.*euler*(-0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12)) - 
           2.*Log(dS)*(-0.5*beta*euler - 0.5*beta*Log(beta) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12)))*Sin(th) - 
        2.*(-0.5*beta*euler - 0.5*beta*Log(beta) - 0.5*beta*Log(dS) + beta*(1.039720770839918 + 0.5*Log(pi)) + 0.5*beta*Log((muR*muR)/s12))*Log(Sin(th))*
         Sin(th) + 0.5*beta*Power(Log(Sin(th)),2)*Sin(th)))/Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.79012345679012346*Alfas*Alfas2*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.79012345679012346*Alfas*Alfas2*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.39506172839506173*Alfas*Alfas2*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.049382716049382716*Alfas*Alfas2*Power(-1. + beta*beta,4)*Power(s12,3)*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.098765432098765432*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.098765432098765432*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.79012345679012346*Alfas*Alfas2*(-1. + beta*beta)*(MGl2*MGl2)*s12*(1. + beta*Cos(th))*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.79012345679012346*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.39506172839506173*Alfas*Alfas2*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.049382716049382716*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.098765432098765432*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.098765432098765432*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.098765432098765432*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),3)*
      (25.944462751582877*beta*Sin(th) - 24.953298500158031*beta*euler*Sin(th) + 6.*beta*(euler*euler)*Sin(th) - 1.*beta*(pi*pi)*Sin(th) - 
        24.953298500158031*beta*Log(beta)*Sin(th) + 12.*beta*euler*Log(beta)*Sin(th) + 6.*beta*Power(Log(beta),2)*Sin(th) - 
        24.953298500158031*beta*Log(dS)*Sin(th) + 12.*beta*euler*Log(dS)*Sin(th) + 12.*beta*Log(beta)*Log(dS)*Sin(th) + 6.*beta*Power(Log(dS),2)*Sin(th) + 
        24.953298500158031*beta*Log(pi)*Sin(th) - 12.*beta*euler*Log(pi)*Sin(th) - 12.*beta*Log(beta)*Log(pi)*Sin(th) - 12.*beta*Log(dS)*Log(pi)*Sin(th) + 
        6.*beta*Power(Log(pi),2)*Sin(th) + 24.953298500158031*beta*Log((muR*muR)/s12)*Sin(th) - 12.*beta*euler*Log((muR*muR)/s12)*Sin(th) - 
        12.*beta*Log(beta)*Log((muR*muR)/s12)*Sin(th) - 12.*beta*Log(dS)*Log((muR*muR)/s12)*Sin(th) + 12.*beta*Log(pi)*Log((muR*muR)/s12)*Sin(th) + 
        6.*beta*Power(Log((muR*muR)/s12),2)*Sin(th) - 24.953298500158031*beta*Log(Sin(th))*Sin(th) + 12.*beta*euler*Log(Sin(th))*Sin(th) + 
        12.*beta*Log(beta)*Log(Sin(th))*Sin(th) + 12.*beta*Log(dS)*Log(Sin(th))*Sin(th) - 12.*beta*Log(pi)*Log(Sin(th))*Sin(th) - 
        12.*beta*Log((muR*muR)/s12)*Log(Sin(th))*Sin(th) + 6.*beta*Power(Log(Sin(th)),2)*Sin(th)))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. - 1.*beta*Cos(th),3)*
      (-0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. - 1.*beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. - 1.*beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. - 1.*beta*Cos(th),3)*
      (-0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. - 1.*beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. - 1.*beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*Power(1. - 1.*beta*Cos(th),2)*
      (0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*Power(1. + beta*Cos(th),2)*
      (0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*Power(1. - 1.*beta*Cos(th),2)*
      (0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*Power(1. + beta*Cos(th),2)*
      (0.125*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(0.51986038541995898 - 0.25*euler - 0.25*Log(beta) + 0.25*Log(pi) + 0.25*Log((muR*muR)/s12)) + 
           0.125*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(1. - 1.*beta*Cos(th))*
      (-0.25*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(1.039720770839918 - 0.5*euler - 0.5*Log(beta) + 0.5*Log(pi) + 0.5*Log((muR*muR)/s12)) + 
           0.25*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(1. + beta*Cos(th))*
      (-0.25*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(1.039720770839918 - 0.5*euler - 0.5*Log(beta) + 0.5*Log(pi) + 0.5*Log((muR*muR)/s12)) + 
           0.25*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(1. - 1.*beta*Cos(th))*
      (-0.25*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(1.039720770839918 - 0.5*euler - 0.5*Log(beta) + 0.5*Log(pi) + 0.5*Log((muR*muR)/s12)) + 
           0.25*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(1. + beta*Cos(th))*
      (-0.25*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(1.039720770839918 - 0.5*euler - 0.5*Log(beta) + 0.5*Log(pi) + 0.5*Log((muR*muR)/s12)) + 
           0.25*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*Power(1. - 1.*(beta*beta),2)*s12*
      (0.25*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(1.039720770839918 - 0.5*euler - 0.5*Log(beta) + 0.5*Log(pi) + 0.5*Log((muR*muR)/s12)) + 
           0.25*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*Power(1. - 1.*(beta*beta),2)*s12*
      (0.25*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(1.039720770839918 - 0.5*euler - 0.5*Log(beta) + 0.5*Log(pi) + 0.5*Log((muR*muR)/s12)) + 
           0.25*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(-0.5*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(2.0794415416798359 - 1.*euler - 1.*Log(beta) + Log(pi) + Log((muR*muR)/s12)) + 
           0.5*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*s12*(-0.5*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) - 
        1.*(-0.5*Log((1. + beta)/(1. - 1.*beta))*(2.0794415416798359 - 1.*euler - 1.*Log(beta) + Log(pi) + Log((muR*muR)/s12)) + 
           0.5*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. - 1.*beta*Cos(th))*
      (0.5*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(2.0794415416798359 - 1.*euler - 1.*Log(beta) + Log(pi) + Log((muR*muR)/s12)) + 
           0.5*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(0.5*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(2.0794415416798359 - 1.*euler - 1.*Log(beta) + Log(pi) + Log((muR*muR)/s12)) + 
           0.5*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. - 1.*beta*Cos(th))*
      (0.5*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(2.0794415416798359 - 1.*euler - 1.*Log(beta) + Log(pi) + Log((muR*muR)/s12)) + 
           0.5*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(0.5*Log((1. + beta)/(1. - 1.*beta))*Log(Sin(th))*Sin(th) + 
        (-0.5*Log((1. + beta)/(1. - 1.*beta))*(2.0794415416798359 - 1.*euler - 1.*Log(beta) + Log(pi) + Log((muR*muR)/s12)) + 
           0.5*(-0.25*Power(Log((1. + beta)/(1. - 1.*beta)),2) + Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 1.*PolyLog(2.,(2.*beta)/(1. + beta))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-0.03125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.03125*beta*euler*s12 + 0.03125*beta*s12*Log(beta) - 0.03125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) + 
           0.015625*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.015625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.03125*beta*s12*Power(Log(beta),2) - 
           0.03125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) - 
           1.*beta*(s12*(0.13512741016449415 - 0.12996509635498975*Log(dS) + 0.03125*Power(Log(dS),2) + 2.*(0.064982548177494873 - 0.03125*Log(dS))*Log(pi) + 
                 0.03125*Power(Log(pi),2) - 1.*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi))*Log(1/s12) + 0.0078125*Power(Log(1/s12),2)) - 
              1.*s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12))*Log(s12) + 0.0078125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.03125*beta*s12*Log(beta) + beta*
               (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           2.*euler*(-0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           (-0.03125*beta*euler*s12 - 0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.015625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
      (-0.03125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.03125*beta*euler*s12 + 0.03125*beta*s12*Log(beta) - 0.03125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) + 
           0.015625*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.015625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.03125*beta*s12*Power(Log(beta),2) - 
           0.03125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) - 
           1.*beta*(s12*(0.13512741016449415 - 0.12996509635498975*Log(dS) + 0.03125*Power(Log(dS),2) + 2.*(0.064982548177494873 - 0.03125*Log(dS))*Log(pi) + 
                 0.03125*Power(Log(pi),2) - 1.*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi))*Log(1/s12) + 0.0078125*Power(Log(1/s12),2)) - 
              1.*s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12))*Log(s12) + 0.0078125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.03125*beta*s12*Log(beta) + beta*
               (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           2.*euler*(-0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           (-0.03125*beta*euler*s12 - 0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.015625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-0.03125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.03125*beta*euler*s12 + 0.03125*beta*s12*Log(beta) - 0.03125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) + 
           0.015625*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.015625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.03125*beta*s12*Power(Log(beta),2) - 
           0.03125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) - 
           1.*beta*(s12*(0.13512741016449415 - 0.12996509635498975*Log(dS) + 0.03125*Power(Log(dS),2) + 2.*(0.064982548177494873 - 0.03125*Log(dS))*Log(pi) + 
                 0.03125*Power(Log(pi),2) - 1.*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi))*Log(1/s12) + 0.0078125*Power(Log(1/s12),2)) - 
              1.*s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12))*Log(s12) + 0.0078125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.03125*beta*s12*Log(beta) + beta*
               (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           2.*euler*(-0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           (-0.03125*beta*euler*s12 - 0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.015625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
      (-0.03125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.03125*beta*euler*s12 + 0.03125*beta*s12*Log(beta) - 0.03125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) + 
           0.015625*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.015625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.03125*beta*s12*Power(Log(beta),2) - 
           0.03125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) - 
           1.*beta*(s12*(0.13512741016449415 - 0.12996509635498975*Log(dS) + 0.03125*Power(Log(dS),2) + 2.*(0.064982548177494873 - 0.03125*Log(dS))*Log(pi) + 
                 0.03125*Power(Log(pi),2) - 1.*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi))*Log(1/s12) + 0.0078125*Power(Log(1/s12),2)) - 
              1.*s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12))*Log(s12) + 0.0078125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.03125*beta*s12*Log(beta) + beta*
               (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           2.*euler*(-0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           (-0.03125*beta*euler*s12 - 0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.015625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 
        2.*Log(Sin(th))*(0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*Power(1. + beta*Cos(th),2)*
      (0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.77777777777777778*Alfas*Alfas2*(1. - 1.*(beta*beta))*(-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.037037037037037037*Alfas*Alfas2*Power(1. - 1.*beta*Cos(th),2)*
      (0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 2.*Log(Sin(th))*
         (0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 
           0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. + beta*Cos(th))*(-0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(-0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(0.25*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 
        2.*Log(Sin(th))*(0.25*beta*euler*s12 + 0.25*beta*s12*Log(beta) - 0.25*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) + 
           0.125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.25*beta*s12*Power(Log(beta),2) - 0.25*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) - 
           1.*beta*(s12*(1.0810192813159532 - 1.039720770839918*Log(dS) + 0.25*Power(Log(dS),2) + 2.*(0.51986038541995898 - 0.25*Log(dS))*Log(pi) + 
                 0.25*Power(Log(pi),2) - 1.*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi))*Log(1/s12) + 0.0625*Power(Log(1/s12),2)) - 
              1.*s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12))*Log(s12) + 0.0625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.25*beta*s12*Log(beta) + beta*
               (s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           2.*euler*(-0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           (-0.25*beta*euler*s12 - 0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(-0.25*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.25*beta*euler*s12 + 0.25*beta*s12*Log(beta) - 0.25*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) + 
           0.125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.25*beta*s12*Power(Log(beta),2) - 0.25*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) - 
           1.*beta*(s12*(1.0810192813159532 - 1.039720770839918*Log(dS) + 0.25*Power(Log(dS),2) + 2.*(0.51986038541995898 - 0.25*Log(dS))*Log(pi) + 
                 0.25*Power(Log(pi),2) - 1.*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi))*Log(1/s12) + 0.0625*Power(Log(1/s12),2)) - 
              1.*s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12))*Log(s12) + 0.0625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.25*beta*s12*Log(beta) + beta*
               (s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           2.*euler*(-0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           (-0.25*beta*euler*s12 - 0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-0.25*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.25*beta*euler*s12 + 0.25*beta*s12*Log(beta) - 0.25*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) + 
           0.125*beta*s12*Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.25*beta*s12*Power(Log(beta),2) - 0.25*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) - 
           1.*beta*(s12*(1.0810192813159532 - 1.039720770839918*Log(dS) + 0.25*Power(Log(dS),2) + 2.*(0.51986038541995898 - 0.25*Log(dS))*Log(pi) + 
                 0.25*Power(Log(pi),2) - 1.*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi))*Log(1/s12) + 0.0625*Power(Log(1/s12),2)) - 
              1.*s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12))*Log(s12) + 0.0625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.25*beta*s12*Log(beta) + beta*
               (s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           2.*euler*(-0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           (-0.25*beta*euler*s12 - 0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)))*
            Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
      (0.03125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 2.*Log(Sin(th))*
         (0.03125*beta*euler*s12 + 0.03125*beta*s12*Log(beta) - 0.03125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) + 
           0.015625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.015625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.03125*beta*s12*Power(Log(beta),2) - 
           0.03125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) - 
           1.*beta*(s12*(0.13512741016449415 - 0.12996509635498975*Log(dS) + 0.03125*Power(Log(dS),2) + 2.*(0.064982548177494873 - 0.03125*Log(dS))*Log(pi) + 
                 0.03125*Power(Log(pi),2) - 1.*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi))*Log(1/s12) + 0.0078125*Power(Log(1/s12),2)) - 
              1.*s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12))*Log(s12) + 0.0078125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.03125*beta*s12*Log(beta) + beta*
               (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           2.*euler*(-0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           (-0.03125*beta*euler*s12 - 0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.015625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (0.03125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 2.*Log(Sin(th))*
         (0.03125*beta*euler*s12 + 0.03125*beta*s12*Log(beta) - 0.03125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) + 
           0.015625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.015625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.03125*beta*s12*Power(Log(beta),2) - 
           0.03125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) - 
           1.*beta*(s12*(0.13512741016449415 - 0.12996509635498975*Log(dS) + 0.03125*Power(Log(dS),2) + 2.*(0.064982548177494873 - 0.03125*Log(dS))*Log(pi) + 
                 0.03125*Power(Log(pi),2) - 1.*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi))*Log(1/s12) + 0.0078125*Power(Log(1/s12),2)) - 
              1.*s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12))*Log(s12) + 0.0078125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.03125*beta*s12*Log(beta) + beta*
               (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           2.*euler*(-0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           (-0.03125*beta*euler*s12 - 0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.015625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
      (-0.03125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.03125*beta*euler*s12 + 0.03125*beta*s12*Log(beta) - 0.03125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) + 
           0.015625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.015625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.03125*beta*s12*Power(Log(beta),2) - 
           0.03125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) - 
           1.*beta*(s12*(0.13512741016449415 - 0.12996509635498975*Log(dS) + 0.03125*Power(Log(dS),2) + 2.*(0.064982548177494873 - 0.03125*Log(dS))*Log(pi) + 
                 0.03125*Power(Log(pi),2) - 1.*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi))*Log(1/s12) + 0.0078125*Power(Log(1/s12),2)) - 
              1.*s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12))*Log(s12) + 0.0078125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.03125*beta*s12*Log(beta) + beta*
               (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           2.*euler*(-0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           (-0.03125*beta*euler*s12 - 0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.015625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-0.03125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.03125*beta*euler*s12 + 0.03125*beta*s12*Log(beta) - 0.03125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) + 
           0.015625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.015625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.03125*beta*s12*Power(Log(beta),2) - 
           0.03125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)) - 
           1.*beta*(s12*(0.13512741016449415 - 0.12996509635498975*Log(dS) + 0.03125*Power(Log(dS),2) + 2.*(0.064982548177494873 - 0.03125*Log(dS))*Log(pi) + 
                 0.03125*Power(Log(pi),2) - 1.*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi))*Log(1/s12) + 0.0078125*Power(Log(1/s12),2)) - 
              1.*s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12))*Log(s12) + 0.0078125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.03125*beta*s12*Log(beta) + beta*
               (s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           2.*euler*(-0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12))) + 
           (-0.03125*beta*euler*s12 - 0.03125*beta*s12*Log(beta) + 0.03125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.064982548177494873 - 0.03125*Log(dS) + 0.03125*Log(pi) - 0.015625*Log(1/s12)) - 0.015625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.015625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 
        2.*Log(Sin(th))*(0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.25925925925925926*Alfas*Alfas2*Power(1. - 1.*beta*Cos(th),3)*
      (-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*Power(1. + beta*Cos(th),2)*(-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*Power(1. + beta*Cos(th),2)*(-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-0.0625*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.0625*beta*euler*s12 + 0.0625*beta*s12*Log(beta) - 0.0625*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) + 
           0.03125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.03125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.0625*beta*s12*Power(Log(beta),2) - 
           0.0625*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)) - 
           1.*beta*(s12*(0.2702548203289883 - 0.25993019270997949*Log(dS) + 0.0625*Power(Log(dS),2) + 2.*(0.12996509635498975 - 0.0625*Log(dS))*Log(pi) + 
                 0.0625*Power(Log(pi),2) - 1.*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi))*Log(1/s12) + 0.015625*Power(Log(1/s12),2)) - 
              1.*s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12))*Log(s12) + 0.015625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.0625*beta*s12*Log(beta) + beta*
               (s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           2.*euler*(-0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12))) + 
           (-0.0625*beta*euler*s12 - 0.0625*beta*s12*Log(beta) + 0.0625*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.12996509635498975 - 0.0625*Log(dS) + 0.0625*Log(pi) - 0.03125*Log(1/s12)) - 0.03125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.03125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*(beta*beta))*(0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) + 
        2.*Log(Sin(th))*(0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) - 
        1.*(-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 
           0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 2.*beta*Log(beta)*
            (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(-0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(-0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(-0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-0.125*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.125*beta*euler*s12 + 0.125*beta*s12*Log(beta) - 0.125*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) + 
           0.0625*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.0625*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.125*beta*s12*Power(Log(beta),2) - 0.125*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)) - 
           1.*beta*(s12*(0.5405096406579766 - 0.51986038541995898*Log(dS) + 0.125*Power(Log(dS),2) + 2.*(0.25993019270997949 - 0.125*Log(dS))*Log(pi) + 
                 0.125*Power(Log(pi),2) - 1.*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi))*Log(1/s12) + 0.03125*Power(Log(1/s12),2)) - 
              1.*s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12))*Log(s12) + 0.03125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.125*beta*s12*Log(beta) + beta*
               (s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           2.*euler*(-0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12))) + 
           (-0.125*beta*euler*s12 - 0.125*beta*s12*Log(beta) + 0.125*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.25993019270997949 - 0.125*Log(dS) + 0.125*Log(pi) - 0.0625*Log(1/s12)) - 0.0625*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.0625*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(-0.25*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.25*beta*euler*s12 + 0.25*beta*s12*Log(beta) - 0.25*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) + 
           0.125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.25*beta*s12*Power(Log(beta),2) - 0.25*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) - 
           1.*beta*(s12*(1.0810192813159532 - 1.039720770839918*Log(dS) + 0.25*Power(Log(dS),2) + 2.*(0.51986038541995898 - 0.25*Log(dS))*Log(pi) + 
                 0.25*Power(Log(pi),2) - 1.*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi))*Log(1/s12) + 0.0625*Power(Log(1/s12),2)) - 
              1.*s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12))*Log(s12) + 0.0625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.25*beta*s12*Log(beta) + beta*
               (s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           2.*euler*(-0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           (-0.25*beta*euler*s12 - 0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.25925925925925926*Alfas*Alfas2*Power(1. - 1.*beta*Cos(th),2)*
      (-0.25*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.25*beta*euler*s12 + 0.25*beta*s12*Log(beta) - 0.25*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) + 
           0.125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.25*beta*s12*Power(Log(beta),2) - 0.25*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) - 
           1.*beta*(s12*(1.0810192813159532 - 1.039720770839918*Log(dS) + 0.25*Power(Log(dS),2) + 2.*(0.51986038541995898 - 0.25*Log(dS))*Log(pi) + 
                 0.25*Power(Log(pi),2) - 1.*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi))*Log(1/s12) + 0.0625*Power(Log(1/s12),2)) - 
              1.*s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12))*Log(s12) + 0.0625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.25*beta*s12*Log(beta) + beta*
               (s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           2.*euler*(-0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           (-0.25*beta*euler*s12 - 0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.037037037037037037*Alfas*Alfas2*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-0.25*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 2.*Log(Sin(th))*
         (0.25*beta*euler*s12 + 0.25*beta*s12*Log(beta) - 0.25*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) + 
           0.125*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.125*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.25*beta*s12*Power(Log(beta),2) - 0.25*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)) - 
           1.*beta*(s12*(1.0810192813159532 - 1.039720770839918*Log(dS) + 0.25*Power(Log(dS),2) + 2.*(0.51986038541995898 - 0.25*Log(dS))*Log(pi) + 
                 0.25*Power(Log(pi),2) - 1.*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi))*Log(1/s12) + 0.0625*Power(Log(1/s12),2)) - 
              1.*s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12))*Log(s12) + 0.0625*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.25*beta*s12*Log(beta) + beta*
               (s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           2.*euler*(-0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12))) + 
           (-0.25*beta*euler*s12 - 0.25*beta*s12*Log(beta) + 0.25*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(0.51986038541995898 - 0.25*Log(dS) + 0.25*Log(pi) - 0.125*Log(1/s12)) - 0.125*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.125*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.25925925925925926*Alfas*Alfas2*(-0.5*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.5*beta*euler*s12 + 0.5*beta*s12*Log(beta) - 0.5*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12)) + 
           0.25*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.25*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.5*beta*s12*Power(Log(beta),2) - 0.5*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12)) - 
           1.*beta*(s12*(2.1620385626319064 - 2.0794415416798359*Log(dS) + 0.5*Power(Log(dS),2) + 2.*(1.039720770839918 - 0.5*Log(dS))*Log(pi) + 
                 0.5*Power(Log(pi),2) - 1.*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi))*Log(1/s12) + 0.125*Power(Log(1/s12),2)) - 
              1.*s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12))*Log(s12) + 0.125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.5*beta*s12*Log(beta) + beta*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12))) + 
           2.*euler*(-0.5*beta*s12*Log(beta) + 0.5*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12))) + 
           (-0.5*beta*euler*s12 - 0.5*beta*s12*Log(beta) + 0.5*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.25*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.25925925925925926*Alfas*Alfas2*(1. - 1.*(beta*beta))*(-0.5*beta*s12*Power(Log(Sin(th)),2)*Sin(th) - 
        2.*Log(Sin(th))*(0.5*beta*euler*s12 + 0.5*beta*s12*Log(beta) - 0.5*beta*s12*Log((muR*muR)/s12) - 
           1.*beta*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12)) + 
           0.25*beta*s12*Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))))*Sin(th) + 
        (-0.25*beta*(2.*(euler*euler) - 0.33333333333333333*(pi*pi))*s12 - 0.5*beta*s12*Power(Log(beta),2) - 0.5*beta*s12*Power(Log((muR*muR)/s12),2) + 
           2.*beta*Log(beta)*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12)) - 
           1.*beta*(s12*(2.1620385626319064 - 2.0794415416798359*Log(dS) + 0.5*Power(Log(dS),2) + 2.*(1.039720770839918 - 0.5*Log(dS))*Log(pi) + 
                 0.5*Power(Log(pi),2) - 1.*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi))*Log(1/s12) + 0.125*Power(Log(1/s12),2)) - 
              1.*s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12))*Log(s12) + 0.125*s12*Power(Log(s12),2)) - 
           2.*Log((muR*muR)/s12)*(-0.5*beta*s12*Log(beta) + beta*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12))) + 
           2.*euler*(-0.5*beta*s12*Log(beta) + 0.5*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12))) + 
           (-0.5*beta*euler*s12 - 0.5*beta*s12*Log(beta) + 0.5*beta*s12*Log((muR*muR)/s12) + 
              beta*(s12*(1.039720770839918 - 0.5*Log(dS) + 0.5*Log(pi) - 0.25*Log(1/s12)) - 0.25*s12*Log(s12)))*
            Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
           0.25*beta*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
              0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
                  (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
              2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
              2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
                 (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th)))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)));
    
//    std::complex<double> sigma = (-42.66666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
//   (42.66666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
//   (21.33333333333333*Alfas*Alfas2*beta*MGl2*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
//   (21.33333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
//   (10.66666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
//   (21.33333333333333*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
//   (10.66666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
//   (2.37037037037037*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) - 
//   (0.5925925925925926*Alfas*Alfas2*beta*(-1. + beta*beta)*(s12*s12)*(-1. + beta*Cos(th))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
//   (1.185185185185185*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
//   (0.5925925925925926*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.5925925925925926*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (1.185185185185185*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (42.66666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
//   (42.66666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
//   (21.33333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
//   (10.66666666666667*Alfas*Alfas2*beta*(s12*s12)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) - 
//   (10.66666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
//   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.07407407407407407*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.5925925925925926*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.5925925925925926*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.07407407407407407*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.1481481481481481*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.03703703703703704*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.07407407407407407*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.07407407407407407*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.07407407407407407*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.03703703703703704*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.07407407407407407*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.07407407407407407*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (1.185185185185185*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (1.185185185185185*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.5925925925925926*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.07407407407407407*Alfas*Alfas2*beta*Power(-1. + beta*beta,4)*Power(s12,3)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)
//     - (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)
//     - (9.481481481481481*Alfas*Alfas2*(-1. + beta*beta)*(MGl2*MGl2)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*
//      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (9.481481481481481*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*
//      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (4.740740740740741*Alfas*Alfas2*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*
//      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (0.5925925925925926*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*
//      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (1.185185185185185*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (1.185185185185185*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (1.185185185185185*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (9.481481481481481*Alfas*Alfas2*(MGl2*MGl2)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (9.481481481481481*Alfas*Alfas2*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (4.740740740740741*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.5925925925925926*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (1.185185185185185*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (1.185185185185185*Alfas*Alfas2*Power(s12,3)*Power(-1. + beta*Cos(th),3)*(1. + beta*Cos(th))*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (1.185185185185185*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),2)*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (1.185185185185185*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (1.185185185185185*Alfas*Alfas2*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),3)*
//      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.07407407407407407*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
//      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.07407407407407407*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
//      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
//      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.07407407407407407*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.07407407407407407*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
//    + (0.01851851851851852*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
//    + (0.01851851851851852*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
//    - (0.01851851851851852*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
//    - (0.03703703703703704*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
//    - (0.01851851851851852*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
//      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
//      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
//    + (0.03703703703703704*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
//    - (0.03703703703703704*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
//    - (1.185185185185185*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (1.185185185185185*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (0.5925925925925926*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (0.07407407407407407*Alfas*Alfas2*beta*Power(-1. + beta*beta,4)*Power(s12,3)*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (1.185185185185185*Alfas*Alfas2*beta*(-1. + beta*beta)*(MGl2*MGl2)*s12*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (1.185185185185185*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.5925925925925926*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (0.07407407407407407*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
//   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),3)*
//      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
//    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
//   (0.07407407407407407*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
//        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.07407407407407407*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
//        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.03703703703703704*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.01851851851851852*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.01851851851851852*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.07407407407407407*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
//        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*(-1. + beta*beta)*s12*(1. + beta*Cos(th))*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.01851851851851852*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.01851851851851852*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.01851851851851852*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
//   (0.01851851851851852*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
//   (0.07407407407407407*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
//        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.07407407407407407*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
//        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.03703703703703704*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.01851851851851852*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.01851851851851852*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.07407407407407407*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
//        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
//    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.03703703703703704*Alfas*Alfas2*(-1. + beta*beta)*s12*(1. + beta*Cos(th))*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.01851851851851852*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.01851851851851852*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.01851851851851852*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
//   (0.01851851851851852*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
//      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
//        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
//   (0.002314814814814815*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),2)*
//      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
//        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
//        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
//        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
//      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
//        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
//        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
//        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
//   (0.001157407407407407*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
//      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
//        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
//        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
//        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
//   (0.005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.00462962962962963*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.002314814814814815*Alfas*Alfas2*beta*(1. + beta*Cos(th))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
//        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
//        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
//   (0.04050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
//        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
//        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
//   (0.004050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
//   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
//   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) - 
//   (0.03240740740740741*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
//   (0.004050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
//   (0.003472222222222222*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
//      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
//        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
//        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
//        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
//   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
//      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
//        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
//        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
//        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
//   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
//      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
//        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
//        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
//        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
//   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
//      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
//        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
//        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
//        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
//   (0.02430555555555556*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
//        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
//        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.0162037037037037*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
//        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
//        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
//   (0.0162037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
//   (0.03240740740740741*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
//   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
//   (0.004050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
//   (0.008101851851851852*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),3)*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
//   (0.004050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
//   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.008101851851851852*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
//   (0.06481481481481481*Alfas*Alfas2*beta*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
//        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
//        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.0162037037037037*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
//        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
//        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.06481481481481481*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
//        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
//        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
//   (0.0162037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
//   (0.03240740740740741*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),2)*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
//   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.008101851851851852*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
//   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
//   (0.00462962962962963*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
//      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
//        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
//        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
//         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
//        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
//           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
//               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
//           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
//           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
//              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
//    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2);
    
    //std::cout << sigma.imag() << '\n';
    
    
    
    ff[0] = to_fb * sigma.real() 
            * pdf_nlo->xfxQ(2, x1, 1500.)/x1 
            * pdf_nlo->xfxQ(2, x2, 1500.)/x2;
    
    // jakobian
    ff[0] *= pi*Power(-4*m_sqr + S,2)*xx[0] / 
            (S*(-4*m_sqr*(-1 + xx[0]) + S*xx[0]));
    return 0;
}

// @todo if muR != muF one needs one more term here

int XSection_SC::integrand_c1(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

    double m_sqr = pow( 1500., 2);
    double MGl2 = pow ( gluino_mass, 2 );
    double muF = 1500.;
    double muR = 1500.;
    
    // integration variables
    double x1 = 4. * m_sqr/S + (1 - 4. * m_sqr/S ) * xx[0];
	double x2 = 4. * m_sqr /(S * x1) + (1 - 4. * m_sqr/(S * x1)) * xx[1];
    double th = pi * xx[2];
    
    double s12 = x1 * x2 * S;
    double beta = sqrt( 1 - 4 * m_sqr/s12 );
    
    double Alfas = pdf_nlo->alphasQ(1500.);
    double Alfas2 = pow( Alfas, 2);
    
    ff[0] = 2 * to_fb
            * pdf_nlo->xfxQ(2, x1, 1500.)/x1 * pdf_nlo->xfxQ(2, x2, 1500.)/x2
            * ((-4*Alfas*Alfas2*Power(beta,3)*s12*(Power(beta,4)*(s12*s12) + 4*(beta*beta)*s12*(2*MGl2 + s12) + Power(4*MGl2 + s12,2) + 2*(beta*beta)*(s12*s12)*cos(2*th))*
     (-Log(Power(muR,4)/(s12*(muF*muF))) + 2*(euler + Log((beta*Sin(th))/(8.*pi))))*Power(Sin(th),3))/
   (9.*Power(Power(4*MGl2 + s12 + beta*beta*s12,2) - 4*(beta*beta)*(s12*s12)*Power(cos(th),2),2)))
            * 4./3. * (2 * log(dS) + 3./2.);
    
    ff[0] *= (pi*Power(-4*m_sqr + S,2)*xx[0])/(S*(-4*m_sqr*(-1 + xx[0]) + S*xx[0]));
    
    return 0;
}

int XSection_SC::integrand_c2(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
    
    double m_sqr = pow( 1500., 2);
    double MGl2 = pow ( gluino_mass, 2 );
    double muF = 1500.;
    
    // integration variables
    double x1 = 4. * m_sqr/S + (1 - 4. * m_sqr/S ) * xx[0];
	double x2 = 4. * m_sqr /(S * x1) + (1 - 4. * m_sqr/(S * x1)) * xx[1];
    double y = x1 + ( 1 - dS - x1 ) * xx[2];
    
    double s12 = x1 * x2 * S;
    double beta = sqrt( 1 - 4 * m_sqr/s12 );
    
    double Alfas = pdf_nlo->alphasQ(1500.);
    double Alfas2 = pow( Alfas, 2);
    
    ff[0] = 2 * to_fb 
            * ((4*(1 - y))/3. + (4*(1 + Power(y,2))*Log((dC*s12*(1 - y))/(Power(muF,2)*y)))/(3.*(1 - y))) 
            * 1./y * pdf_nlo->xfxQ(2, std::min(x1/y, 1.), 1500.)/std::min(x1/y, 1.) 
            * pdf_nlo->xfxQ(2, x2, 1500.)/x2 
            * Alfas/two_pi 
            // 2->2 cross-section
            * (Alfas2*pi*(-8*beta*s12 + 2*(4*MGl2 + s12 + Power(beta,2)*s12)*Log((4*MGl2 + Power(1 + beta,2)*s12)/
            (4*MGl2 + Power(-1 + beta,2)*s12))))/(9.*Power(s12,2));
    
    // multiply by jakobian of integration variable transformation
    ff[0] *= (Power(-4*m_sqr + S,2)*xx[0]*(4*m_sqr*(-1 + xx[0]) - S*(-1 + dS + xx[0])))/
        (Power(S,2)*(-4*m_sqr*(-1 + xx[0]) + S*xx[0]));
    
  return 0;
}
