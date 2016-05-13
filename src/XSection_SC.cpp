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
  constexpr long long int neval_max { 1000000000 }; 
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

  cubareal integral_c[ncomp], error_c[ncomp], prob_c[ncomp];
  llVegas( ndim, ncomp, integrand_c, NULL, 1,
           accuracy_rel_c, accuracy_abs, 8 | 1, 0,
           neval_min, neval_max, nstart, nincrease, nbatch,
           gridno, state_file, NULL,
           &neval, &fail, integral_c, error_c, prob_c );
  
  std::array <double, 3> result_finite { 
      integral_sc[0] + integral_c[0], 
      sqrt( pow( error_sc[0], 2) + pow ( error_c[0], 2) ),
      prob_sc[0] + prob_c[0] 
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
    
    std::complex<double> sigma = (-42.66666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (42.66666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.33333333333333*Alfas*Alfas2*beta*MGl2*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.33333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (10.66666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) + 
   (21.33333333333333*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
   (10.66666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),4) - 
   (2.37037037037037*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) - 
   (0.5925925925925926*Alfas*Alfas2*beta*(-1. + beta*beta)*(s12*s12)*(-1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
   (1.185185185185185*Alfas*Alfas2*beta*(s12*s12)*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),3) + 
   (0.5925925925925926*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.5925925925925926*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (1.185185185185185*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (42.66666666666667*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (42.66666666666667*Alfas*Alfas2*beta*MGl2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (21.33333333333333*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (10.66666666666667*Alfas*Alfas2*beta*(s12*s12)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) - 
   (10.66666666666667*Alfas*Alfas2*beta*(s12*s12)*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),4) + 
   (0.2962962962962963*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.07407407407407407*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.5925925925925926*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.5925925925925926*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),2)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.07407407407407407*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.1481481481481481*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.1481481481481481*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.03703703703703704*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.07407407407407407*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.07407407407407407*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.07407407407407407*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.03703703703703704*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.1481481481481481*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.07407407407407407*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.07407407407407407*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (1.185185185185185*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.185185185185185*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.5925925925925926*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.07407407407407407*Alfas*Alfas2*beta*Power(-1. + beta*beta,4)*Power(s12,3)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(1. + beta*Cos(th),2)*(2.*Log(dS) - 1.*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)
     - (0.2962962962962963*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.2962962962962963*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)
     - (9.481481481481481*Alfas*Alfas2*(-1. + beta*beta)*(MGl2*MGl2)*s12*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (9.481481481481481*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (4.740740740740741*Alfas*Alfas2*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.5925925925925926*Alfas*Alfas2*Power(-1. + beta*beta,3)*Power(s12,3)*(-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*
      Sin(th))/(Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.185185185185185*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.185185185185185*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.185185185185185*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (9.481481481481481*Alfas*Alfas2*(MGl2*MGl2)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (9.481481481481481*Alfas*Alfas2*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (4.740740740740741*Alfas*Alfas2*(-1. + beta*beta)*MGl2*(s12*s12)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.5925925925925926*Alfas*Alfas2*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.185185185185185*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.185185185185185*Alfas*Alfas2*Power(s12,3)*Power(-1. + beta*Cos(th),3)*(1. + beta*Cos(th))*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.185185185185185*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.185185185185185*Alfas*Alfas2*(-1. + beta*beta)*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.185185185185185*Alfas*Alfas2*Power(s12,3)*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),3)*
      (-1.*Log((1. + beta)/(1. - 1.*beta)) + 2.*beta*Log(dS) - 1.*beta*Log((muR*muR)/s12))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.07407407407407407*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.07407407407407407*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),3)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.07407407407407407*Alfas*Alfas2*beta*s12*(1. + beta*Cos(th))*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.07407407407407407*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    + (0.01851851851851852*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    + (0.01851851851851852*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.01851851851851852*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.03703703703703704*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)\
    - (0.01851851851851852*Alfas*Alfas2*beta*s12*Power(1. + beta*Cos(th),3)*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*s12*(4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*
      Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.03703703703703704*Alfas*Alfas2*beta*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    + (0.03703703703703704*Alfas*Alfas2*beta*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    - (0.03703703703703704*Alfas*Alfas2*beta*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)\
    - (1.185185185185185*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*(MGl2*MGl2)*s12*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (1.185185185185185*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.5925925925925926*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*MGl2*(s12*s12)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.07407407407407407*Alfas*Alfas2*beta*Power(-1. + beta*beta,4)*Power(s12,3)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(-1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.185185185185185*Alfas*Alfas2*beta*(-1. + beta*beta)*(MGl2*MGl2)*s12*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (1.185185185185185*Alfas*Alfas2*beta*(-1. + beta*beta)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.5925925925925926*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*MGl2*(s12*s12)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.07407407407407407*Alfas*Alfas2*beta*Power(-1. + beta*beta,3)*Power(s12,3)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) - 
   (0.1481481481481481*Alfas*Alfas2*beta*Power(-1. + beta*beta,2)*Power(s12,3)*(-1. + beta*Cos(th))*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*Power(s12,3)*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.1481481481481481*Alfas*Alfas2*beta*(-1. + beta*beta)*Power(s12,3)*Power(1. + beta*Cos(th),3)*
      (4.*Power(Log(dS),2) - 4.*Log(dS)*Log((muR*muR)/s12) + Power(Log((muR*muR)/s12),2))*Sin(th))/
    (Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2)*Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2)) + 
   (0.07407407407407407*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.07407407407407407*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.03703703703703704*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.01851851851851852*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.01851851851851852*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.07407407407407407*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*(-1. + beta*beta)*s12*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.01851851851851852*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.01851851851851852*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.01851851851851852*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) - 
   (0.01851851851851852*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 - 2.*beta*s12*Cos(th),2) + 
   (0.07407407407407407*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*Power(-1. + beta*beta,2)*s12*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.07407407407407407*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.03703703703703704*Alfas*Alfas2*(-1. + beta*beta)*s12*(-1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.01851851851851852*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(-1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.01851851851851852*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.07407407407407407*Alfas*Alfas2*s12*(1. + beta*Cos(th))*(-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 
        2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/
    Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.03703703703703704*Alfas*Alfas2*(-1. + beta*beta)*s12*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.01851851851851852*Alfas*Alfas2*s12*Power(-1. + beta*Cos(th),2)*(1. + beta*Cos(th))*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.01851851851851852*Alfas*Alfas2*(-1. + beta*beta)*s12*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.01851851851851852*Alfas*Alfas2*s12*(-1. + beta*Cos(th))*Power(1. + beta*Cos(th),2)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) - 
   (0.01851851851851852*Alfas*Alfas2*s12*Power(1. + beta*Cos(th),3)*
      (-1.*Power(Log((1. + beta)/(1. - 1.*beta)),2) + 4.*Log((1. + beta)/(1. - 1.*beta))*Log(dS) - 2.*Log((1. + beta)/(1. - 1.*beta))*Log((muR*muR)/s12) - 
        4.*PolyLog(2.,(2.*beta)/(1. + beta)))*Sin(th))/Power(4.*MGl2 + s12 + beta*beta*s12 + 2.*beta*s12*Cos(th),2) + 
   (0.002314814814814815*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.001157407407407407*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.00462962962962963*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.002314814814814815*Alfas*Alfas2*beta*(1. + beta*Cos(th))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.04050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.004050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) - 
   (0.03240740740740741*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) - 
   (0.004050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. - 1.*beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(-1.*beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. - 1.*beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(-1.*beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. - 1.*beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.003472222222222222*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.0005787037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) + 0.5*s12*Power(Log((muR*muR)/s12),2) - 
        1.*s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) + 0.5*s12*Power(Log(s12),2) + Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) - 
        1.*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) - 
        1.*s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2) + 
   (0.02430555555555556*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0162037037037037*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.0162037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.03240740740740741*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.004050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*Power(1. - 1.*beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.008101851851851852*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),3)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.004050925925925926*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.008101851851851852*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. - 1.*beta*Cos(th)),2) - 
   (0.06481481481481481*Alfas*Alfas2*beta*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.0162037037037037*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.06481481481481481*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 
        0.5*s12*Power(Log((muR*muR)/s12),2) + s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 
        1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.0162037037037037*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) + 
   (0.03240740740740741*Alfas*Alfas2*beta*Power(1. - 1.*beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2)) - 
   (0.008101851851851852*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.008101851851851852*Alfas*Alfas2*beta*Power(1. + beta*Cos(th),2)*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(MGl2 - 0.25*(1. - 1.*(beta*beta))*s12 + 0.5*s12*(1. - 1.*beta*Cos(th)),2) + 
   (0.001157407407407407*Alfas*Alfas2*beta*(1. - 1.*(beta*beta))*(1. - 1.*beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    ((1. + beta*Cos(th))*Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2)) + 
   (0.00462962962962963*Alfas*Alfas2*beta*(1. - 1.*beta*Cos(th))*(1. + beta*Cos(th))*
      (-1.*s12*(2.*Power(Log(dS),2) + 2.*Log(dS)*Log(1/s12) + 0.5*Power(Log(1/s12),2)) - 0.5*s12*Power(Log((muR*muR)/s12),2) + 
        s12*(-2.*Log(dS) - 1.*Log(1/s12))*Log(s12) - 0.5*s12*Power(Log(s12),2) - 1.*Log((muR*muR)/s12)*(s12*(-2.*Log(dS) - 1.*Log(1/s12)) - 1.*s12*Log(s12)) + 
        (s12*(-2.*Log(dS) - 1.*Log(1/s12)) + s12*Log((muR*muR)/s12) - 1.*s12*Log(s12))*
         Log(Power(1. + beta*Cos(th),2)/(1. - 1.*(beta*beta)*Power(Cos(th),2) - 1.*(beta*beta)*Power(Sin(th),2))) + 
        s12*(-1.*Power(Log((1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))),2) + 
           0.5*Power(Log((1. + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/
               (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))),2) + 
           2.*PolyLog(2.,(beta*Cos(th) - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2)))/(1. + beta*Cos(th))) - 
           2.*PolyLog(2.,(-1.*(beta*Cos(th) + Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))/
              (1. - 1.*Sqrt(beta*beta*Power(Cos(th),2) + beta*beta*Power(Sin(th),2))))))*Sin(th))/
    Power(-1.*MGl2 + 0.25*(1. - 1.*(beta*beta))*s12 - 0.5*s12*(1. + beta*Cos(th)),2);
    
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

int XSection_SC::integrand_c(const int *ndim, const cubareal xx[],
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
