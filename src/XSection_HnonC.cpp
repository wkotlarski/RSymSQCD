#include "XSection_HnonC.hpp"
#include <cassert>
#include <iostream>

// @todo this is absolutely neeed but I don't know why
using namespace std;

std::array<double, 3> XSection_HnonC::integrate() {

   //  integral dimension, number of integrands
   constexpr int ndim { 7 }, ncomp { 1 };
   //  accuraccy
   double accuracy_rel { prec_hnc }, 
      accuracy_abs { 1e-12 };

   constexpr int neval_min = 10'000;
   long long int neval;
   constexpr long long int neval_max { 1'000'000'000'000 }; 

   // technical (Vegas specific) stuff
   constexpr int nstart = 200'000;
   constexpr int nincrease = 10000;
   constexpr int nbatch = 1000;
   constexpr int gridno = 0;
   const char* state_file = "";
   int nregions, fail;

   cubareal integral[ncomp], error[ncomp], prob[ncomp];
   llVegas( ndim, ncomp, integrand, NULL, 1,
      accuracy_rel, accuracy_abs, 1, 0,
      neval_min, neval_max, nstart, nincrease, nbatch,
      gridno, state_file, NULL,
      &neval, &fail, integral, error, prob );
   
   std::array <double, 3> result_finite
      { integral[0], error[0], prob[0] };
   return result_finite;
}

int XSection_HnonC::integrand(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

   double m_sqr = m1 * m1;
  
   /*
   * 3-body phase space parametrization based on
   * http://www.t39.ph.tum.de/T39_files/T39_people_files/duell_files/Dipl-MultiPion.pdf
   */

   // failsafe (this should never happen)
   // but sometimes does for suave
   if (
        xx[0] < 0 || xx[0] > 1        // gluon energy
        || xx[1] < 0 || xx[1] > 1     // sgluon energy
        || xx[2] < 0 || xx[2] > 1     // angle
        || xx[3] < 0 || xx[3] > 1     // angle
        || xx[4] < 0 || xx[4] > 1     // angle
        || xx[5] < 0 || xx[5] > 1   // Bjorken x
        || xx[6] < 0 || xx[6] > 1   // Bjorken x
        ) {
      ff[0] = 0;
      return 0;
   }

	double x1 = 4. * m_sqr/S + (1. - 4. * m_sqr/S ) * xx[5];
	double x2 = 4. * m_sqr /(S * x1) + (1. - 4. * m_sqr/(S * x1)) * xx[6];
	double shat = x1 * x2 * S;
	double shat_sqrt = sqrt( shat );

	double Ej_max = shat_sqrt/2. - 2. * m_sqr/shat_sqrt;
   
   // dS < 0 signs that no soft cut should be applied
	if ( Ej_max < dS * shat_sqrt/2. ) {
	  ff[0] = 0;
	  return 0;
	}
   
	double Ej = dS * shat_sqrt/2. + ( Ej_max - dS * shat_sqrt/2.) * xx[0];  
   
	double c = shat - 2. * shat_sqrt * Ej;
	// Eq. 4.5
	double E1_max = ( shat_sqrt - Ej) * c + Ej * sqrt( (c - 2. * m_sqr) 
    * (c - 2 * m_sqr) - 4. * m_sqr * m_sqr );
	E1_max /= 2 * c;
   double E1_min = ( shat_sqrt - Ej) * c - Ej * sqrt( (c - 2 * m_sqr) 
      * (c - 2. * m_sqr) - 4. * m_sqr * m_sqr );
   E1_min /= 2. * c;
	double E1 = E1_min + (E1_max - E1_min) * xx[1];

	// Eq. 4.2 with E2 = Ej
   double cosx = (shat - 2 * shat_sqrt * ( E1 + Ej ) + 2. * Ej * E1 )/
      (2. * Ej * sqrt( ( E1 - m1 ) * ( E1 + m1 ) ) );

   // check if due to numerics abs(cosx) is not > 1
   if ( cosx > 1 || cosx < -1)  {
      std::cout << "Warning! cos(x) = " << cosx << " out of bounds. Setting it +/- 1." <<  std::endl;
      cosx = cosx > 0 ? 1.0 : -1.0;
   }

   // initialize matrix of particles momenta
   std::vector< double* > p;

   // incoming parton 1 momenta
   p.push_back(new double[4]);
   p[0][0] = shat_sqrt/2;
   p[0][1] = 0;
   p[0][2] = 0;
   p[0][3] = shat_sqrt/2;

   // incoming parton 2 momenta
   p.push_back(new double[4]);
   p[1][0] = shat_sqrt/2;
   p[1][1] = 0;
   p[1][2] = 0;
   p[1][3] = - shat_sqrt/2;

   // 1st sgluons momenta
   p.push_back(new double[4]);
   p[2][0] = E1;
   double p1 = sqrt( ( E1 - m1 ) * ( E1 + m1 ) );
   p[2][1] = p1 * sin( pi * xx[2] ) * cos( two_pi * xx[3] );
   p[2][2] = p1 * sin( pi * xx[2] ) * sin( two_pi * xx[3] );
   p[2][3] = p1 * cos( pi * xx[2] );

   geom3::UnitVector3 rotation_axis(
      sin( pi * xx[2] ) * cos( two_pi * xx[3] ),
      sin( pi * xx[2] ) * sin( two_pi * xx[3] ),
      cos( pi * xx[2] )
   );

   // construct rotation by angle xx[4] * two_pi
   // around sgluon momentum
   geom3::Rotation3 rot( rotation_axis, xx[4] * two_pi);

   /*
    *  kinematics was solved for the cosx of angle
    *  between parton and SUSY particle
    *
    *  probably it would be also ok to skip it hoping that
    *  periodicity of sin and cos would solve the thing
    */
   double parton_theta, parton_phi;
   if( pi * xx[2] + acos(cosx) < pi ) {
      parton_theta = pi * xx[2] + acos(cosx);
      parton_phi = two_pi * xx[3];
   }
   else {
      parton_theta = two_pi - pi * xx[2] - acos(cosx);
      parton_phi = two_pi * xx[3] + pi;
   }

   geom3::Vector3 p_parton(
      Ej * sin( parton_theta ) * cos( parton_phi ),
      Ej * sin( parton_theta ) * sin( parton_phi ),
      Ej * cos( parton_theta )
   );

   // rotate parton momentum
   geom3::Vector3 p_temp = rot.rotate(p_parton);

   // set parton momenta
   double*  p_temp_2 = new double [4];
   p_temp_2[0] = Ej;
   p_temp_2[1] = p_temp.x();
   p_temp_2[2] = p_temp.y();
   p_temp_2[3] = p_temp.z();

   // 2nd sgluon momenta
   p.push_back(new double[4]);
   p[3][0] = shat_sqrt - E1 - Ej;
   for( int i = 1; i < 4; ++i) p[3][i] = - p[2][i] - p_temp_2[i];
  
   assert( abs(
      (pow(p[2][0], 2) - pow(p[2][1], 2) - pow(p[2][2], 2) - pow(p[2][3], 2))/(m1 * m1) - 1) < 1e-10 
         && p[2][0] >= m1
   );
   assert( abs(
      (pow(p[3][0], 2) - pow(p[3][1], 2) - pow(p[3][2], 2) - pow(p[3][3], 2))/(m2 * m2) - 1) < 1e-10 
      && p[3][0] >= m2
   );
  
   // write parton momentum to momentum matrix p
   p.push_back(new double[4]);
   for(int i = 0; i < 4; ++i) p[4][i] = p_temp_2[i];
   delete[] p_temp_2;
      
   double t15 = p[0][0] * p[4][0] - p[0][3] * p[4][3];
   t15 = - 2. * t15;
   double t25 = p[1][0] * p[4][0] - p[1][3] * p[4][3];
   t25 = - 2. * t25;

   // check if we are not in the collinear region
   // if yes, return
   if( -t15 < dC * shat_sqrt * Ej || -t25 < dC * shat_sqrt * Ej ) {
      ff[0] = 0;
      return 0;
	}
        
   double ME2 = (processID->*processID->matrixelementReal_HnonC)(p);
   assert( ME2 >= 0 );
   
   // delete (otherwise causes memory leak)
   for(std::vector<double*>::iterator i = p.begin(); i != p.end(); ++i) {
      delete (*i);
      *i = 0;
   }
   p.clear();
   p.shrink_to_fit();
  
   // some final factors
   ME2 *= to_fb;
   ME2 /=  2 * shat;

   ME2 *= sin( xx[2] * pi ) * 4;
   ME2 /= 256 * pi_sqr;

   double pdf_flux = 0.0;
   for (const auto& f : processID->flav) {
      pdf_flux += f.at(2) * pdf->xfxQ( f.at(0), x1, mu_f ) * pdf->xfxQ( f.at(1), x2, mu_f );
   }
   pdf_flux /= x1 * x2;

   ME2 *=  pdf_flux;
   double xx0 = xx[5];
   double xx1 = xx[6];
   double xx2 = xx[0];

   double jacobian =  (xx0*(-4*dS*pow(m1,2) + (-1 + dS)*xx0*xx1*(-S + 4*pow(m1,2)))*
     (xx0*xx1*xx2*(S - 4*pow(m1,2)) + 
       dS*(-1 + xx2)*(-(S*xx0*xx1) + 4*(-1 + xx0*xx1)*pow(m1,2)))*pow(S,-1)*
     pow(S - 4*pow(m1,2),2)*pow(S*xx0 - 4*(-1 + xx0)*pow(m1,2),-1)*
     pow(S*xx0*xx1 + (4 - 4*xx0*xx1)*pow(m1,2),-1)*
     pow((-1 + dS)*S*xx0*xx1*(-1 + xx2) + 
       4*(1 + xx0*xx1*(-1 + xx2) - dS*(-1 + xx0*xx1)*(-1 + xx2))*pow(m1,2),-1)*
     pow((-1 + xx2)*(S*xx0*xx1*(-1 + dS + xx2 - dS*xx2) + 
         4*(-1 + xx0*xx1 + dS*(-1 + xx0*xx1)*(-1 + xx2) - xx0*xx1*xx2)*pow(m1,2))*
       (-4*dS*pow(m1,2) + (-1 + dS)*xx0*xx1*(-S + 4*pow(m1,2))),0.5))/4.;
   
   ff[0] = ME2 * abs(jacobian);

	return 0;
}
