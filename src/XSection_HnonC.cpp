#include "XSection_HnonC.hpp"
#include <cassert>
#include <iostream>

// @todo this is absolutely neeed but I don't know why
using namespace std;
/*
extern "C" {
   void reshuffle_momenta_2__ (double[], double*, double*, double*, double*, double[]);
}
*/
std::array<double, 3> XSection_HnonC::integrate() {

   //  integral dimension, number of integrands
   constexpr int ndim { 7 }, ncomp { 1 };
   //  accuraccy
   const double accuracy_rel { pow(10., -vm["precision-hard"].as<int>()) };
   constexpr accuracy_abs { 1e-12 };

   constexpr int neval_min = 10'000;
   long long int neval;
   constexpr long long int neval_max { 1'000'000'000'000 }; 

   // technical (Vegas specific) stuff
   constexpr int nstart = 200'000;
   constexpr int nincrease = 10000;
   constexpr int nbatch = 1000;
   constexpr int gridno = 0;
   const int flags = vm["verbosity-hard"].as<int>();
   constexpr int seed = 0;
   const char* state_file = "";
   int nregions, fail;

   cubareal integral[ncomp], error[ncomp], prob[ncomp];
   llVegas( ndim, ncomp, integrand, NULL, 1,
      accuracy_rel, accuracy_abs, flags, seed,
      neval_min, neval_max, nstart, nincrease, nbatch,
      gridno, state_file, NULL,
      &neval, &fail, integral, error, prob );
   
   std::array <double, 3> result_finite
      { integral[0], error[0], prob[0] };
   return result_finite;
}

int XSection_HnonC::integrand(const int *ndim, const cubareal xx[],
   const int *ncomp, cubareal ff[], void *userdata) {

   /*
   *  3-body phase space parametrization based on
   *  http://www.t39.ph.tum.de/T39_files/T39_people_files/duell_files/Dipl-MultiPion.pdf
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

   const double m_sqr = m1 * m1;

	const double x1 = 4. * m_sqr/S + (1. - 4. * m_sqr/S ) * xx[5];
	const double x2 = 4. * m_sqr /(S * x1) + (1. - 4. * m_sqr/(S * x1)) * xx[6];
	const double shat = x1 * x2 * S;
	const double shat_sqrt = sqrt( shat );

   const double s35_min = m1*m1;
   const double s35_max = pow(shat_sqrt - m2, 2);
   const double s35 = s35_min + (s35_max - s35_min) * xx[0];
	const double E2 = -0.5 * (s35 - shat - m1*m1)/shat_sqrt;
   
	const double c = shat - 2 * shat_sqrt * E2 + m1*m1 + m2*m2;
	// Eq. 4.5
	double E1_max = ( shat_sqrt - E2) * c + sqrt(E2*E2-m2*m2) * sqrt( (c - 2. * m_sqr) 
    * (c - 2 * m_sqr) );
	E1_max /= 2 * (c - m1*m1);
   double E1_min = ( shat_sqrt - E2) * c - sqrt(E2*E2-m2*m2) * sqrt( (c - 2 * m_sqr) 
      * (c - 2. * m_sqr)  );
   E1_min /= 2. * (c - m1*m1);

   const double s45_min = shat + m2*m2 - 2 * shat_sqrt * E1_max;
   const double s45_max = shat + m2*m2 - 2 * shat_sqrt * E1_min;
   const double s45 = s45_min + (s45_max - s45_min) * xx[1];
   const double E1 = -0.5 * (s45 - shat - m2*m2)/shat_sqrt;

   if ( shat_sqrt - E1 - E2 < 0.5 * dS * shat_sqrt) { ff[0] = 0.; return 0; }

	// eq. 4.2
   const double cosx = (shat - 2 * shat_sqrt * ( E1 + E2 ) + 2 * E2 * E1 + m1*m1 + m2*m2 )/
      (2. * sqrt( ( E1 - m1 ) * ( E1 + m1 ) ) * sqrt( ( E2 - m2 ) * ( E2 + m2 ) ) );

   // check if due to numerics |cos(x)| is not > 1
   // if yes, return 0 and continue
   if ( cosx > 1 || cosx < -1)  {
      std::cout << "Warning! 1 - |cos(x)| = " << 1 - abs(cosx) << "  - Skipping the phase space point.\n";
      ff[0] = 0;
      return 0;
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

   double p2 = sqrt ((E2-m2) * (E2+m2));
   geom3::Vector3 p_parton(
      p2 * sin( parton_theta ) * cos( parton_phi ),
      p2 * sin( parton_theta ) * sin( parton_phi ),
      p2 * cos( parton_theta )
   );

   // rotate parton momentum
   geom3::Vector3 p_temp = rot.rotate(p_parton);

   // 2nd sgluon momenta
   p.push_back(new double[4]);
   p[3][0] = E2;
   p[3][1] = p_temp.x();
   p[3][2] = p_temp.y();
   p[3][3] = p_temp.z();

   // parton
   p.push_back(new double[4]);
   p[4][0] = shat_sqrt - E1 - E2;
   for( int i = 1; i < 4; ++i) p[4][i] = - p[2][i] - p[3][i];
  
   assert( abs(
      (pow(p[2][0], 2) - pow(p[2][1], 2) - pow(p[2][2], 2) - pow(p[2][3], 2))/(m1 * m1) - 1) < 1e-10 
         && p[2][0] >= m1
   );
   if( abs(
      (pow(p[3][0], 2) - pow(p[3][1], 2) - pow(p[3][2], 2) - pow(p[3][3], 2))/(m2 * m2) - 1) > 1e-10 
      || p[3][0] < m2 ) {
      std::cout << "Error in kinematics. " <<
              abs(
      (pow(p[3][0], 2) - pow(p[3][1], 2) - pow(p[3][2], 2) - pow(p[3][3], 2))/(m2 * m2) - 1) << " " << p[3][0] << '\n';
   }
  
      
   const double t15 = -2. * (p[0][0] * p[4][0] - p[0][3] * p[4][3]);
   const double t25 = -2. * (p[1][0] * p[4][0] - p[1][3] * p[4][3]);

   // check if we are not in the collinear region
   // if yes, return
   if( -t15 < dC * shat_sqrt * p[4][0] || -t25 < dC * shat_sqrt * p[4][0] ) {
      ff[0] = 0;
      return 0;
	}
        
   double ME2 = (processID->*processID->matrixelementReal_HnonC)(p);

   // unleas using scheme like DRII the ME2 should be always positive
   if( ME2 < 0 || std::isnan(ME2) ) {
      std::cout << "Warning, negative ME2 " << ME2 << '\n';
      ff[0] = 0;
      return 0;      
   }
  
   // some final factors
   double fac = to_fb;
   fac /=  2 * shat;

   fac *= sin( xx[2] * pi ) * 4;
   fac /= 256 * pi_sqr;

   double pdf_flux = 0.0;
   for (const auto& f : processID->flav) {
      pdf_flux += f.at(2) * pdf->xfxQ( f.at(0), x1, mu_f ) * pdf->xfxQ( f.at(1), x2, mu_f );
   }
   pdf_flux /= x1 * x2;

   fac *=  pdf_flux;
   double xx0 = xx[5];
   double xx1 = xx[6];
   double xx2 = xx[0];

   /*    The integration over s35 and s45 is mapped on unit square spanned by [xx0, xx1]
    *
    *    Jakobian of transformations:
    *       xx0 -> s35min + (s35max - s35min) xx0
    *       xx1 -> s45min + (s45max - s45min) xx1
    */       
      double m {m1};
   double J = (Power(-4*Power(m,2) + S,2)*xx0*Sqrt(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1))*Power(-2*m + Sqrt(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1)),2)*xx2*
     Sqrt((-1 + xx2)*(S*xx0*xx1*(-1 + xx2) - 4*m*Sqrt(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1))*xx2 + Power(m,2)*(-4*xx0*xx1*(-1 + xx2) + 8*xx2))))/
   (4.*S*(-4*Power(m,2)*(-1 + xx0) + S*xx0)*(S*xx0*xx1*xx2 - 2*m*Sqrt(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1))*xx2 + Power(m,2)*(1 + (4 - 4*xx0*xx1)*xx2)));
   ME2 *= abs(J);

   double MassGlu = 2e+3;
   if (processID->matrixelementReal_HnonC_CSub1 != nullptr && shat_sqrt > m1 + MassGlu) {
      /*
      double mfort_1 = 1500;
      double mfort_2 = 1500;
      double mfort_3 = 0;
      double mfort_4 = 2000.;
      double p_to_fortran[20] = { 
         p[0][0], p[0][1], p[0][2], p[0][3],
         p[1][0], p[1][1], p[1][2], p[1][3],
         p[2][0], p[2][1], p[2][2], p[2][3],
         p[3][0], p[3][1], p[3][2], p[3][3],
         p[4][0], p[4][1], p[4][2], p[4][3],
      };
      std::vector< double* > p_mapped;
      double p_from_fortran[20];
      reshuffle_momenta_2__(p_to_fortran, &mfort_1, &mfort_2, &mfort_3, &mfort_4, p_from_fortran);
      p_mapped.push_back(new double[4]);
      p_mapped.push_back(new double[4]);
      p_mapped.push_back(new double[4]);
      p_mapped.push_back(new double[4]);
      p_mapped.push_back(new double[4]);
      for(int i=0; i<20; ++i) p_mapped[floor(i/4.)][i % 4] = p_from_fortran[i];
      */
      double J_s35mapped = -((m - MassGlu)*(m + MassGlu)*Power(-4*Power(m,2) + S,2)*xx0*Sqrt(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1))*
      Power(-2*m + Sqrt(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1)),2)*
      Sqrt(2*(m - MassGlu)*(m + MassGlu) + (-4*Power(m,2) + S)*xx0*xx1 - Power(Power(m,2) - Power(MassGlu,2),2)/(-(S*xx0*xx1) + 4*Power(m,2)*(-1 + xx0*xx1))))/
   (4.*Power(MassGlu,2)*S*(-4*Power(m,2)*(-1 + xx0) + S*xx0)*(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1) - 2*m*Sqrt(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1))));
      ME2 -= (processID->*processID->matrixelementReal_HnonC_CSub1)(p) * abs(J_s35mapped);
   }

   if(processID->matrixelementReal_HnonC_CSub2 != nullptr) {
      double J_s45mapped = -(xx0*xx2*(S*xx0*xx1 + (4 - 4*xx0*xx1)*pow(m1,2))*pow(S,-1)*pow(S - 4*pow(m1,2),2)*pow(S*xx0 - 4*(-1 + xx0)*pow(m1,2),-1)*
     (2*m1 - pow(S*xx0*xx1 + (4 - 4*xx0*xx1)*pow(m1,2),0.5))*
     (S*xx0*xx1 + (4 - 4*xx0*xx1)*pow(m1,2) - 2*m1*pow(S*xx0*xx1 + (4 - 4*xx0*xx1)*pow(m1,2),0.5))*
     pow((-1 + xx2)*(S*xx0*xx1*(-1 + xx2) + (-4*xx0*xx1*(-1 + xx2) + 8*xx2)*pow(m1,2) - 
         4*m1*xx2*pow(S*xx0*xx1 + (4 - 4*xx0*xx1)*pow(m1,2),0.5)),0.5)*
     pow(S*xx0*xx1*xx2 + (1 + (4 - 4*xx0*xx1)*xx2)*pow(m1,2) - 2*m1*xx2*pow(S*xx0*xx1 + (4 - 4*xx0*xx1)*pow(m1,2),0.5),-1));
      std::cout << " shouldn't be here\n";
      ME2 -= (processID->*processID->matrixelementReal_HnonC_CSub2)(p) * abs(J);
   }
   // delete (otherwise causes memory leak)
   for(std::vector<double*>::iterator i = p.begin(); i != p.end(); ++i) {
      delete (*i);
      *i = 0;
   }
   p.clear();
   p.shrink_to_fit();

   ff[0] = fac * ME2;

	return 0;
}
