#include "XSection_HnonC.h"
#include <cassert>
#include "Process_uu_ulurg.cc"
#include <iostream>
Process_uu_ulurg process;

std::array<double, 3> XSection_HnonC::integrate() {

  process.initProc("/Users/Navir/Fizyka/Programy/"
  "MG5_aMC/PROC_SA_MRSSMQCD_UFO_from_philip_0/Cards/param_card.dat");
  
  //  integral dimension, number of integrands
  constexpr int ndim { 7 }, ncomp { 1 };
  //  accuraccy
  constexpr double accuracy_rel { 1e-3 }, 
          accuracy_abs { 1e-12 };

  constexpr int neval_min = 10000;
  long long int neval;
  constexpr long long int neval_max { 1000000000 }; 
    // @TODO: read from external source strtoll( "1e+3", NULL, 10 );

  // technical (Vegas specific) stuff
  constexpr int nstart = 200000;
  constexpr int nincrease = 100;
  constexpr int nbatch = 1000;
  constexpr int gridno = 0;
  const char* state_file = "";
  int nregions, fail;

  cubareal integral[ncomp], error[ncomp], prob[ncomp];
  llVegas( ndim, ncomp, integrand, NULL, 1,
           accuracy_rel, accuracy_abs, 8 | 1, 0,
           neval_min, neval_max, nstart, nincrease, nbatch,
           gridno, state_file, NULL,
           &neval, &fail, integral, error, prob );

  std::array <double, 3> result_finite
    { integral[0], error[0], prob[0] };
  return result_finite;
}

int XSection_HnonC::integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
  double m = squark_mass.at(1).at(1);
  double m_sqr = m*m;

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

	double x1 = 4 * m_sqr/S + (1 - 4 * m_sqr/S ) * xx[5];
	double x2 = 4 * m_sqr /(S * x1) + (1 - 4 * m_sqr/(S * x1)) * xx[6];
	double shat = x1 * x2 * S;
	double shat_sqrt = sqrt( shat );

	double Ej_max = shat_sqrt/2 - 2 * m_sqr/shat_sqrt;

	if ( Ej_max < dS * shat_sqrt/2) {
	  ff[0] = 0;
	  return 1;
	}

	double Ej = dS * shat_sqrt/2 + ( Ej_max - dS * shat_sqrt/2) * xx[0];

	double c = shat - 2 * shat_sqrt * Ej;
	// Eq. 4.5
	double E1_max = ( shat_sqrt - Ej) * c + Ej * sqrt( (c - 2 * m_sqr) 
    * (c - 2 * m_sqr) - 4 * m_sqr * m_sqr );
	E1_max /= 2 * c;
  double E1_min = ( shat_sqrt - Ej) * c - Ej * sqrt( (c - 2 * m_sqr) 
    * (c - 2 * m_sqr) - 4 * m_sqr * m_sqr );
  E1_min /= 2 * c;
	double E1 = E1_min + (E1_max - E1_min) * xx[1];

	// Eq. 4.2 with E2 = Ej
  double cosx = (shat - 2 * shat_sqrt * ( E1 + Ej ) + 2 * Ej * E1 )/
      (2 * Ej * sqrt( ( E1 - m ) * ( E1 + m ) ) );

  // check if due to numerics abs(cosx) is not > 1
  if ( cosx > 1 || cosx < -1)  {
    //cout << "Warning! cos(x) = " << cosx << " out of bounds. Setting it +/- 1." <<  endl;
    cosx = cosx > 0 ? 1.0 : -1.0;
  }

  // initialize matrix of particles momenta
  std::vector< double* > p(1, new double [4]);

  // incoming parton 1 momenta
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
  double p1 = sqrt( ( E1 - m ) * ( E1 + m ) );
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
  if ( pi * xx[2] + acos(cosx) < pi ) {
    parton_theta = pi * xx[2] + acos(cosx);
    parton_phi = two_pi * xx[3];
  }
  else if (pi * xx[2] - acos(cosx) > 0 ) {
    parton_theta = pi * xx[2] - acos(cosx);
    parton_phi = two_pi * xx[3];
  }
  else {
    parton_theta = abs(pi * xx[2] - acos(cosx));
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
  std::vector< double* > p_temp_2(1, new double [4]);
  p_temp_2[0][0] = Ej;
  p_temp_2[0][1] = p_temp.x();
  p_temp_2[0][2] = p_temp.y();
  p_temp_2[0][3] = p_temp.z();

  // 2nd sgluon momenta
  p.push_back(new double[4]);
  p[3][0] = shat_sqrt - E1 - Ej;
  for ( int i = 1; i < 4; ++i) {
    p[3][i] = - p[2][i] - p_temp_2[0][i];
  }

  // write parton momentum to momentum matrix p
  p.push_back(new double[4]);
  p[4] = p_temp_2[0];    
      
  double t15 = p[0][0] * p[4][0] - p[0][3] * p[4][3];
  t15 = - 2 * t15;
  double t25 = p[1][0] * p[4][0] - p[1][3] * p[4][3];
  t25 = - 2 * t25;

  // check if we are not in the collinear region
  // if yes, return
  if ( -t15 < dC * shat || -t25 < dC * shat ) {
	  ff[0] = 0;
	  return 0;
	}

  process.setMomenta(p);	// Set momenta for this event
  process.sigmaKin();		// Evaluate matrix element
  const double* matrix_elements = process.getMatrixElements();

  // some final factors
  double temp = to_fb * matrix_elements[0];
  temp /=  2 * shat;

  temp *= sin( xx[2] * pi ) * 4;
  temp /= 256 * pi_sqr;

  // choose gg or qqbar initial state
  double temp2 = 0;
  //temp2 = pdf_nlo->xfxQ(21, x1, m)/x1 * pdf_nlo->xfxQ(21, x2, m)/x2;
  for ( int flav = 1; flav < 6; ++flav ) {
    //temp2 += 2 * pdf->xfxQ(flav, x1, m)/x1 * pdf->xfxQ(-flav, x2, m)/x2;
  }
  temp2 = pdf_nlo->xfxQ(2, x1, m)/x1 * pdf_nlo->xfxQ(2, x2, m)/x2;
  //temp2 = 2 * pdf->xfxQ(21, x1, m)/x1 * pdf->xfxQ(2, x2, m)/x2;

  temp *=  temp2;
  double xx0 = xx[5];
  double xx1 = xx[6];
  double xx2 = xx[0];

  // TODO: simplify this
  double jacobian =  ((Power(-4*Power(m,2) + S,2)*xx0*((-4*Power(m,2) + S)*xx0*xx1 +
        dS*(-(S*xx0*xx1) + 4*Power(m,2)*(-1 + xx0*xx1)))*
	      (dS*(-(S*xx0*xx1) + 4*Power(m,2)*(-1 + xx0*xx1))*(-1 + xx2) +
	        (-4*Power(m,2) + S)*xx0*xx1*xx2)*
	      Sqrt(((-4*Power(m,2) + S)*xx0*xx1 +
	          dS*(-(S*xx0*xx1) + 4*Power(m,2)*(-1 + xx0*xx1)))*(-1 + xx2)*
	        (S*xx0*xx1*(-1 + dS + xx2 - dS*xx2) +
	          4*Power(m,2)*(-1 + xx0*xx1 + dS*(-1 + xx0*xx1)*(-1 + xx2) - xx0*xx1*xx2))))/
	    (4.*S*(-4*Power(m,2)*(-1 + xx0) + S*xx0)*(S*xx0*xx1 + Power(m,2)*(4 - 4*xx0*xx1))*
	      ((-1 + dS)*S*xx0*xx1*(-1 + xx2) -
	        4*Power(m,2)*(-1 + xx0*xx1 + dS*(-1 + xx0*xx1)*(-1 + xx2) - xx0*xx1*xx2))));

    ff[0] = temp * abs(jacobian);

	return 0;
}
