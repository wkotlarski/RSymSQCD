//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "Process_uu_ulurg.h"
#include "HelAmps_MRSSMQCD_UFO_from_philip.h"
#include <iostream>
using namespace MG5_MRSSMQCD_UFO_from_philip; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: u u > ul ur g NP=1 WEIGHTED=3 @1

//--------------------------------------------------------------------------
// Initialize process.

void Process_uu_ulurg::initProc(string param_card_name, double mass_squark,
                                double mass_gluino, double alphaS) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MRSSMQCD_UFO_from_philip::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 

  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  /* begin hard code */
  pars->mdl_msul = mass_squark;
  pars->mdl_msur = mass_squark;
  pars->mdl_msdl = mass_squark;
  pars->mdl_msdr = mass_squark;
  pars->mdl_MD3  = mass_gluino;
  pars->aS = alphaS;
  /* end hard code */
  pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  //mME.push_back(pars->mdl_msul); 
  //mME.push_back(pars->mdl_msur); 
  /* begin hard code */
  mME.push_back(mass_squark); 
  mME.push_back(mass_squark); 
  /* end hard code */
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[4]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Process_uu_ulurg::sigmaKin() 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
    pars->printDependentParameters(); 
    pars->printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 4; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 8; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, 0, 0, -1}, {-1, -1,
      0, 0, 1}, {-1, 1, 0, 0, -1}, {-1, 1, 0, 0, 1}, {1, -1, 0, 0, -1}, {1, -1,
      0, 0, 1}, {1, 1, 0, 0, -1}, {1, 1, 0, 0, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {36}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_1_uu_ulurg(); 

        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_1_uu_ulurg(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }
  
  for (int i = 0; i < nprocesses; i++ ) {
    matrix_element[i] /= denominators[i];
  }



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double Process_uu_ulurg::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 2 && id2 == 2)
  {
    // Add matrix elements for processes with beams (2, 2)
    return matrix_element[0]; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void Process_uu_ulurg::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  sxxxxx(p[perm[2]], +1, w[2]); 
  sxxxxx(p[perm[3]], +1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  FFV1_2(w[0], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[5]); 
  FFS3C1_1(w[1], w[2], pars->GC_14, pars->mdl_MD3, pars->mdl_wdglu, w[6]); 
  FFS5C1_1(w[1], w[3], pars->GC_15, pars->mdl_MD3, pars->mdl_wdglu, w[7]); 
  VSS1_2(w[4], w[2], pars->GC_11, pars->mdl_msul, pars->ZERO, w[8]); 
  FFS5_2(w[0], w[3], pars->GC_15, pars->mdl_MD3, pars->mdl_wdglu, w[9]); 
  VSS1_2(w[4], w[3], pars->GC_11, pars->mdl_msur, pars->ZERO, w[10]); 
  FFS3C1_2(w[0], w[2], pars->GC_14, pars->mdl_MD3, pars->mdl_wdglu, w[11]); 
  FFV1C1_1(w[1], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[12]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFS5_0(w[5], w[6], w[3], pars->GC_15, amp[0]); 
  FFS3C1_0(w[5], w[7], w[2], pars->GC_14, amp[1]); 
  FFS3C1_0(w[9], w[1], w[8], pars->GC_14, amp[2]); 
  FFS3C1_0(w[0], w[7], w[8], pars->GC_14, amp[3]); 
  FFS5C1_0(w[11], w[1], w[10], pars->GC_15, amp[4]); 
  FFS5_0(w[0], w[6], w[10], pars->GC_15, amp[5]); 
  FFS5C1_0(w[11], w[12], w[3], pars->GC_15, amp[6]); 
  FFS3C1_0(w[9], w[12], w[2], pars->GC_14, amp[7]); 
  FFV1C1_0(w[11], w[7], w[4], pars->GC_13, amp[8]); 
  FFV1_0(w[9], w[6], w[4], pars->GC_13, amp[9]); 

}
double Process_uu_ulurg::matrix_1_uu_ulurg() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 10; 
  const int ncolor = 4; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1, 1, 1, 1}; 
  static const double cf[ncolor][ncolor] = {{12, 4, 4, 0}, {4, 12, 0, 4}, {4,
      0, 12, 4}, {0, 4, 4, 12}};

  // Calculate color flows
  jamp[0] = +1./2. * (-1./3. * amp[4] + amp[5] - 1./3. * amp[6] + amp[7] -
      std::complex<double> (0, 1) * amp[9]);
  jamp[1] = +1./2. * (-1./3. * amp[0] + amp[1] + amp[4] - 1./3. * amp[5] -
      std::complex<double> (0, 1) * amp[8]);
  jamp[2] = +1./2. * (-1./3. * amp[2] + amp[3] + amp[6] - 1./3. * amp[7] +
      std::complex<double> (0, 1) * amp[8]);
  jamp[3] = +1./2. * (+amp[0] - 1./3. * amp[1] + amp[2] - 1./3. * amp[3] +
      std::complex<double> (0, 1) * amp[9]);

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



