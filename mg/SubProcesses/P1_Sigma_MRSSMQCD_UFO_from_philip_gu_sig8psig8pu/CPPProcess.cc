//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "CPPProcess.h"
#include "HelAmps_MRSSMQCD_UFO_from_philip.h"

using namespace MG5_MRSSMQCD_UFO_from_philip; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g u > sig8p sig8p u NP=1 WEIGHTED=3 @1
// Process: g c > sig8p sig8p c NP=1 WEIGHTED=3 @1
// Process: g d > sig8p sig8p d NP=1 WEIGHTED=3 @1
// Process: g s > sig8p sig8p s NP=1 WEIGHTED=3 @1
// Process: g b > sig8p sig8p b NP=1 WEIGHTED=3 @1

//--------------------------------------------------------------------------
// Initialize process.

void CPPProcess::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_MRSSMQCD_UFO_from_philip::getInstance(); 
  SLHAReader slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
  pars->printIndependentParameters(); 
  pars->printIndependentCouplings(); 
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->mdl_mOp); 
  mME.push_back(pars->mdl_mOp); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[6]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void CPPProcess::sigmaKin() 
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
  for(int i = 0; i < 6; i++ )
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
  const int denominators[nprocesses] = {192, 192}; 

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
        t[0] = matrix_1_gu_sig8psig8pu(); 
        // Mirror initial state momenta for mirror process
        perm[0] = 1; 
        perm[1] = 0; 
        // Calculate wavefunctions
        calculate_wavefunctions(perm, helicities[ihel]); 
        // Mirror back
        perm[0] = 0; 
        perm[1] = 1; 
        // Calculate matrix elements
        t[1] = matrix_1_gu_sig8psig8pu(); 
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
      t[0] = matrix_1_gu_sig8psig8pu(); 
      // Mirror initial state momenta for mirror process
      perm[0] = 1; 
      perm[1] = 0; 
      // Calculate wavefunctions
      calculate_wavefunctions(perm, helicities[ihel]); 
      // Mirror back
      perm[0] = 0; 
      perm[1] = 1; 
      // Calculate matrix elements
      t[1] = matrix_1_gu_sig8psig8pu(); 
      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double CPPProcess::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 1 && id2 == 21)
  {
    // Add matrix elements for processes with beams (1, 21)
    return matrix_element[1]; 
  }
  else if(id1 == 21 && id2 == 1)
  {
    // Add matrix elements for processes with beams (21, 1)
    return matrix_element[0]; 
  }
  else if(id1 == 2 && id2 == 21)
  {
    // Add matrix elements for processes with beams (2, 21)
    return matrix_element[1]; 
  }
  else if(id1 == 21 && id2 == 2)
  {
    // Add matrix elements for processes with beams (21, 2)
    return matrix_element[0]; 
  }
  else if(id1 == 21 && id2 == 3)
  {
    // Add matrix elements for processes with beams (21, 3)
    return matrix_element[0]; 
  }
  else if(id1 == 21 && id2 == 4)
  {
    // Add matrix elements for processes with beams (21, 4)
    return matrix_element[0]; 
  }
  else if(id1 == 21 && id2 == 5)
  {
    // Add matrix elements for processes with beams (21, 5)
    return matrix_element[0]; 
  }
  else if(id1 == 4 && id2 == 21)
  {
    // Add matrix elements for processes with beams (4, 21)
    return matrix_element[1]; 
  }
  else if(id1 == 3 && id2 == 21)
  {
    // Add matrix elements for processes with beams (3, 21)
    return matrix_element[1]; 
  }
  else if(id1 == 5 && id2 == 21)
  {
    // Add matrix elements for processes with beams (5, 21)
    return matrix_element[1]; 
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

void CPPProcess::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  ixxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]); 
  sxxxxx(p[perm[2]], +1, w[2]); 
  sxxxxx(p[perm[3]], +1, w[3]); 
  oxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  FFV1_2(w[1], w[0], pars->GC_12, pars->ZERO, pars->ZERO, w[5]); 
  VSS1P0_1(w[2], w[3], pars->GC_13, pars->ZERO, pars->ZERO, w[6]); 
  VSS1_2(w[0], w[2], pars->GC_13, pars->mdl_mOp, pars->mdl_wOp, w[7]); 
  FFV1P0_3(w[1], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[8]); 
  VSS1_2(w[0], w[3], pars->GC_13, pars->mdl_mOp, pars->mdl_wOp, w[9]); 
  FFV1_1(w[4], w[0], pars->GC_12, pars->ZERO, pars->ZERO, w[10]); 
  VVSS1P0_1(w[0], w[2], w[3], pars->GC_17, pars->ZERO, pars->ZERO, w[11]); 
  VVSS1P0_1(w[0], w[2], w[3], pars->GC_17, pars->ZERO, pars->ZERO, w[12]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[5], w[4], w[6], pars->GC_12, amp[0]); 
  VSS1_0(w[8], w[7], w[3], pars->GC_13, amp[1]); 
  VSS1_0(w[8], w[9], w[2], pars->GC_13, amp[2]); 
  FFV1_0(w[1], w[10], w[6], pars->GC_12, amp[3]); 
  VVV1_0(w[0], w[8], w[6], pars->GC_10, amp[4]); 
  FFV1_0(w[1], w[4], w[11], pars->GC_12, amp[5]); 
  FFV1_0(w[1], w[4], w[12], pars->GC_12, amp[6]); 

}
double CPPProcess::matrix_1_gu_sig8psig8pu() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 7; 
  const int ncolor = 6; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {9, 9, 9, 9, 9, 9}; 
  static const double cf[ncolor][ncolor] = {{64, -8, -8, 1, 1, 10}, {-8, 64, 1,
      10, -8, 1}, {-8, 1, 64, -8, 10, 1}, {1, 10, -8, 64, 1, -8}, {1, -8, 10,
      1, 64, -8}, {10, 1, 1, -8, -8, 64}};

  // Calculate color flows
  jamp[0] = +amp[1] - std::complex<double> (0, 1) * amp[3] + amp[4] + amp[6]; 
  jamp[1] = +amp[2] + std::complex<double> (0, 1) * amp[3] - amp[4] + amp[5]; 
  jamp[2] = -amp[1] - amp[2] - amp[6] - amp[5]; 
  jamp[3] = -std::complex<double> (0, 1) * amp[0] + amp[2] - amp[4] + amp[5]; 
  jamp[4] = -amp[1] - amp[2] - amp[6] - amp[5]; 
  jamp[5] = +std::complex<double> (0, 1) * amp[0] + amp[1] + amp[4] + amp[6]; 

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



