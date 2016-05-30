//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "Process_gg_s8ps8pg.h"
#include "HelAmps_MRSSMQCD_UFO_from_philip.h"

using namespace MG5_MRSSMQCD_UFO_from_philip; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: g g > sig8p sig8p g NP=1 WEIGHTED=3 @1

//--------------------------------------------------------------------------
// Initialize process.

void Process_gg_s8ps8pg::initProc(string param_card_name) 
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
  jamp2[0] = new double[24]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void Process_gg_s8ps8pg::sigmaKin() 
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
  for(int i = 0; i < 24; i++ )
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
  const int denominators[nprocesses] = {512}; 

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
        t[0] = matrix_1_gg_sig8psig8pg(); 

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
      t[0] = matrix_1_gg_sig8psig8pg(); 

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

double Process_gg_s8ps8pg::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 21 && id2 == 21)
  {
    // Add matrix elements for processes with beams (21, 21)
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

void Process_gg_s8ps8pg::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  sxxxxx(p[perm[2]], +1, w[2]); 
  sxxxxx(p[perm[3]], +1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  VVV1P0_1(w[0], w[1], pars->GC_10, pars->ZERO, pars->ZERO, w[5]); 
  VSS1P0_1(w[2], w[3], pars->GC_13, pars->ZERO, pars->ZERO, w[6]); 
  VSS1_2(w[4], w[2], pars->GC_13, pars->mdl_mOp, pars->mdl_wOp, w[7]); 
  VSS1_2(w[4], w[3], pars->GC_13, pars->mdl_mOp, pars->mdl_wOp, w[8]); 
  VSS1_2(w[0], w[2], pars->GC_13, pars->mdl_mOp, pars->mdl_wOp, w[9]); 
  VSS1_2(w[1], w[3], pars->GC_13, pars->mdl_mOp, pars->mdl_wOp, w[10]); 
  VVV1P0_1(w[1], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[11]); 
  VSS1_2(w[0], w[3], pars->GC_13, pars->mdl_mOp, pars->mdl_wOp, w[12]); 
  VSS1_2(w[1], w[2], pars->GC_13, pars->mdl_mOp, pars->mdl_wOp, w[13]); 
  VVV1P0_1(w[0], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[14]); 
  VVSS1_3(w[0], w[1], w[2], pars->GC_17, pars->mdl_mOp, pars->mdl_wOp, w[15]); 
  VVSS1_3(w[0], w[1], w[2], pars->GC_17, pars->mdl_mOp, pars->mdl_wOp, w[16]); 
  VVSS1_3(w[0], w[1], w[3], pars->GC_17, pars->mdl_mOp, pars->mdl_wOp, w[17]); 
  VVSS1_3(w[0], w[1], w[3], pars->GC_17, pars->mdl_mOp, pars->mdl_wOp, w[18]); 
  VVVV1P0_1(w[0], w[1], w[4], pars->GC_17, pars->ZERO, pars->ZERO, w[19]); 
  VVVV3P0_1(w[0], w[1], w[4], pars->GC_17, pars->ZERO, pars->ZERO, w[20]); 
  VVVV4P0_1(w[0], w[1], w[4], pars->GC_17, pars->ZERO, pars->ZERO, w[21]); 
  VVSS1P0_1(w[0], w[2], w[3], pars->GC_17, pars->ZERO, pars->ZERO, w[22]); 
  VVSS1P0_1(w[0], w[2], w[3], pars->GC_17, pars->ZERO, pars->ZERO, w[23]); 
  VVSS1_3(w[0], w[4], w[2], pars->GC_17, pars->mdl_mOp, pars->mdl_wOp, w[24]); 
  VVSS1_3(w[0], w[4], w[2], pars->GC_17, pars->mdl_mOp, pars->mdl_wOp, w[25]); 
  VVSS1_3(w[0], w[4], w[3], pars->GC_17, pars->mdl_mOp, pars->mdl_wOp, w[26]); 
  VVSS1_3(w[0], w[4], w[3], pars->GC_17, pars->mdl_mOp, pars->mdl_wOp, w[27]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVV1_0(w[5], w[6], w[4], pars->GC_10, amp[0]); 
  VSS1_0(w[5], w[7], w[3], pars->GC_13, amp[1]); 
  VSS1_0(w[5], w[2], w[8], pars->GC_13, amp[2]); 
  VVSS1_0(w[4], w[5], w[2], w[3], pars->GC_17, amp[3]); 
  VVSS1_0(w[4], w[5], w[2], w[3], pars->GC_17, amp[4]); 
  VSS1_0(w[4], w[9], w[10], pars->GC_13, amp[5]); 
  VSS1_0(w[11], w[9], w[3], pars->GC_13, amp[6]); 
  VSS1_0(w[1], w[9], w[8], pars->GC_13, amp[7]); 
  VVSS1_0(w[1], w[4], w[3], w[9], pars->GC_17, amp[8]); 
  VVSS1_0(w[1], w[4], w[3], w[9], pars->GC_17, amp[9]); 
  VSS1_0(w[4], w[12], w[13], pars->GC_13, amp[10]); 
  VSS1_0(w[11], w[12], w[2], pars->GC_13, amp[11]); 
  VSS1_0(w[1], w[12], w[7], pars->GC_13, amp[12]); 
  VVSS1_0(w[1], w[4], w[2], w[12], pars->GC_17, amp[13]); 
  VVSS1_0(w[1], w[4], w[2], w[12], pars->GC_17, amp[14]); 
  VSS1_0(w[14], w[13], w[3], pars->GC_13, amp[15]); 
  VSS1_0(w[14], w[10], w[2], pars->GC_13, amp[16]); 
  VVV1_0(w[14], w[1], w[6], pars->GC_10, amp[17]); 
  VVSS1_0(w[1], w[14], w[2], w[3], pars->GC_17, amp[18]); 
  VVSS1_0(w[1], w[14], w[2], w[3], pars->GC_17, amp[19]); 
  VSS1_0(w[0], w[13], w[8], pars->GC_13, amp[20]); 
  VSS1_0(w[0], w[10], w[7], pars->GC_13, amp[21]); 
  VVV1_0(w[0], w[11], w[6], pars->GC_10, amp[22]); 
  VSS1_0(w[4], w[3], w[15], pars->GC_13, amp[23]); 
  VSS1_0(w[4], w[3], w[16], pars->GC_13, amp[24]); 
  VSS1_0(w[4], w[2], w[17], pars->GC_13, amp[25]); 
  VSS1_0(w[4], w[2], w[18], pars->GC_13, amp[26]); 
  VSS1_0(w[19], w[2], w[3], pars->GC_13, amp[27]); 
  VSS1_0(w[20], w[2], w[3], pars->GC_13, amp[28]); 
  VSS1_0(w[21], w[2], w[3], pars->GC_13, amp[29]); 
  VVV1_0(w[1], w[4], w[22], pars->GC_10, amp[30]); 
  VVV1_0(w[1], w[4], w[23], pars->GC_10, amp[31]); 
  VSS1_0(w[1], w[3], w[24], pars->GC_13, amp[32]); 
  VSS1_0(w[1], w[3], w[25], pars->GC_13, amp[33]); 
  VSS1_0(w[1], w[2], w[26], pars->GC_13, amp[34]); 
  VSS1_0(w[1], w[2], w[27], pars->GC_13, amp[35]); 

}
double Process_gg_s8ps8pg::matrix_1_gg_sig8psig8pg() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 36; 
  const int ncolor = 24; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {108, 108, 108, 108, 108, 108, 108, 108,
      108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108,
      108, 108};
  static const double cf[ncolor][ncolor] = {{455, -58, -58, 14, 14, 68, -58,
      -4, 14, -58, 5, -4, 14, 5, 68, -4, 14, 68, -58, -4, -4, 68, 68, -40},
      {-58, 455, 14, 68, -58, 14, -4, -58, 5, -4, 14, -58, -58, -4, -4, 68, 68,
      -40, 14, 5, 68, -4, 14, 68}, {-58, 14, 455, -58, 68, 14, 14, 5, 68, -4,
      14, 68, -58, -4, 14, -58, 5, -4, -4, -58, 68, -40, -4, 68}, {14, 68, -58,
      455, 14, -58, -58, -4, -4, 68, 68, -40, -4, -58, 5, -4, 14, -58, 5, 14,
      14, 68, 68, -4}, {14, -58, 68, 14, 455, -58, 5, 14, 14, 68, 68, -4, -4,
      -58, 68, -40, -4, 68, -58, -4, 14, -58, 5, -4}, {68, 14, 14, -58, -58,
      455, -4, -58, 68, -40, -4, 68, 5, 14, 14, 68, 68, -4, -4, -58, 5, -4, 14,
      -58}, {-58, -4, 14, -58, 5, -4, 455, -58, -58, 14, 14, 68, 68, -4, 14, 5,
      68, 14, -4, 68, -58, -4, -40, 68}, {-4, -58, 5, -4, 14, -58, -58, 455,
      14, 68, -58, 14, -4, 68, -58, -4, -40, 68, 68, -4, 14, 5, 68, 14}, {14,
      5, 68, -4, 14, 68, -58, 14, 455, -58, 68, 14, 14, -58, -58, -4, -4, 5,
      68, -40, -4, -58, 68, -4}, {-58, -4, -4, 68, 68, -40, 14, 68, -58, 455,
      14, -58, 5, -4, -4, -58, -58, 14, 14, 68, 5, 14, -4, 68}, {5, 14, 14, 68,
      68, -4, 14, -58, 68, 14, 455, -58, 68, -40, -4, -58, 68, -4, 14, -58,
      -58, -4, -4, 5}, {-4, -58, 68, -40, -4, 68, 68, 14, 14, -58, -58, 455,
      14, 68, 5, 14, -4, 68, 5, -4, -4, -58, -58, 14}, {14, -58, -58, -4, -4,
      5, 68, -4, 14, 5, 68, 14, 455, -58, -58, 14, 14, 68, 68, -4, -40, 68,
      -58, -4}, {5, -4, -4, -58, -58, 14, -4, 68, -58, -4, -40, 68, -58, 455,
      14, 68, -58, 14, -4, 68, 68, 14, 14, 5}, {68, -4, 14, 5, 68, 14, 14, -58,
      -58, -4, -4, 5, -58, 14, 455, -58, 68, 14, -40, 68, 68, -4, -4, -58},
      {-4, 68, -58, -4, -40, 68, 5, -4, -4, -58, -58, 14, 14, 68, -58, 455, 14,
      -58, 68, 14, -4, 68, 5, 14}, {14, 68, 5, 14, -4, 68, 68, -40, -4, -58,
      68, -4, 14, -58, 68, 14, 455, -58, -58, 14, -4, 5, -58, -4}, {68, -40,
      -4, -58, 68, -4, 14, 68, 5, 14, -4, 68, 68, 14, 14, -58, -58, 455, -4, 5,
      -58, 14, -4, -58}, {-58, 14, -4, 5, -58, -4, -4, 68, 68, 14, 14, 5, 68,
      -4, -40, 68, -58, -4, 455, -58, -58, 14, 14, 68}, {-4, 5, -58, 14, -4,
      -58, 68, -4, -40, 68, -58, -4, -4, 68, 68, 14, 14, 5, -58, 455, 14, 68,
      -58, 14}, {-4, 68, 68, 14, 14, 5, -58, 14, -4, 5, -58, -4, -40, 68, 68,
      -4, -4, -58, -58, 14, 455, -58, 68, 14}, {68, -4, -40, 68, -58, -4, -4,
      5, -58, 14, -4, -58, 68, 14, -4, 68, 5, 14, 14, 68, -58, 455, 14, -58},
      {68, 14, -4, 68, 5, 14, -40, 68, 68, -4, -4, -58, -58, 14, -4, 5, -58,
      -4, 14, -58, 68, 14, 455, -58}, {-40, 68, 68, -4, -4, -58, 68, 14, -4,
      68, 5, 14, -4, 5, -58, 14, -4, -58, 68, 14, 14, -58, -58, 455}};

  // Calculate color flows
  jamp[0] = +2. * (+std::complex<double> (0, 1) * amp[0] + std::complex<double>
      (0, 1) * amp[2] - std::complex<double> (0, 1) * amp[4] +
      std::complex<double> (0, 1) * amp[15] - std::complex<double> (0, 1) *
      amp[17] + std::complex<double> (0, 1) * amp[18] - std::complex<double>
      (0, 1) * amp[20] + std::complex<double> (0, 1) * amp[23] +
      std::complex<double> (0, 1) * amp[29] + std::complex<double> (0, 1) *
      amp[28] - std::complex<double> (0, 1) * amp[34]);
  jamp[1] = +2. * (+std::complex<double> (0, 1) * amp[1] - std::complex<double>
      (0, 1) * amp[2] + std::complex<double> (0, 1) * amp[4] +
      std::complex<double> (0, 1) * amp[3] - std::complex<double> (0, 1) *
      amp[10] - std::complex<double> (0, 1) * amp[12] + std::complex<double>
      (0, 1) * amp[14] + std::complex<double> (0, 1) * amp[13] +
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[23] + std::complex<double> (0, 1) * amp[26] + std::complex<double>
      (0, 1) * amp[35] + std::complex<double> (0, 1) * amp[34]);
  jamp[2] = +2. * (-std::complex<double> (0, 1) * amp[0] - std::complex<double>
      (0, 1) * amp[1] - std::complex<double> (0, 1) * amp[3] +
      std::complex<double> (0, 1) * amp[16] + std::complex<double> (0, 1) *
      amp[17] + std::complex<double> (0, 1) * amp[19] - std::complex<double>
      (0, 1) * amp[21] + std::complex<double> (0, 1) * amp[25] -
      std::complex<double> (0, 1) * amp[29] - std::complex<double> (0, 1) *
      amp[28] - std::complex<double> (0, 1) * amp[32]);
  jamp[3] = +2. * (+std::complex<double> (0, 1) * amp[1] - std::complex<double>
      (0, 1) * amp[2] + std::complex<double> (0, 1) * amp[4] +
      std::complex<double> (0, 1) * amp[3] - std::complex<double> (0, 1) *
      amp[5] - std::complex<double> (0, 1) * amp[7] + std::complex<double> (0,
      1) * amp[9] + std::complex<double> (0, 1) * amp[8] + std::complex<double>
      (0, 1) * amp[21] + std::complex<double> (0, 1) * amp[24] -
      std::complex<double> (0, 1) * amp[25] + std::complex<double> (0, 1) *
      amp[33] + std::complex<double> (0, 1) * amp[32]);
  jamp[4] = +2. * (-std::complex<double> (0, 1) * amp[0] - std::complex<double>
      (0, 1) * amp[1] - std::complex<double> (0, 1) * amp[3] -
      std::complex<double> (0, 1) * amp[11] + std::complex<double> (0, 1) *
      amp[12] - std::complex<double> (0, 1) * amp[14] + std::complex<double>
      (0, 1) * amp[22] - std::complex<double> (0, 1) * amp[26] -
      std::complex<double> (0, 1) * amp[29] + std::complex<double> (0, 1) *
      amp[27] - std::complex<double> (0, 1) * amp[30]);
  jamp[5] = +2. * (+std::complex<double> (0, 1) * amp[0] + std::complex<double>
      (0, 1) * amp[2] - std::complex<double> (0, 1) * amp[4] -
      std::complex<double> (0, 1) * amp[6] + std::complex<double> (0, 1) *
      amp[7] - std::complex<double> (0, 1) * amp[9] - std::complex<double> (0,
      1) * amp[22] - std::complex<double> (0, 1) * amp[24] +
      std::complex<double> (0, 1) * amp[29] - std::complex<double> (0, 1) *
      amp[27] - std::complex<double> (0, 1) * amp[31]);
  jamp[6] = +2. * (+std::complex<double> (0, 1) * amp[5] + std::complex<double>
      (0, 1) * amp[7] - std::complex<double> (0, 1) * amp[9] -
      std::complex<double> (0, 1) * amp[8] - std::complex<double> (0, 1) *
      amp[15] - std::complex<double> (0, 1) * amp[16] - std::complex<double>
      (0, 1) * amp[19] - std::complex<double> (0, 1) * amp[18] +
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[24] - std::complex<double> (0, 1) * amp[23] - std::complex<double>
      (0, 1) * amp[33] + std::complex<double> (0, 1) * amp[34]);
  jamp[7] = +2. * (+std::complex<double> (0, 1) * amp[6] - std::complex<double>
      (0, 1) * amp[7] + std::complex<double> (0, 1) * amp[9] +
      std::complex<double> (0, 1) * amp[10] + std::complex<double> (0, 1) *
      amp[11] - std::complex<double> (0, 1) * amp[13] - std::complex<double>
      (0, 1) * amp[20] + std::complex<double> (0, 1) * amp[24] +
      std::complex<double> (0, 1) * amp[23] + std::complex<double> (0, 1) *
      amp[31] + std::complex<double> (0, 1) * amp[30] - std::complex<double>
      (0, 1) * amp[35] - std::complex<double> (0, 1) * amp[34]);
  jamp[8] = +2. * (-std::complex<double> (0, 1) * amp[5] - std::complex<double>
      (0, 1) * amp[6] + std::complex<double> (0, 1) * amp[8] +
      std::complex<double> (0, 1) * amp[16] + std::complex<double> (0, 1) *
      amp[17] + std::complex<double> (0, 1) * amp[19] - std::complex<double>
      (0, 1) * amp[22] - std::complex<double> (0, 1) * amp[28] -
      std::complex<double> (0, 1) * amp[27] - std::complex<double> (0, 1) *
      amp[31] + std::complex<double> (0, 1) * amp[33]);
  jamp[9] = +2. * (-std::complex<double> (0, 1) * amp[0] - std::complex<double>
      (0, 1) * amp[2] + std::complex<double> (0, 1) * amp[4] +
      std::complex<double> (0, 1) * amp[6] - std::complex<double> (0, 1) *
      amp[7] + std::complex<double> (0, 1) * amp[9] + std::complex<double> (0,
      1) * amp[22] + std::complex<double> (0, 1) * amp[24] -
      std::complex<double> (0, 1) * amp[29] + std::complex<double> (0, 1) *
      amp[27] + std::complex<double> (0, 1) * amp[31]);
  jamp[10] = +2. * (-std::complex<double> (0, 1) * amp[5] -
      std::complex<double> (0, 1) * amp[6] + std::complex<double> (0, 1) *
      amp[8] - std::complex<double> (0, 1) * amp[11] + std::complex<double> (0,
      1) * amp[12] - std::complex<double> (0, 1) * amp[14] +
      std::complex<double> (0, 1) * amp[21] - std::complex<double> (0, 1) *
      amp[26] - std::complex<double> (0, 1) * amp[25] - std::complex<double>
      (0, 1) * amp[31] - std::complex<double> (0, 1) * amp[30] +
      std::complex<double> (0, 1) * amp[33] + std::complex<double> (0, 1) *
      amp[32]);
  jamp[11] = +2. * (-std::complex<double> (0, 1) * amp[1] +
      std::complex<double> (0, 1) * amp[2] - std::complex<double> (0, 1) *
      amp[4] - std::complex<double> (0, 1) * amp[3] + std::complex<double> (0,
      1) * amp[5] + std::complex<double> (0, 1) * amp[7] - std::complex<double>
      (0, 1) * amp[9] - std::complex<double> (0, 1) * amp[8] -
      std::complex<double> (0, 1) * amp[21] - std::complex<double> (0, 1) *
      amp[24] + std::complex<double> (0, 1) * amp[25] - std::complex<double>
      (0, 1) * amp[33] - std::complex<double> (0, 1) * amp[32]);
  jamp[12] = +2. * (+std::complex<double> (0, 1) * amp[10] +
      std::complex<double> (0, 1) * amp[12] - std::complex<double> (0, 1) *
      amp[14] - std::complex<double> (0, 1) * amp[13] - std::complex<double>
      (0, 1) * amp[15] - std::complex<double> (0, 1) * amp[16] -
      std::complex<double> (0, 1) * amp[19] - std::complex<double> (0, 1) *
      amp[18] + std::complex<double> (0, 1) * amp[21] - std::complex<double>
      (0, 1) * amp[26] - std::complex<double> (0, 1) * amp[25] +
      std::complex<double> (0, 1) * amp[32] - std::complex<double> (0, 1) *
      amp[35]);
  jamp[13] = +2. * (+std::complex<double> (0, 1) * amp[5] +
      std::complex<double> (0, 1) * amp[6] - std::complex<double> (0, 1) *
      amp[8] + std::complex<double> (0, 1) * amp[11] - std::complex<double> (0,
      1) * amp[12] + std::complex<double> (0, 1) * amp[14] -
      std::complex<double> (0, 1) * amp[21] + std::complex<double> (0, 1) *
      amp[26] + std::complex<double> (0, 1) * amp[25] + std::complex<double>
      (0, 1) * amp[31] + std::complex<double> (0, 1) * amp[30] -
      std::complex<double> (0, 1) * amp[33] - std::complex<double> (0, 1) *
      amp[32]);
  jamp[14] = +2. * (-std::complex<double> (0, 1) * amp[10] -
      std::complex<double> (0, 1) * amp[11] + std::complex<double> (0, 1) *
      amp[13] + std::complex<double> (0, 1) * amp[15] - std::complex<double>
      (0, 1) * amp[17] + std::complex<double> (0, 1) * amp[18] +
      std::complex<double> (0, 1) * amp[22] + std::complex<double> (0, 1) *
      amp[28] + std::complex<double> (0, 1) * amp[27] - std::complex<double>
      (0, 1) * amp[30] + std::complex<double> (0, 1) * amp[35]);
  jamp[15] = +2. * (+std::complex<double> (0, 1) * amp[0] +
      std::complex<double> (0, 1) * amp[1] + std::complex<double> (0, 1) *
      amp[3] + std::complex<double> (0, 1) * amp[11] - std::complex<double> (0,
      1) * amp[12] + std::complex<double> (0, 1) * amp[14] -
      std::complex<double> (0, 1) * amp[22] + std::complex<double> (0, 1) *
      amp[26] + std::complex<double> (0, 1) * amp[29] - std::complex<double>
      (0, 1) * amp[27] + std::complex<double> (0, 1) * amp[30]);
  jamp[16] = +2. * (-std::complex<double> (0, 1) * amp[6] +
      std::complex<double> (0, 1) * amp[7] - std::complex<double> (0, 1) *
      amp[9] - std::complex<double> (0, 1) * amp[10] - std::complex<double> (0,
      1) * amp[11] + std::complex<double> (0, 1) * amp[13] +
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[24] - std::complex<double> (0, 1) * amp[23] - std::complex<double>
      (0, 1) * amp[31] - std::complex<double> (0, 1) * amp[30] +
      std::complex<double> (0, 1) * amp[35] + std::complex<double> (0, 1) *
      amp[34]);
  jamp[17] = +2. * (-std::complex<double> (0, 1) * amp[1] +
      std::complex<double> (0, 1) * amp[2] - std::complex<double> (0, 1) *
      amp[4] - std::complex<double> (0, 1) * amp[3] + std::complex<double> (0,
      1) * amp[10] + std::complex<double> (0, 1) * amp[12] -
      std::complex<double> (0, 1) * amp[14] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[20] + std::complex<double>
      (0, 1) * amp[23] - std::complex<double> (0, 1) * amp[26] -
      std::complex<double> (0, 1) * amp[35] - std::complex<double> (0, 1) *
      amp[34]);
  jamp[18] = +2. * (+std::complex<double> (0, 1) * amp[10] +
      std::complex<double> (0, 1) * amp[11] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[15] + std::complex<double>
      (0, 1) * amp[17] - std::complex<double> (0, 1) * amp[18] -
      std::complex<double> (0, 1) * amp[22] - std::complex<double> (0, 1) *
      amp[28] - std::complex<double> (0, 1) * amp[27] + std::complex<double>
      (0, 1) * amp[30] - std::complex<double> (0, 1) * amp[35]);
  jamp[19] = +2. * (+std::complex<double> (0, 1) * amp[5] +
      std::complex<double> (0, 1) * amp[6] - std::complex<double> (0, 1) *
      amp[8] - std::complex<double> (0, 1) * amp[16] - std::complex<double> (0,
      1) * amp[17] - std::complex<double> (0, 1) * amp[19] +
      std::complex<double> (0, 1) * amp[22] + std::complex<double> (0, 1) *
      amp[28] + std::complex<double> (0, 1) * amp[27] + std::complex<double>
      (0, 1) * amp[31] - std::complex<double> (0, 1) * amp[33]);
  jamp[20] = +2. * (-std::complex<double> (0, 1) * amp[10] -
      std::complex<double> (0, 1) * amp[12] + std::complex<double> (0, 1) *
      amp[14] + std::complex<double> (0, 1) * amp[13] + std::complex<double>
      (0, 1) * amp[15] + std::complex<double> (0, 1) * amp[16] +
      std::complex<double> (0, 1) * amp[19] + std::complex<double> (0, 1) *
      amp[18] - std::complex<double> (0, 1) * amp[21] + std::complex<double>
      (0, 1) * amp[26] + std::complex<double> (0, 1) * amp[25] -
      std::complex<double> (0, 1) * amp[32] + std::complex<double> (0, 1) *
      amp[35]);
  jamp[21] = +2. * (+std::complex<double> (0, 1) * amp[0] +
      std::complex<double> (0, 1) * amp[1] + std::complex<double> (0, 1) *
      amp[3] - std::complex<double> (0, 1) * amp[16] - std::complex<double> (0,
      1) * amp[17] - std::complex<double> (0, 1) * amp[19] +
      std::complex<double> (0, 1) * amp[21] - std::complex<double> (0, 1) *
      amp[25] + std::complex<double> (0, 1) * amp[29] + std::complex<double>
      (0, 1) * amp[28] + std::complex<double> (0, 1) * amp[32]);
  jamp[22] = +2. * (-std::complex<double> (0, 1) * amp[5] -
      std::complex<double> (0, 1) * amp[7] + std::complex<double> (0, 1) *
      amp[9] + std::complex<double> (0, 1) * amp[8] + std::complex<double> (0,
      1) * amp[15] + std::complex<double> (0, 1) * amp[16] +
      std::complex<double> (0, 1) * amp[19] + std::complex<double> (0, 1) *
      amp[18] - std::complex<double> (0, 1) * amp[20] + std::complex<double>
      (0, 1) * amp[24] + std::complex<double> (0, 1) * amp[23] +
      std::complex<double> (0, 1) * amp[33] - std::complex<double> (0, 1) *
      amp[34]);
  jamp[23] = +2. * (-std::complex<double> (0, 1) * amp[0] -
      std::complex<double> (0, 1) * amp[2] + std::complex<double> (0, 1) *
      amp[4] - std::complex<double> (0, 1) * amp[15] + std::complex<double> (0,
      1) * amp[17] - std::complex<double> (0, 1) * amp[18] +
      std::complex<double> (0, 1) * amp[20] - std::complex<double> (0, 1) *
      amp[23] - std::complex<double> (0, 1) * amp[29] - std::complex<double>
      (0, 1) * amp[28] + std::complex<double> (0, 1) * amp[34]);

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



