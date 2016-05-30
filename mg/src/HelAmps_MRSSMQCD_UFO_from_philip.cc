//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.3.3, 2015-10-25
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include "HelAmps_MRSSMQCD_UFO_from_philip.h"
#include <complex> 
#include <cmath> 
#include <iostream> 
#include <cstdlib> 
using namespace std; 

namespace MG5_MRSSMQCD_UFO_from_philip 
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}

void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], pow((pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2)), 0.5)); 
    if (pp == 0.0)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.0), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * pow(2.0 * p[0], 0.5), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}

void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
    if (pp == 0.000)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[abs(ip)]; 
      fo[3] = ip * nsf * sqm[abs(ip)]; 
      fo[4] = im * nsf * sqm[abs(im)]; 
      fo[5] = ip * sqm[abs(im)]; 
    }
    else
    {
      pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.00), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * pow(2.0 * p[0], 0.5); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}


void txxxxx(double p[4], double tmass, int nhel, int nst, complex<double>
    tc[18])
{
  complex<double> ft[6][4], ep[4], em[4], e0[4]; 
  double pt, pt2, pp, pzpt, emp, sqh, sqs; 
  int i, j; 

  sqh = pow(0.5, 0.5); 
  sqs = pow(0.5/3, 0.5); 

  pt2 = p[1] * p[1] + p[2] * p[2]; 
  pp = min(p[0], pow(pt2 + p[3] * p[3], 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 

  ft[4][0] = complex<double> (p[0] * nst, p[3] * nst); 
  ft[5][0] = complex<double> (p[1] * nst, p[2] * nst); 

  // construct eps+
  if(nhel >= 0)
  {
    if(pp == 0)
    {
      ep[0] = complex<double> (0, 0); 
      ep[1] = complex<double> (-sqh, 0); 
      ep[2] = complex<double> (0, nst * sqh); 
      ep[3] = complex<double> (0, 0); 
    }
    else
    {
      ep[0] = complex<double> (0, 0); 
      ep[3] = complex<double> (pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = p[3]/(pp * pt) * sqh; 
        ep[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        ep[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        ep[1] = complex<double> (-sqh, 0); 
        ep[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }

  }

  // construct eps-
  if(nhel <= 0)
  {
    if(pp == 0)
    {
      em[0] = complex<double> (0, 0); 
      em[1] = complex<double> (sqh, 0); 
      em[2] = complex<double> (0, nst * sqh); 
      em[3] = complex<double> (0, 0); 
    }
    else
    {
      em[0] = complex<double> (0, 0); 
      em[3] = complex<double> (-pt/pp * sqh, 0); 

      if(pt != 0)
      {
        pzpt = -p[3]/(pp * pt) * sqh; 
        em[1] = complex<double> (-p[1] * pzpt, -nst * p[2]/pt * sqh); 
        em[2] = complex<double> (-p[2] * pzpt, nst * p[1]/pt * sqh); 
      }
      else
      {
        em[1] = complex<double> (sqh, 0); 
        em[2] = complex<double> (0, nst * Sgn(sqh, p[3])); 
      }
    }
  }

  // construct eps0
  if(fabs(nhel) <= 1)
  {
    if(pp == 0)
    {
      e0[0] = complex<double> (0, 0); 
      e0[1] = complex<double> (0, 0); 
      e0[2] = complex<double> (0, 0); 
      e0[3] = complex<double> (1, 0); 
    }
    else
    {
      emp = p[0]/(tmass * pp); 
      e0[0] = complex<double> (pp/tmass, 0); 
      e0[3] = complex<double> (p[3] * emp, 0); 

      if(pt != 0)
      {
        e0[1] = complex<double> (p[1] * emp, 0); 
        e0[2] = complex<double> (p[2] * emp, 0); 
      }
      else
      {
        e0[1] = complex<double> (0, 0); 
        e0[2] = complex<double> (0, 0); 
      }
    }
  }

  if(nhel == 2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = ep[i] * ep[j]; 
    }
  }
  else if(nhel == -2)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = em[i] * em[j]; 
    }
  }
  else if(tmass == 0)
  {
    for(j = 0; j < 4; j++ )
    {
      for(i = 0; i < 4; i++ )
        ft[i][j] = 0; 
    }
  }
  else if(tmass != 0)
  {
    if(nhel == 1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (ep[i] * e0[j] + e0[i] * ep[j]); 
      }
    }
    else if(nhel == 0)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqs * (ep[i] * em[j] + em[i] * ep[j]
         + 2.0 * e0[i] * e0[j]); 
      }
    }
    else if(nhel == -1)
    {
      for(j = 0; j < 4; j++ )
      {
        for(i = 0; i < 4; i++ )
          ft[i][j] = sqh * (em[i] * e0[j] + e0[i] * em[j]); 
      }
    }
    else
    {
      std::cerr <<  "Invalid helicity in txxxxx.\n"; 
      std::exit(1); 
    }
  }

  tc[0] = ft[4][0]; 
  tc[1] = ft[5][0]; 

  for(j = 0; j < 4; j++ )
  {
    for(i = 0; i < 4; i++ )
      tc[j * 4 + i + 2] = ft[j][i]; 
  }
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = pow(0.5, 0.5); 
  hel = double(nhel); 
  nsvahl = nsv * abs(hel); 
  pt2 = pow(p[1], 2) + pow(p[2], 2); 
  pp = min(p[0], pow(pt2 + pow(p[3], 2), 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = pow(pow(p[1], 2) + pow(p[2], 2), 0.5); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void FFV1P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5]
      * F2[3]);
  V3[3] = denom * - cI * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3]
      * F2[4]);
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] *
      F2[4] + F1[4] * F2[3]));
  V3[5] = denom * - cI * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5]
      * F2[3]);
}


void VSS1_0(complex<double> V1[], complex<double> S2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> TMP8; 
  P2[0] = S2[0].real(); 
  P2[1] = S2[1].real(); 
  P2[2] = S2[1].imag(); 
  P2[3] = S2[0].imag(); 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP8 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  vertex = COUP * S2[2] * S3[2] * (-cI * (TMP7) + cI * (TMP8)); 
}


void FFV1_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * (+cI *
      (V3[4]) - V3[3]))));
  F2[3] = denom * - cI * (F1[2] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] * (V3[2] +
      V3[5]))));
  F2[4] = denom * - cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] +
      cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * - 1. *
      (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] *
      (V3[2] - V3[5]))));
}


void VSS1_2(complex<double> V1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> S2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> denom; 
  complex<double> TMP8; 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  S2[0] = +V1[0] + S3[0]; 
  S2[1] = +V1[1] + S3[1]; 
  P2[0] = -S2[0].real(); 
  P2[1] = -S2[1].real(); 
  P2[2] = -S2[1].imag(); 
  P2[3] = -S2[0].imag(); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP8 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  S2[2] = denom * S3[2] * (-cI * (TMP8) + cI * (TMP7)); 
}


void VVV1_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> TMP5; 
  complex<double> TMP4; 
  complex<double> TMP3; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP8 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP5 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP4 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP7 = (V1[2] * P2[0] - V1[3] * P2[1] - V1[4] * P2[2] - V1[5] * P2[3]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP1 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP0 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP3 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP2 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  vertex = COUP * (TMP1 * (-cI * (TMP0) + cI * (TMP2)) + (TMP3 * (-cI * (TMP5)
      + cI * (TMP4)) + TMP6 * (-cI * (TMP7) + cI * (TMP8))));
}


void FFV1C1_1(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * - cI * (F1[2] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] -
      cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[3] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] * - 1. *
      (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + M2 * (F1[4] * (V3[2] + V3[5]) + F1[5] *
      (V3[3] + cI * (V3[4])))));
  F2[3] = denom * cI * (F1[2] * (P2[0] * (+cI * (V3[4]) - V3[3]) + (P2[1] *
      (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) + P2[3] *
      (V3[3] - cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      - 1. * (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[4] * (+cI * (V3[4]) - V3[3]) + F1[5] *
      (V3[5] - V3[2]))));
  F2[4] = denom * cI * (F1[4] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * (+cI *
      (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) - P2[3] *
      (V3[2] + V3[5])))) + (F1[5] * (P2[0] * (V3[3] + cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[2]) + cI * (V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + M2 * (F1[2] * (V3[5] - V3[2]) + F1[3] *
      (V3[3] + cI * (V3[4])))));
  F2[5] = denom * - cI * (F1[4] * (P2[0] * (+cI * (V3[4]) - V3[3]) + (P2[1] *
      (V3[2] + V3[5]) + (P2[2] * - 1. * (+cI * (V3[2] + V3[5])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + (F1[5] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] + cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[2] * (+cI * (V3[4]) - V3[3]) + F1[3] * (V3[2] +
      V3[5]))));
}


void VVVV4P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP10; 
  complex<double> denom; 
  complex<double> TMP9; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP9 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  TMP10 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V3[2] * TMP9) + cI * (V2[2] * TMP10)); 
  V1[3] = denom * (-cI * (V3[3] * TMP9) + cI * (V2[3] * TMP10)); 
  V1[4] = denom * (-cI * (V3[4] * TMP9) + cI * (V2[4] * TMP10)); 
  V1[5] = denom * (-cI * (V3[5] * TMP9) + cI * (V2[5] * TMP10)); 
}


void FFS3C1_1(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * F1[2] * S3[2] * M2; 
  F2[3] = denom * cI * F1[3] * S3[2] * M2; 
  F2[4] = denom * cI * S3[2] * (F1[2] * (P2[3] - P2[0]) + F1[3] * (P2[1] + cI *
      (P2[2])));
  F2[5] = denom * - cI * S3[2] * (F1[2] * (+cI * (P2[2]) - P2[1]) + F1[3] *
      (P2[0] + P2[3]));
}


void FFS5_0(complex<double> F1[], complex<double> F2[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP14; 
  TMP14 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  vertex = COUP * - cI * TMP14 * S3[2]; 
}


void VVSS1_3(complex<double> V1[], complex<double> V2[], complex<double> S4[],
    complex<double> COUP, double M3, double W3, complex<double> S3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  double P3[4]; 
  complex<double> denom; 
  S3[0] = +V1[0] + V2[0] + S4[0]; 
  S3[1] = +V1[1] + V2[1] + S4[1]; 
  P3[0] = -S3[0].real(); 
  P3[1] = -S3[1].real(); 
  P3[2] = -S3[1].imag(); 
  P3[3] = -S3[0].imag(); 
  TMP1 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  S3[2] = denom * cI * TMP1 * S4[2]; 
}


void FFS5C1_1(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * - cI * S3[2] * (F1[4] * (P2[0] + P2[3]) + F1[5] * (P2[1] + cI
      * (P2[2])));
  F2[3] = denom * cI * S3[2] * (F1[4] * (+cI * (P2[2]) - P2[1]) + F1[5] *
      (P2[3] - P2[0]));
  F2[4] = denom * cI * F1[4] * S3[2] * M2; 
  F2[5] = denom * cI * F1[5] * S3[2] * M2; 
}


void VVVV3P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP10; 
  complex<double> TMP6; 
  complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP10 = (V3[2] * V4[2] - V3[3] * V4[3] - V3[4] * V4[4] - V3[5] * V4[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (TMP6 * V4[2]) + cI * (V2[2] * TMP10)); 
  V1[3] = denom * (-cI * (TMP6 * V4[3]) + cI * (V2[3] * TMP10)); 
  V1[4] = denom * (-cI * (TMP6 * V4[4]) + cI * (V2[4] * TMP10)); 
  V1[5] = denom * (-cI * (TMP6 * V4[5]) + cI * (V2[5] * TMP10)); 
}


void VSS1P0_1(complex<double> S2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> denom; 
  P2[0] = S2[0].real(); 
  P2[1] = S2[1].real(); 
  P2[2] = S2[1].imag(); 
  P2[3] = S2[0].imag(); 
  P3[0] = S3[0].real(); 
  P3[1] = S3[1].real(); 
  P3[2] = S3[1].imag(); 
  P3[3] = S3[0].imag(); 
  V1[0] = +S2[0] + S3[0]; 
  V1[1] = +S2[1] + S3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * S2[2] * S3[2] * (-cI * (P2[0]) + cI * (P3[0])); 
  V1[3] = denom * S2[2] * S3[2] * (-cI * (P2[1]) + cI * (P3[1])); 
  V1[4] = denom * S2[2] * S3[2] * (-cI * (P2[2]) + cI * (P3[2])); 
  V1[5] = denom * S2[2] * S3[2] * (-cI * (P2[3]) + cI * (P3[3])); 
}


void FFV1_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP11; 
  TMP11 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  vertex = COUP * - cI * TMP11; 
}


void FFS5_2(complex<double> F1[], complex<double> S3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + S3[0]; 
  F2[1] = +F1[1] + S3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * S3[2] * (F1[4] * (P2[0] - P2[3]) + F1[5] * (+cI *
      (P2[2]) - P2[1]));
  F2[3] = denom * - cI * S3[2] * (F1[4] * (P2[1] + cI * (P2[2])) - F1[5] *
      (P2[0] + P2[3]));
  F2[4] = denom * cI * F1[4] * S3[2] * M2; 
  F2[5] = denom * cI * F1[5] * S3[2] * M2; 
}


void FFV1C1_0(complex<double> F2[], complex<double> F1[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP12; 
  TMP12 = -1. * (F1[2] * (F2[4] * (V3[2] - V3[5]) + F2[5] * (+cI * (V3[4]) -
      V3[3])) + (F1[3] * (F2[4] * - 1. * (V3[3] + cI * (V3[4])) + F2[5] *
      (V3[2] + V3[5])) + (F1[4] * (F2[2] * (V3[2] + V3[5]) + F2[3] * (V3[3] -
      cI * (V3[4]))) + F1[5] * (F2[2] * (V3[3] + cI * (V3[4])) + F2[3] * (V3[2]
      - V3[5])))));
  vertex = COUP * - cI * TMP12; 
}


void FFS5C1_0(complex<double> F2[], complex<double> F1[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP14; 
  TMP14 = (F1[4] * F2[4] + F1[5] * F2[5]); 
  vertex = COUP * - cI * TMP14 * S3[2]; 
}


void FFS3C1_2(complex<double> F2[], complex<double> S3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + S3[0]; 
  F1[1] = +F2[1] + S3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * F2[2] * S3[2] * M1; 
  F1[3] = denom * cI * F2[3] * S3[2] * M1; 
  F1[4] = denom * - cI * S3[2] * (F2[2] * - 1. * (P1[0] + P1[3]) + F2[3] * (+cI
      * (P1[2]) - P1[1]));
  F1[5] = denom * cI * S3[2] * (F2[2] * (P1[1] + cI * (P1[2])) + F2[3] * (P1[0]
      - P1[3]));
}


void FFV1_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI
      * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * - 1. * (V3[2] +
      V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] + cI *
      (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])))));
  F1[3] = denom * - cI * (F2[2] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] *
      - 1. * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] *
      (V3[5] - V3[2]))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (V3[5] - V3[2]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] +
      V3[5]))));
}


void VVSS1_0(complex<double> V1[], complex<double> V2[], complex<double> S3[],
    complex<double> S4[], complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  TMP1 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  vertex = COUP * - cI * TMP1 * S4[2] * S3[2]; 
}


void FFS3C1_0(complex<double> F2[], complex<double> F1[], complex<double> S3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP13; 
  TMP13 = (F1[2] * F2[2] + F1[3] * F2[3]); 
  vertex = COUP * - cI * TMP13 * S3[2]; 
}


void VVVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP6; 
  complex<double> denom; 
  complex<double> TMP9; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP9 = (V2[2] * V4[2] - V2[3] * V4[3] - V2[4] * V4[4] - V2[5] * V4[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (TMP6 * V4[2]) + cI * (V3[2] * TMP9)); 
  V1[3] = denom * (-cI * (TMP6 * V4[3]) + cI * (V3[3] * TMP9)); 
  V1[4] = denom * (-cI * (TMP6 * V4[4]) + cI * (V3[4] * TMP9)); 
  V1[5] = denom * (-cI * (TMP6 * V4[5]) + cI * (V3[5] * TMP9)); 
}


void VVSS1P0_1(complex<double> V2[], complex<double> S3[], complex<double>
    S4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  V1[0] = +V2[0] + S3[0] + S4[0]; 
  V1[1] = +V2[1] + S3[1] + S4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * - cI * V2[2] * S4[2] * S3[2]; 
  V1[3] = denom * - cI * V2[3] * S4[2] * S3[2]; 
  V1[4] = denom * - cI * V2[4] * S4[2] * S3[2]; 
  V1[5] = denom * - cI * V2[5] * S4[2] * S3[2]; 
}


void VVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  double P1[4]; 
  complex<double> TMP0; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> TMP5; 
  complex<double> TMP4; 
  complex<double> denom; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP5 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP4 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP6 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP0 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP2 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP6 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI *
      (TMP0) + cI * (TMP2)) + V3[2] * (-cI * (TMP5) + cI * (TMP4))));
  V1[3] = denom * (TMP6 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI *
      (TMP0) + cI * (TMP2)) + V3[3] * (-cI * (TMP5) + cI * (TMP4))));
  V1[4] = denom * (TMP6 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI *
      (TMP0) + cI * (TMP2)) + V3[4] * (-cI * (TMP5) + cI * (TMP4))));
  V1[5] = denom * (TMP6 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI *
      (TMP0) + cI * (TMP2)) + V3[5] * (-cI * (TMP5) + cI * (TMP4))));
}


}  // end namespace $(namespace)s_MRSSMQCD_UFO_from_philip

