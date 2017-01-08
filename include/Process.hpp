#ifndef PROCESS_H_
#define PROCESS_H_

#include <cmath>
#include <complex>
#include "clooptools.h"
#include <string>
#include <array>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <mathematica_wrapper.hpp>
#include <constants.hpp>
#include "LHAPDF/LHAPDF.h"

class Process {
   private:
      // particle masses
      double MasssigmaO, MassphiO, MassGlu, MassTop, MassSq,
         MassSuL, MassSuR, MassSdL, MassSdR, MassSsL, MassSsR,
         MassScL, MassScR, MassSbL, MassSbR, MassStL, MassStR, 
         mu_r, mu_f, dS;
      const LHAPDF::PDF* pdf;
       
      /*
       *    for the moment use those functions for missing ME
       */
      inline double matrix_virt_stub( double, double, double, double, int ) { return 0.; };
      inline double matrix_soft_stub( double, double )  { return 0.; };
      inline double matrix_tree_stub( double, double )  { return 0.; };
      inline double matrix_hard_stub( std::vector< double* >& )  { return 0.; };
      inline double matrix_xsec_stub( double )  { return 0.; };
      
      std::array<double, 2> Pqq( double z) { return std::array<double, 2> {
         CF * (1 + z * z)/(1 - z), - CF * (1 - z) }; 
      }
      std::array<double, 2> Pgq( double z) { return Pqq(1 - z); }; 
      
      std::array<double, 2> Pgg( double z) { return std::array<double, 2> {
         2 * CA * (z/(1 - z) + (1 - z)/z + z * (1 - z)), 0 }; 
      }

      std::array<double, 2> Pqg( double z) { 
         return std::array<double, 2> { 1/2. * ( z*z + (1 - z)*(1 - z) ), -z * (1-z) }; 
      }
      
      // tree-level matrix elements
      double sigmaMSSMTree_uu_suLsuR( double );
      double sigmaMRSSMTree_uubar_suLsuLdagger( double );
      double sigmaMRSSMTree_ddbar_suLsuLdagger( double );
      double sigmaMRSSMTree_gg_suLsuLdagger( double );
      double sigmaMRSSMTree_uu_suLsuR( double );
      double sigmaMRSSMTree_uubar_OO( double );
      double matrixMSSMTree_uu_suLsuR( double, double );
      double matrixMSSMTree_uu_suLsuL( double, double );
      double matrixMSSMTree_ud_suLsdR( double, double );
      double matrixMSSMTree_ud_suLsdL( double, double );

      double matrixMRSSMTree_ddbar_suLsuLdagger( double, double );
      double matrixMRSSMTree_uubar_suLsuLdagger( double, double );
      double matrixMSSMTree_uubar_suLsuRdagger( double, double );
      double matrixMRSSMTree_udbar_suLsdLdagger( double, double );
      double matrixMRSSMTree_GG_suLsuLdagger( double, double );
      
      inline double matrixSgluonTree_qqbar_OO(double);
      inline double matrixSgluonTree_gg_OO(double);

      
      // loop-level matrix elements
      double matrixMSSMVirt_uu_suLsuR(double, double, double, double, int);
      double matrixMRSSMVirt_uu_suLsuR(double, double, double, double, int);
      double matrixMSSMVirt_ud_suLsdR(double, double, double, double, int);
      double matrixMRSSMVirt_ud_suLsdR(double, double, double, double, int);
      double matrixMRSSMVirt_uubar_suLsuLdagger(double, double, double, double, int);
      double matrixMRSSMVirt_ddbar_suLsuLdagger(double, double, double, double, int);
      double matrixMRSSMVirt_GG_suLsuLdagger(double, double, double, double, int);
      
      // soft matrix elements
      double matrixMRSSMSoft_uu_suLsuRg(double, double);
      double matrixMRSSMHard_uu_suLsuRg( std::vector< double* >& );
      //double matrixMRSSMSoft_ud_suLsdRg(double, double);
      
      // pp > suL suL*
      double matrixMRSSMSoft_gg_suLsuLdaggerg( double, double );
      double matrixMRSSMSoft_ddbar_suLsuLdaggerg( double, double );
      double matrixMRSSMSoft_uubar_suLsuLdaggerg( double, double );
      double matrixSimplifiedSoft_uubar_OOg( double, double );
      double matrixSimplifiedHard_uubar_OOg( std::vector< double* >& );
      double matrixSimplifiedSoft_gg_OOg( double, double );
      double matrixSimplifiedHard_gg_OOg( std::vector< double* >& );
      double matrixMRSSMHard_gg_suLsuLdaggerg( std::vector< double* >& );
      double matrixMRSSMHard_uubar_suLsuLdaggerg( std::vector< double* >& );
      double matrixMRSSMHard_ddbar_suLsuLdaggerg( std::vector< double* >& );
      double matrixMRSSMHard_gd_suLsuLdaggerd( std::vector< double* >& );
      double matrixMRSSMHard_gu_suLsuLdaggeru( std::vector< double* >& );
      double matrixMRSSMHard_gu_suLsuRubar( std::vector< double* >& );
      
      // pp > OO
      //double matrixMRSSMSoft_qqbar_OOg( double, double );
      //double matrixMRSSMSoft_gg_OOg( double, double );
      
      double det(std::vector< double* >&, int, int, int, int);
      double determinant( boost::numeric::ublas::matrix<double>& );
      int determinant_sign(const boost::numeric::ublas::permutation_matrix<std::size_t>& );
      
   public: 
      Process(std::string,  boost::property_tree::ptree); 

      double (Process::* matrixelementTree)(double, double); // matrix element squared 
      double (Process::* sigmaPartTree1)(double); // partonic cross section
      double (Process::* sigmaPartTree2)(double);
      std::array<double, 2> (Process::* splitting_kernel1)(double);
      std::array<double, 2> (Process::* splitting_kernel2)(double);
      double (Process::* matrixelementVirt)(double, double, double, 
         double, int); // Re[M^1L M^(B star)]

      double (Process::* matrixelementReal_SC)(double, double);
      double (Process::* matrixelementReal_HC1)(double, double, double);
      double (Process::* matrixelementReal_HC2)(double, double, double);
      double (Process::* matrixelementReal_HnonC)(std::vector< double* >& );
      
      double m1, m2;// masses of final state particle 1 and 2, respectively 
      double k;     // 1/k = average over initial state colors and helicities
      double h;     // h = sum over initial and final state helicities of fermions (_hel = 0 in FormCalc)
      bool partonic;
      std::vector< std::vector<int> > flav;  // each entry is a vector containing the following 3 numbers: 
                                             // first two = initial state flavors, third = how many times does
                                             // initial state occur, e.g. uu is 1, but ud is 2 as there is also du
};
#endif /* PROCESS_H_ */
