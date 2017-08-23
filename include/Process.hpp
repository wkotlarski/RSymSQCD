#ifndef PROCESS_H_
#define PROCESS_H_

#include <cmath>
#include <complex>
#include "clooptools.h"
#include <string>
#include <array>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <mathematica_wrapper.hpp>
#include <constants.hpp>

#include "MRSSM.h"

class Process : public MRSSM{

protected:
    // particle masses
   double MassSq;
   const double mu_r;
   const double mu_f;
   const double dS;
   const double WidthGlu;
   const double delta;
   const double eta_sign;

private:
      // gauge vector for DR ME
      std::array<double,4> eta; 

      /*
       *    for the moment use those functions for missing ME
       */
      inline double matrix_virt_stub( double, double, double, double, int ) { return 0.; };
      inline double matrix_soft_stub( double, double )  { return 0.; };
      inline double matrix_tree_stub( double, double ) const  { return 0.; };
      inline double matrix_hard_stub( std::vector< double* >& )  { return 0.; };
      inline double matrix_xsec_stub( double )  { return 0.; };
      
      std::array<double, 2> Pqq( double z) const noexcept { 
         return std::array<double, 2> { CF * (1 + z * z)/(1 - z), - CF * (1 - z) }; }

      std::array<double, 2> Pgq( double z) const noexcept { return Pqq(1 - z); }; 
      
      std::array<double, 2> Pgg( double z) const noexcept { return std::array<double, 2> {
         2 * CA * (z/(1 - z) + (1 - z)/z + z * (1 - z)), 0 }; 
      }

      std::array<double, 2> Pqg( double z) const noexcept { 
         return std::array<double, 2> { 1/2. * ( z*z + (1 - z)*(1 - z) ), -z * (1-z) }; 
      }
      
      // tree-level xsections
      double sigmaMSSMTree_uu_suLsuR( double );
      double sigmaMRSSMTree_uubar_suLsuLdagger( double );
      double sigmaMRSSMTree_ddbar_suLsuLdagger( double );
      double sigmaMRSSMTree_gg_suLsuLdagger( double );
      double sigmaMRSSMTree_uu_suLsuR( double );
      double sigmaMRSSMTree_uubar_OO( double );

      // tree-level matrix elements
      double matrixMSSMTree_uu_suLsuR( double, double ) const;
      double matrixMSSMTree_uu_suLsuL( double, double ) const;
      double matrixMSSMTree_ud_suLsdR( double, double ) const;
      double matrixMSSMTree_ud_suLsdL( double, double ) const;
      double matrixMRSSMTree_ddbar_suLsuLdagger( double, double ) const;
      double matrixMRSSMTree_uubar_suLsuLdagger( double, double ) const;
      double matrixMSSMTree_uubar_suLsuRdagger( double, double ) const;
      double matrixMRSSMTree_udbar_suLsdLdagger( double, double ) const;
      double matrixMRSSMTree_GG_suLsuLdagger( double, double ) const;
      double matrixSgluonTree_qqbar_OO(double) const;
      double matrixSgluonTree_gg_OO(double) const;
      
      // loop-level matrix elements
      double matrixMSSMVirt_uu_suLsuR(double, double, double, double, int);
      double matrixMRSSMVirt_uu_suLsuR(double, double, double, double, int);
      double matrixMSSMVirt_ud_suLsdR(double, double, double, double, int);
      double matrixMRSSMVirt_ud_suLsdR(double, double, double, double, int);
      double matrixMRSSMVirt_uubar_suLsuLdagger(double, double, double, double, int);
      double matrixMRSSMVirt_ddbar_suLsuLdagger(double, double, double, double, int);
      double matrixMRSSMVirt_GG_suLsuLdagger(double, double, double, double, int);
      // @todo: move to separate class
      double matrixSMVirt_eebar_ttbar (double, double, double, double, int);
      double matrixSMHard_eebar_ttbarg (const std::vector< double* >&) const;
      double matrixSMTree_eebar_ttbar( double, double ) const;

      // soft matrix elements
      double matrixMRSSMSoft_uu_suLsuRg(double, double);
      double matrixMRSSMHard_uu_suLsuRg( const std::vector< double* >& ) const ;
      
      // pp > suL suL*
      double matrixMRSSMSoft_gg_suLsuLdaggerg( double, double );
      double matrixMRSSMSoft_ddbar_suLsuLdaggerg( double, double );
      double matrixMRSSMSoft_uubar_suLsuLdaggerg( double, double );
      double matrixSimplifiedSoft_uubar_OOg( double, double );
      double matrixSimplifiedHard_uubar_OOg( const std::vector< double* >& ) const ;
      double matrixSimplifiedSoft_gg_OOg( double, double );
      double matrixSimplifiedHard_gg_OOg( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gg_suLsuLdaggerg( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_uubar_suLsuLdaggerg( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_ddbar_suLsuLdaggerg( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gd_suLsuLdaggerd( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuLdaggeru_DR( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuLdaggeru_DR2( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuLdaggeru_DR_wEta( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuLdaggeru_DS( const std::vector< double* >& ) const;
      double matrixMRSSMHard_gu_suLsuLdaggeru_DS_CSub2( const std::vector< double* >& ) const noexcept;
      double matrixMRSSMHard_gubar_suLsuLdaggerubar_DS( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gubar_suLsuLdaggerubar_DS_CSub1( const std::vector< double* >& ) const noexcept;
      double matrixMRSSMHard_gu_suLsuLdaggeru( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuRubar( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuRubar_DR( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuRubar_DR_wEtaDep( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuRubar_DS( const std::vector< double* >& ) const ;
      double matrixMRSSMHard_gu_suLsuRubar_DS_CSub1( const std::vector< double* >& ) const noexcept;
      double matrixMRSSMHard_gu_suLsuRubar_DS_CSub2( const std::vector< double* >& ) const noexcept;
      
   public: 
      Process(std::string,  boost::property_tree::ptree);

      double (Process::* matrixelementTree)(double, double) const; // matrix element squared 
      double (Process::* sigmaPartTree1) (double); // partonic cross section
      double (Process::* sigmaPartTree2) (double);
      std::array<double, 2> (Process::* splitting_kernel1) (double) const noexcept;
      std::array<double, 2> (Process::* splitting_kernel2) (double) const noexcept;
      double (Process::* matrixelementVirt)(double, double, double, 
         double, int); // Re[M^1L M^(B star)]

      double (Process::* matrixelementReal_SC)(double, double);
      double (Process::* matrixelementReal_HC1)(double, double, double);
      double (Process::* matrixelementReal_HC2)(double, double, double);
      double (Process::* matrixelementReal_HnonC)(const std::vector< double* >& ) const ;
      double (Process::* matrixelementReal_HnonC_CSub1)(const std::vector< double* >& ) const noexcept;
      double (Process::* matrixelementReal_HnonC_CSub2)(const std::vector< double* >& ) const noexcept;
      
      double m1, m2;// masses of final state particle 1 and 2, respectively 
      double c1, c2, c3, c4;
      double k;     // 1/k = average over initial state colors and helicities
      double h;     // h = sum over initial and final state helicities of fermions (_hel = 0 in FormCalc)
      bool partonic;
      std::vector< std::vector<int> > flav;  // each entry is a vector containing the following 3 numbers: 
                                             // first two = initial state flavors, third = how many times does
                                             // initial state occur, e.g. uu is 1, but ud is 2 as there is also du
};
#endif /* PROCESS_H_ */
