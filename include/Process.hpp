#ifndef PROCESS_H_
#define PROCESS_H_

#include "constants.hpp"
#include "mathematica_wrapper.hpp"

#include <boost/property_tree/ptree.hpp>

#include <array>
#include <string>
#include <vector>

class Process {
   public:
      // gauge vector for DR ME
      std::array<double, 4> eta;

      // particle masses
      double MasssigmaO, MassphiO, MassGlu, MassTop, MassSq, MassSquarks, WidthGlu;

      /*
       *    for the moment use those functions for missing ME
       */
      inline double matrix_virt_stub(double, double, double, double, double, int, double) { return 0.; };
      inline double matrix_soft_stub(double, double, double, double, double)  { return 0.; };
      inline double matrix_tree_stub(double, double, double) const { return 0.; };
      inline double matrix_hard_stub( std::vector< double* >& )  { return 0.; };
      inline double matrix_xsec_stub(double, double )  { return 0.; };

      // tree-level matrix elements
      double sigmaMSSMTree_uu_suLsuR(double, double);
      double sigmaMRSSMTree_uubar_suLsuLdagger(double, double);
      double sigmaMRSSMTree_ddbar_suLsuLdagger(double, double);
      double sigmaMRSSMTree_gg_suLsuLdagger(double, double );
      double sigmaMRSSMTree_uu_suLsuR(double, double );
      double sigmaMRSSMTree_uubar_OO(double, double );
      double matrixMSSMTree_uu_suLsuR(double, double, double) const;
      double matrixMSSMTree_uu_suLsuL(double, double, double) const;
      double matrixMSSMTree_ud_suLsdR(double, double, double) const;
      double matrixMSSMTree_ud_suLsdL(double, double, double) const;

      double matrixMRSSMTree_ddbar_suLsuLdagger(double, double, double) const;
      double matrixMRSSMTree_uubar_suLsuLdagger(double, double, double) const;
      double matrixMSSMTree_uubar_suLsuRdagger(double, double, double) const;
      double matrixMRSSMTree_udbar_suLsdLdagger(double, double, double) const;
      double matrixMRSSMTree_GG_suLsuLdagger(double, double, double) const;

      inline double matrixSgluonTree_qqbar_OO(double, double);
      inline double matrixSgluonTree_gg_OO(double, double);


      // loop-level matrix elements
      double matrixMSSMVirt_uu_suLsuR(double, double, double, double, double, int, double);
      double matrixMRSSMVirt_uu_suLsuR(double, double, double, double, double, int, double);
      double matrixMSSMVirt_ud_suLsdR(double, double, double, double, double, int, double);
      double matrixMRSSMVirt_ud_suLsdR(double, double, double, double, double, int, double);
      double matrixMRSSMVirt_uubar_suLsuLdagger(double, double, double, double, double, int, double);
      double matrixMRSSMVirt_ddbar_suLsuLdagger(double, double, double, double, double, int, double);
      double matrixMRSSMVirt_GG_suLsuLdagger(double, double, double, double, double, int, double);

      // soft matrix elements
      double matrixMRSSMSoft_uu_suLsuRg(double, double, double, double, double);
      double matrixMRSSMHard_uu_suLsuRg(double, std::array<std::array<double, 4>, 5> const&) const;
      //double matrixMRSSMSoft_ud_suLsdRg(double, double);

      // pp > suL suL*
      double matrixMRSSMSoft_gg_suLsuLdaggerg(double, double, double, double, double);
      double matrixMRSSMSoft_ddbar_suLsuLdaggerg(double, double, double, double, double);
      double matrixMRSSMSoft_uubar_suLsuLdaggerg(double, double, double, double, double);
      // double matrixSimplifiedSoft_uubar_OOg( double, double );
      // double matrixSimplifiedHard_uubar_OOg(std::array<std::array<double, 4>, 5> const&) const;
      // double matrixSimplifiedSoft_gg_OOg( double, double );
      // double matrixSimplifiedHard_gg_OOg(std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gg_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_uubar_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_ddbar_suLsuLdaggerg(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gd_suLsuLdaggerd(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gu_suLsuLdaggeru_DR(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gu_suLsuLdaggeru_DR_wEta(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gu_suLsuLdaggeru_DS(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gu_suLsuLdaggeru(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gu_suLsuRubar(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gu_suLsuRubar_DR(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gu_suLsuRubar_DR_wEtaDep(double, std::array<std::array<double, 4>, 5> const&) const;
      double matrixMRSSMHard_gu_suLsuRubar_DS(double, std::array<std::array<double, 4>, 5> const&) const;

      // pp > OO
      //double matrixMRSSMSoft_qqbar_OOg( double, double );
      //double matrixMRSSMSoft_gg_OOg( double, double );

      Process() = default;
      Process(std::string const&, boost::property_tree::ptree const&);

      double (Process::* matrixelementTree)(double, double, double) const; // matrix element squared
      double (Process::* sigmaPartTree1)(double, double); // partonic cross section
      double (Process::* sigmaPartTree2)(double, double);
      std::array<double, 2> (*splitting_kernel1)(double);
      std::array<double, 2> (*splitting_kernel2)(double);
      double (Process::* matrixelementVirt)(double, double, double, double,
         double, int, double); // Re[M^1L M^(B star)]

      double (Process::* matrixelementReal_SC)(double, double, double, double, double);
      double (Process::* matrixelementReal_HC1)(double, double, double);
      double (Process::* matrixelementReal_HC2)(double, double, double);
      double (Process::* matrixelementReal_HnonC)(std::array<std::array<double, 4>, 5> const&) const;

      double m1, m2;// masses of final state particle 1 and 2, respectively
      // double c1, c2, c3, c4;
      double k;     // 1/k = average over initial state colors and helicities
      double h;     // h = sum over initial and final state helicities of fermions (_hel = 0 in FormCalc)
      bool partonic;
};
#endif /* PROCESS_H_ */
