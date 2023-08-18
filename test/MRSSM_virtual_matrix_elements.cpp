#include "models/MRSSM.hpp"
#include "constants.hpp"

#include "clooptools.h"

#include "gtest/gtest.h"
#include <iomanip>

namespace {

// Appendix C of arXiv:1707.04557
TEST(MRSSMVirtualMatrixElementTest, AppendixC) {
   MRSSMParameters params;
   params.MassTop = 172.;
   params.MassGlu = 1000.;
   params.MasssigmaO = 5000.;
   params.MassSq = 1500.;
   params.eta_sign = -1;
   params.delta = 0.;
   params.WidthGlu = -1.;

   MRSSM mrssm {params};

   static constexpr double alphas {0.1184};
   static constexpr double S = 6000.*6000;
   static constexpr double T = -22208172;

   static constexpr double mu = 1500.;

   const double born_uu_suLsuR = mrssm.matrixTree_uu_suLsuR(alphas, S, T, 0);
   const double born_uubar_suLsuLdagger = mrssm.matrixTree_uubar_suLsuLdagger(alphas, S, T, 0);
   const double born_ddbar_suLsuLdagger = mrssm.matrixTree_ddbar_suLsuLdagger(alphas, S, T, 0);
   const double born_gg_suLsuLdagger = mrssm.matrixTree_gg_suLsuLdagger(alphas, S, T, 0);
   const double alpha_o_2pi = alphas/(2.*pi);

   ltini();

   double FiniteGs = 0.;
   double Dminus4Coeff = 0.;
   int lambda = -2;

   // double poles
   EXPECT_DOUBLE_EQ(mrssm.matrixVirt_uu_suLsuR(alphas, S, T, FiniteGs, Dminus4Coeff, lambda, mu)/(alpha_o_2pi*born_uu_suLsuR),                   -8/3.);
   EXPECT_DOUBLE_EQ(mrssm.matrixVirt_uubar_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4Coeff, lambda, mu)/(alpha_o_2pi*born_uubar_suLsuLdagger), -8/3.);
   EXPECT_DOUBLE_EQ(mrssm.matrixVirt_ddbar_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4Coeff, lambda, mu)/(alpha_o_2pi*born_ddbar_suLsuLdagger), -8/3.);
   EXPECT_NEAR(mrssm.matrixVirt_gg_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4Coeff, lambda, mu)/(alpha_o_2pi*born_gg_suLsuLdagger),            -6.0, 8e-15);

   // single pole
   EXPECT_DOUBLE_EQ(
      (mrssm.matrixVirt_uu_suLsuR(alphas, S, T, FiniteGs, 0, -1, mu) - 2.*mrssm.matrixVirt_uu_suLsuR(alphas, S, T, FiniteGs, 1, -2, mu))/(alpha_o_2pi*born_uu_suLsuR),
      6.3424494456737985
   );
   EXPECT_DOUBLE_EQ(
      (mrssm.matrixVirt_uubar_suLsuLdagger(alphas, S, T, FiniteGs, 0, -1, mu) - 2.*mrssm.matrixVirt_uubar_suLsuLdagger(alphas, S, T, FiniteGs, 1, -2, mu))/(alpha_o_2pi*born_uubar_suLsuLdagger),
       4.8560767580444564
   );
   EXPECT_DOUBLE_EQ(
      (mrssm.matrixVirt_ddbar_suLsuLdagger(alphas, S, T, FiniteGs, 0, -1, mu) -2.*mrssm.matrixVirt_ddbar_suLsuLdagger(alphas, S, T, FiniteGs, 1, -2, mu))/(alpha_o_2pi*born_ddbar_suLsuLdagger),
      4.8369561733416342
   );
   EXPECT_NEAR(
      (mrssm.matrixVirt_gg_suLsuLdagger(alphas, S, T, FiniteGs, 0, -1, mu) -2.*mrssm.matrixVirt_gg_suLsuLdagger(alphas, S, T, FiniteGs, 1, -2, mu))/(alpha_o_2pi*born_gg_suLsuLdagger),
      8.4161500386713985,
      2e-14
   );

   // finite
   EXPECT_NEAR(
      (
         mrssm.matrixVirt_uu_suLsuR(alphas, S, T, 1, 0, 0, mu)
         - 2.*mrssm.matrixVirt_uu_suLsuR(alphas, S, T, 0, 1, -1, mu)
         + 4.*mrssm.matrixVirt_uu_suLsuR(alphas, S, T, 0, 2, -2, mu)
      )/(alpha_o_2pi*born_uu_suLsuR),
      36.720472180005572,
      5e-14
   );
   EXPECT_DOUBLE_EQ(
      (
         mrssm.matrixVirt_uubar_suLsuLdagger(alphas, S, T, 1, 0, 0, mu)
         - 2.*mrssm.matrixVirt_uubar_suLsuLdagger(alphas, S, T, 0, 1, -1, mu)
         + 4.*mrssm.matrixVirt_uubar_suLsuLdagger(alphas, S, T, 0, 2, -2, mu)
      )/(alpha_o_2pi*born_uubar_suLsuLdagger),
      18.078757030780139
   );
   EXPECT_DOUBLE_EQ(
      (
         mrssm.matrixVirt_ddbar_suLsuLdagger(alphas, S, T, 1, 0, 0, mu)
         -2.*mrssm.matrixVirt_ddbar_suLsuLdagger(alphas, S, T, 0, 1, -1, mu)
         +4.*mrssm.matrixVirt_ddbar_suLsuLdagger(alphas, S, T, 0, 2, -2, mu)
      )/(alpha_o_2pi*born_ddbar_suLsuLdagger),
      -4.713602224435669
   );
   EXPECT_NEAR(
      (
         mrssm.matrixVirt_gg_suLsuLdagger(alphas, S, T, 1, 0, 0, mu)
         -2.*mrssm.matrixVirt_gg_suLsuLdagger(alphas, S, T, 0, 1, -1, mu)
         +4.*mrssm.matrixVirt_gg_suLsuLdagger(alphas, S, T, 0, 2, -2, mu)
      )/(alpha_o_2pi*born_gg_suLsuLdagger),
      -3.7965233237193581,
      3e-14
   );
}
}
