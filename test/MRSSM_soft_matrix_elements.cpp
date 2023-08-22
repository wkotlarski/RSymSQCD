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
   static constexpr double alphas_by_2Pi = alphas/(2*pi);
   static constexpr double S = 6000.*6000;
   static constexpr double T = -22208172;
   const double theta = std::acos((S+2*T-2*params.MassSq*params.MassSq)/(std::sqrt(1-4*params.MassSq*params.MassSq/S)*S));

   static constexpr double mu = 1500.;

   static constexpr double dS = 5e-4;
   const double Asc_qq = 4/3.*(2*std::log(dS) + 1.5);
   static constexpr int nf = 6;
   const double Asc_gg = 2*CA*std::log(dS) + (11*CA-2*(nf-1))/6.;

   // uu -> suLsuR
   const double uu_suLsuR_tree = mrssm.matrixTree_uu_suLsuR(alphas, S, T, 0.);
   EXPECT_NEAR(
      mrssm.matrixSoft_uu_suLsuRg_dp(alphas, S, theta)/(alphas_by_2Pi*uu_suLsuR_tree),
      8/3.,
      9e-16
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_uu_suLsuRg_sp(alphas, S, theta, dS, mu)/(alphas_by_2Pi*uu_suLsuR_tree) + 2.*Asc_qq,
      -6.342449445673779,
      3e-14
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_uu_suLsuRg_finite(alphas, S, theta, dS, mu),
      3.4307758460708908,
      2e-15
   );

   // uubar -> suLsuL*
   const double uubar_suLsuLdagger_tree = mrssm.matrixTree_uubar_suLsuLdagger(alphas, S, T, 0.);
   EXPECT_NEAR(
      mrssm.matrixSoft_uubar_suLsuLdaggerg_dp(alphas, S, theta)/(alphas_by_2Pi*uubar_suLsuLdagger_tree),
      8/3.,
      5e-16
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_uubar_suLsuLdaggerg_sp(alphas, S, theta, dS, mu)/(alphas_by_2Pi*uubar_suLsuLdagger_tree) + 2.*Asc_qq,
      -4.8560767580444697,
      3e-14
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_uubar_suLsuLdaggerg_finite(alphas, S, theta, dS, mu),
      1.4397955811133312,
      7e-16
   );

   // ddbar -> suLsuL*
   const double ddbar_suLsuLdagger_tree = mrssm.matrixTree_ddbar_suLsuLdagger(alphas, S, T, 0.);
   EXPECT_NEAR(
      mrssm.matrixSoft_ddbar_suLsuLdaggerg_dp(alphas, S, theta)/(alphas_by_2Pi*ddbar_suLsuLdagger_tree),
      8/3.,
      5e-16
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_ddbar_suLsuLdaggerg_sp(alphas, S, theta, dS, mu)/(alphas_by_2Pi*ddbar_suLsuLdagger_tree) + 2.*Asc_qq,
      -4.8369561733416369,
      4e-15
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_ddbar_suLsuLdaggerg_finite(alphas, S, theta, dS, mu),
      0.52926992388252259,
      2e-16
   );

   // gg -> suLsuLdagger
   const double gg_suLsuLdagger_tree = mrssm.matrixTree_gg_suLsuLdagger(alphas, S, T, 0.);
   EXPECT_NEAR(
      mrssm.matrixSoft_gg_suLsuLdaggerg_dp(alphas, S, theta)/(alphas_by_2Pi*gg_suLsuLdagger_tree),
      6,
      9e-15
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_gg_suLsuLdaggerg_sp(alphas, S, theta, dS, mu)/(alphas_by_2Pi*gg_suLsuLdagger_tree) + 2.*Asc_gg,
      -8.4161500386714039,
      2e-13
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_gg_suLsuLdaggerg_finite(alphas, S, theta, dS, mu),
      0.98374751725525633,
      3e-16
   );
}
}
