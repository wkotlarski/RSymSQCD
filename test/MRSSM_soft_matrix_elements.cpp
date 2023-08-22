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
   const double theta = std::acos((S+2*T-2*params.MassSq*params.MassSq)/(std::sqrt(1-4*params.MassSq*params.MassSq/S)*S));

   static constexpr double mu = 1500.;
   static constexpr double dS = 5e-4;

   EXPECT_NEAR(
      mrssm.matrixSoft_uu_suLsuRg(alphas, S, theta, dS, mu),
      7.4718949947862054e-10,
      3e-25.
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_uubar_suLsuLdaggerg(alphas, S, theta, dS, mu),
      3.1357342708229199e-10,
      3e-25
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_ddbar_suLsuLdaggerg(alphas, S, theta, dS, mu),
      1.1526982445319974e-10,
      2e-26.
   );
   EXPECT_NEAR(
      mrssm.matrixSoft_gg_suLsuLdaggerg(alphas, S, theta, dS, mu),
      2.1425060919474067e-10,
      3e-25
   );
}
}
