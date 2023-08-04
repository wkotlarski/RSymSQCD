#include "models/MRSSM.hpp"
#include "constants.hpp"

#include "clooptools.h"

#include "gtest/gtest.h"
#include <iomanip>

namespace {

TEST(MRSSMVirtualMatrixElementTest, BMP1) {
   MRSSMParameters params;
   params.MassTop = 172.;
   params.MassGlu = 1000.;
   params.MasssigmaO = 5000.;
   params.MassSq = 1500.;
   params.eta_sign = -1;
   params.delta = 0.;
   params.WidthGlu = -1.;

   MRSSM mrssm {params};

   static constexpr double alphas {8.41234775121963707e-02};
   static constexpr double S = 4000.*4000;
   static constexpr double T = -2980004;

   static constexpr double mu = 500.;

   ltini();

   double FiniteGs = 1.;
   double Dminus4 = 0.;
   int divergence = 0.;

   // finite
   EXPECT_NEAR(mrssm.matrixVirt_uu_suLsuR(alphas, S, T, FiniteGs, Dminus4, divergence, mu), 0.062335750753775121, 0.);
   EXPECT_NEAR(mrssm.matrixVirt_uubar_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4, divergence, mu), 0.048420797803705902, 0.);
   EXPECT_NEAR(mrssm.matrixVirt_ddbar_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4, divergence, mu), -0.0046349281931064362, 0.);
   EXPECT_NEAR(mrssm.matrixVirt_GG_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4, divergence, mu), -0.010870674966149739, 0.);

   FiniteGs = 0.;
   Dminus4 = -2.;
   divergence = -1;
   EXPECT_NEAR(mrssm.matrixVirt_uu_suLsuR(alphas, S, T, FiniteGs, Dminus4, divergence, mu), 0.012091082721229563, 0.);
   EXPECT_NEAR(mrssm.matrixVirt_uubar_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4, divergence, mu), 0.016649076944301452, 0.);
   EXPECT_NEAR(mrssm.matrixVirt_ddbar_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4, divergence, mu), 0.00096447407012237567, 0.);
   EXPECT_NEAR(mrssm.matrixVirt_GG_suLsuLdagger(alphas, S, T, FiniteGs, Dminus4, divergence, mu), 0.0049742005949962483, 0.);
}

}
