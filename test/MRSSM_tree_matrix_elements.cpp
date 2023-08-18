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

   const double Dminus4Coeff = 0;
   EXPECT_DOUBLE_EQ(mrssm.matrixTree_uu_suLsuR(alphas, S, T, Dminus4Coeff),          1.1187184131205632);
   EXPECT_DOUBLE_EQ(mrssm.matrixTree_uubar_suLsuLdagger(alphas, S, T, Dminus4Coeff), 0.41567383692610566);
   EXPECT_DOUBLE_EQ(mrssm.matrixTree_ddbar_suLsuLdagger(alphas, S, T, Dminus4Coeff), 0.15281365356525714);
   EXPECT_DOUBLE_EQ(mrssm.matrixTree_gg_suLsuLdagger(alphas, S, T, Dminus4Coeff),    0.11114746957753696);
}
}
