#include "gtest/gtest.h"

#include "models/MRSSM.hpp"

TEST(MRSSMTreeXSectionsTest, BMP1) {
   MRSSMParameters params;
   params.MassTop = 172.;
   params.MassGlu = 1000.;
   params.MasssigmaO = 5000.;
   params.MassSq = 1500.;
   params.eta_sign = -1;
   params.delta = 0.;
   params.WidthGlu = -1.;

   MRSSM mrssm {params};
   static constexpr double alphas {0.1};
   static constexpr double s {16e+6};

   EXPECT_NEAR(mrssm.sigmaTree_uu_suLsuR(alphas, s), 3.9989030508184933e-10, 0.);
   EXPECT_NEAR(mrssm.sigmaTree_uubar_suLsuLdagger(alphas, s), 2.8107819190543564e-10, 0.);
   EXPECT_NEAR(mrssm.sigmaTree_ddbar_suLsuLdagger(alphas, s), 4.2088476688702667e-11, 0.);
   EXPECT_NEAR(mrssm.sigmaTree_gg_suLsuLdagger(alphas, s), 6.8074076794552077e-11, 0.);
}
