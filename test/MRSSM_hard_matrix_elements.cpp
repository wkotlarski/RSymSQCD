#include "models/MRSSM.hpp"

#include "gtest/gtest.h"
#include <iomanip>

namespace {

TEST(MRSSMHardMatrixElementTest, BMP1) {
   static std::array<std::array<double, 4>, 5> p {{
      {{3042.4548058144451, 0, 0, 3042.4548058144451}},
         {{3042.4548058144451, 0, 0, -3042.4548058144451}},
         {{3000.4336882842676, -2407.3740392707309, -296.96044906456916, -932.18401859709491}},
         {{2979.7366113194435, 2419.4097890284415, 299.06182291926109, 828.15975008632654}},
         {{104.73931202517923, -12.035749757710732, -2.1013738546919383, 104.02426851076837}}
   }};

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

   // EXPECT_NEAR(mrssm.matrixHard_uu_suLsuRg(alphas, p), 0.020376409677647898, 5e-17);
   EXPECT_NEAR(mrssm.matrixHard_gu_suLsuRubar(alphas, p), 1.3897879902109231e-07, 6.3e-21);

   // new
   params.WidthGlu = 1e-1 * params.MassGlu;
   std::cout << std::setprecision(17);
   EXPECT_NEAR(mrssm.matrixHard_gu_suLsuRubar_DS(alphas, p), 1.17995509181456260e-04, 4e-22);
   EXPECT_NEAR(mrssm.matrixHard_gu_suLsuRubar_DS_unsimp(alphas, p), 1.17995274638592898e-04, 4e-22);
   EXPECT_EQ(mrssm.matrixHard_gu_suLsuRubar_DS(alphas, p), mrssm.matrixHard_gu_suLsuRubar_DS_unsimp(alphas, p));
   std::cout << std::scientific << std::setprecision(17) <<mrssm.matrixHard_gu_suLsuRubar_DS(alphas, p) << ' ' << mrssm.matrixHard_gu_suLsuRubar_DS_unsimp(alphas, p) << ' ' << 1-mrssm.matrixHard_gu_suLsuRubar_DS(alphas, p)/mrssm.matrixHard_gu_suLsuRubar_DS_unsimp(alphas, p) << std::endl;
}

}
