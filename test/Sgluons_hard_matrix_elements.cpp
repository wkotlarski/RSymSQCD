#include "models/Sgluons.hpp"

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

   SgluonParameters params;
   params.mO = 1000.;

   Sgluons sgluons {params};
   static constexpr double alphas {8.41234775121963707e-02};

}

}
