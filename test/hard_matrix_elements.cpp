#include "Process.hpp"
#include "XSection.hpp"

#include "gtest/gtest.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/program_options.hpp>

namespace {

TEST(HardMatrixElementTest, MRSSM) {
   std::array<std::array<double, 4>, 5> p {{
      {{3042.4548058144451, 0, 0, 3042.4548058144451}},
         {{3042.4548058144451, 0, 0, -3042.4548058144451}},
         {{3000.4336882842676, -2407.3740392707309, -296.96044906456916, -932.18401859709491}},
         {{2979.7366113194435, 2419.4097890284415, 299.06182291926109, 828.15975008632654}},
         {{104.73931202517923, -12.035749757710732, -2.1013738546919383, 104.02426851076837}}
   }};

   boost::property_tree::ptree pt;
   pt.put("collider setup.sqrt_S", 13e+3);
   pt.put("masses.gluino", 1000.);
   pt.put("masses.top", 172.);
   pt.put("masses.pseudoscalar_sgluon", 5000.);
   pt.put("masses.squarks", 1500.);
   pt.put("collider setup.pdf", "MMHT2014nlo68cl");
   pt.put("collider setup.mu_r", 1.5e+3);
   pt.put("collider setup.mu_f", 1.5e+3);
   pt.put("technical parameters.dS", 1.5e+3);
   pt.put("technical parameters.eta_sign", 1.5e+3);
   pt.put("technical parameters.delta", 1.5e+3);
   pt.put("technical parameters.WidthOverMass", 1.5e+3);

   Process process1("MRSSM,uu_suLsuR", pt);
   EXPECT_NEAR((process1.*process1.Process::matrixelementReal_HnonC)(p), 0.020376409677647898, 5e-17);

   Process process2("MRSSM,gu_suLsuR", pt);
   EXPECT_NEAR((process2.*process2.Process::matrixelementReal_HnonC)(p), 1.3897879902109231e-07, 7e-21);
}

}
