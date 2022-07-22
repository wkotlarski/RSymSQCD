#include "gtest/gtest.h"

#include "Process.hpp"

namespace {

TEST(TreeMatrixElementTest, MRSSM) {
   boost::property_tree::ptree pt;
   pt.put("masses.gluino", 2000.);
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

   Process p1 = Process("MRSSM,uu_suLsuR", pt);

   const double ME2_1 = (p1.*p1.Process::matrixelementTree)(76858666.754357532, -35567448.864206761);
   EXPECT_NEAR(ME2_1, 3.6118052738604587, 3e-15);
   const double ME2_2 = (p1.*p1.Process::matrixelementTree)(12231711.723910889, -775506.38876226079);
   EXPECT_NEAR(ME2_2, 0.038728580246769444, 2e-17);

   Process p2 = Process("MRSSM,uubar_suLsuLdagger", pt);
   const double ME2_3 = (p2.*p2.Process::matrixelementTree)(45819697.193701074, -3009110.2033250779);
   EXPECT_NEAR(ME2_3,  5.7603581273692486, 4e-15);
   const double ME2_4 = (p2.*p2.Process::matrixelementTree)(17445312.5, -5714031.0653457958);
   EXPECT_NEAR(ME2_4,  1.7101097316789309, 2e-15);

   Process p3 = Process("MRSSM,GG_suLsuLdagger", pt);
   const double ME2_5 = (p3.*p3.Process::matrixelementTree)(1.14422471782008708e+07, -6.02139997033274081e+06);
   EXPECT_NEAR(ME2_5, 2.52442229108592073e+01, 2e-14);
   const double ME2_6 = (p3.*p3.Process::matrixelementTree)(9.31249764834811352e+06, -1.55330597453783639e+06);
   EXPECT_NEAR(ME2_6, 2.17602531175420566e+01, 2e-14);

   Process p4 = Process("MRSSM,ddbar_suLsuLdagger", pt);
   const double ME2_7 = (p4.*p4.Process::matrixelementTree)(4.12440390444746092e+07, -2.64010504270515181e+07);
   EXPECT_NEAR(ME2_7,  7.04258113126766672e-01, 5e-16);
   const double ME2_8 = (p4.*p4.Process::matrixelementTree)(9.24701406368827634e+06, -1.64344817782115540e+06);
   EXPECT_NEAR(ME2_8, 1.98912627103063578e-03, 9e-19);
}

TEST(TreeMatrixElementTest, MSSM) {

}

}
