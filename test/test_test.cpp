#include "gtest/gtest.h"
#include <cmath>

double square_root(double x) {return sqrt(x);}

TEST(SquareRootTest, PositiveNos) {
    EXPECT_EQ (18.0, square_root (324.0));
    EXPECT_EQ (25.4, square_root (645.16));
    EXPECT_EQ (50.3321, square_root (2533.310224));
}