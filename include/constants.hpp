#ifndef CONSTANTS_HPP_
#define CONSTANTS_HPP_

#include <boost/math/constants/constants.hpp>

// mathematical constants

// http://www.boost.org/doc/libs/1_60_0/libs/math/doc/html/math_toolkit/constants.html
using boost::math::double_constants::pi;
using boost::math::double_constants::pi_sqr;
using boost::math::double_constants::pi_cubed;
using boost::math::double_constants::two_pi;
using boost::math::double_constants::euler;
using boost::math::double_constants::e;

// physics

constexpr double to_fb {3.893793656e+11};
// SU(3) group factors
constexpr double CF {4/3.};
constexpr double CA {3.};

#endif /* CONSTANTS_HPP_ */
