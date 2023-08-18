#ifndef SPLITTING_KERNELS_H_
#define SPLITTING_KERNELS_H_

#include <array>

enum class SplittingKernel {Pqq, Pgq, Pgg, Pqg};

std::array<double, 2> get_sp(SplittingKernel sp, double z);

std::array<double, 2> Pqq (double);
std::array<double, 2> Pgq (double);
std::array<double, 2> Pgg (double);
std::array<double, 2> Pqg (double);

#endif
