![tests](https://github.com/wkotlarski/RSymSQCD/actions/workflows/test.yml/badge.svg)

# RSymSQCD

A C++ code for calculation of NLO (S)QCD corrections in (not only) the Minimal R-symmetric Supersymmetric Standard Model (**MRSSM**).

This code was developed as part of the ongoing research into an R-symmetric SUSY.
If you use it, please cite:

* **P. Diessner, W. Kotlarski, S. Liebschner and D. St&ouml;ckinger** Squark production in R-symmetric SUSY with Dirac gluinos: NLO corrections, *J. High Energ. Phys. (2017) 2017: 142* ([inspirehep](https://inspirehep.net/literature/1610032))

## Getting started

### Prerequisites

* [boost](http://www.boost.org)
* [cmake](https://cmake.org)
* [Cuba](http://www.feynarts.de/cuba)
* [LoopTools](http://www.feynarts.de/looptools)
* [LHAPDF](https://lhapdf.hepforge.org)
* [rk](http://rk.hepforge.org)
* Threading Building Blocks

Some of those dependencies are usually avaliable in the repository of your linux distribution.
On Ubuntu, they can be installed as
```console
sudo apt install -y nlohmann-json3-dev libspdlog-dev libboost-program-options-dev libtbb-dev
```
LoopTools has to be compiled in position independent mode
```console
FFLAGS="-O3 -fPIC" ./configure
```

### Bulding RSymSQCD

To compile RSymSQCD you need a C++17 compiler (icpc, g++ >= 7.1 or clang >= 5.0).

Location of `LoopTools`, `Cuba`, `rk` can be passed to cmake via `LT_PREFIX`, `CUBA_PREFIX`, `RK_PREFIX` variables 

Example:
```console
cmake -Bbuild -DCUBA_PREFIX=YOUR_CUBA_LOCATION/Cuba-4.2.2 -DLT_PREFIX=YOUR_LT_LOCATION/LoopTools-2.16 -DRK_PREFIX=YOUR_LT_LOCATION/rk-1.8
cmake --build build
```

There is one project specific CMake variable: `OPTIMIZE_FOR_NATIVE` which enables (default is `OFF`) CPU specific optimizations (`-march=native` on GCC or Clang and `-xHost` on Intel).
Use it if you plan to run the code on the same (or reasonably similar) machine as the one you are compiling it on.


## Models

* MRSSM
* MSSM
* SM + colour octet scalar

## Authors
Wojciech Kotlarski and Sebastian Liebschner


## Licence
This code is distributed under the terms of GPL3 licence.

## Acknowledgments

The code includes the $Li_2$ function from [polylogarithm](https://github.com/Expander/polylogarithm) package by Alexander Voigt.
