# RSymSQCD

A C++ code for calculation of NLO (S)QCD corrections in (not only) the Minimal R-symmetric Supersymmetric Standard Model (**MRSSM**).

## Getting started

### Prerequisites

* [cmake](https://cmake.org)
* [boost](http://www.boost.org)
* [eigen](https://eigen.tuxfamily.org)
* [Cuba](http://www.feynarts.de/cuba)
* [LoopTools](http://www.feynarts.de/looptools)
* [rk](http://rk.hepforge.org)
* [LHAPDF](https://lhapdf.hepforge.org)

### Bulding RSymSQCD

1. Instalation
   To compile RSymSQCD you need a C++17 compiler. The code was succesufully compiled with intel icpc >= 16 and g++ >= 4.9.

2. Dependencies
   Location of dependecies can be passed to cmake cmake via `CMAKE_PREFIX_PATH`
   viariable. It accepts a `;` separated list of directories. Don't forget to
   escape the semicolon in the shell `\;`!

There is one project specific CMake variable: `OPTIMIZE_FOR_NATIVE` which enables CPU specific optimizations (`-march=native` on GCC or Clang and `-xHost` on Intel)

## Authors
Wojciech Kotlarski and Sebastian Liebschner


## Licence
This code is distributed under the terms of GPL3 licence.
