# RSymSQCD

A C++ code for calculation of NLO (S)QCD corrections in (not only) the Minimal R-symmetric Supersymmetric Standard Model (**MRSSM**).

This code was developed as part of the ongoing research into an R-symmetric SUSY.
If you use it, please cite:

* P. Diessner, W. Kotlarski, S. Liebschner and D. Stöckinger *Squark production in R-symmetric SUSY with Dirac gluinos: NLO corrections*, J. High Energ. Phys. (2017) 2017: 142 ([inspirehep](https://inspirehep.net/literature/1610032))

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

To compile RSymSQCD you need a C++17 compiler (icpc, g++ >= 7.1 or clang >= 5.0).

Location of dependecies can be passed to cmake cmake via `CMAKE_PREFIX_PATH`
   viariable. It accepts a `;` separated list of directories. Don't forget to
   escape the semicolon in the shell `\;`!

There is one project specific CMake variable: `OPTIMIZE_FOR_NATIVE` which enables (default is `OFF`) CPU specific optimizations (`-march=native` on GCC or Clang and `-xHost` on Intel).
Use it if you plan to run the code on the same (or reasonably simalar) machine as the one you are compiling it on.

## Authors
Wojciech Kotlarski and Sebastian Liebschner


## Licence
This code is distributed under the terms of GPL3 licence.
