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

### Installation
2. Instalation
   To compile RSymSQCD you need a C++14 compiler. The code was succesufully compiled with intel icpc >= 16 and g++ >= 4.9.

2.1 Dependencies
   Location of dependecies can be passed to cmake cmake via `CMAKE_PREFIX_PATH`
   viariable. It accepts a `;` separated list of directories. Don't forget to
   escape the semicolon in the shell `\;`!

2.2 Compiler settings.

   If you plan on running the code on the same computer you compile it, add the -xHost (intel) or -march=native (gcc).

## Authors
Wojciech Kotlarski and Sebastian Liebschner


## Licence
This code is distributed under the terms of GPL3 licence.
