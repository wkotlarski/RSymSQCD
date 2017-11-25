# RSymSQCD

A C++ code for calculation of NLO (S)QCD corrections in (not only) the Minimal R-symmetric Supersymmetric Standard Model (**MRSSM**).

## Getting started

### Prerequisites

* [cmake](https://cmake.org) - 
* [boost](http://www.boost.org) -
* [Cuba](http://www.feynarts.de/cuba) -
* [LoopTools](http://www.feynarts.de/looptools) - 
* [rk](http://rk.hepforge.org) - 
* [LHAPDF](https://lhapdf.hepforge.org) -

### Installation
2. Instalation
   To compile RSymSQCD you need a C++14 compiler. The code was succesufully compiled with intel icpc >= 16 and g++ >= 4.9.
   scons

2.1 Dependencies
      If the LHAPDF script called lhapdf-config is in your PATH the scons script will automatically find LHAPDF.

2.2 Compiler settings.

   If possible compile with intel icpc compiler.
   If you plan on running the code on the same computer you compile it, add the -xHost (intel) or -march=native (gcc).

## Authors
Wojciech Kotlarski and Sebastian Liebschner


## Licence
This code is distributed under the terms of GPL3 licence.