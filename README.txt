1. Licence
   ask Alex?

2. Instalation
   To compile RSymSQCD you need a C++14 compiler. The code was succesufully compiled with intel icpc >= 16 and g++ >= 4.9.
   scons

2.1 Dependencies
      If the LHAPDF script called lhapdf-config is in your PATH the scons script will automatically find LHAPDF.
    Cuba, LoopTools, LHAPDF, rk, boost

2.2 Compiler settings.

   If possible compile with intel icpc compiler.
   If you plan on running the code on the same computer you compile it, add the -xHost (intel) or -march=native (gcc).
