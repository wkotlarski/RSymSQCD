#include <valarray>
#include <string>

#include <constants.hpp>
#include "ColorFull/Core/Tree_level_gluon_basis.h"

class CSDipole {
private:
   // to costruct a subtraction counter term we need ME with color infromation
   // and dipole indices i j and spectaror index k
   double (*born_me2_)(std::valarray<std::valarray<double>> const&);
   std::valarray<std::valarray<double>> p_;
   std::valarray<std::valarray<double>> ptilde_;
   std::valarray<double> q_;
   int i_ = -1;
   int j_ = -1;
   int k_ = -1;
   double mi = 0.;
   double mj = 0.;
   double mk = 0.;
   double Alfa = 1/137.;
   double Alfa2 = pow(Alfa, 2);
   double Alfas = 0.1184;

   double x(int);
   std::valarray<std::valarray<double>> ff_mom_shuffle();
   double Mat(const std::string&, const std::string&);
   double dot(const std::valarray<double>&, const std::valarray<double>&);
   double yijk();
   double z(int);
   double Vijk();

public:
   double unsu(const std::valarray<std::valarray<double>>&);
   double Born(const std::valarray<std::valarray<double>>&);
   CSDipole(double (*)(std::valarray<std::valarray<double>> const&), int, int, int, std::valarray<std::valarray<double>> const&);
   double eval_dipole(const std::valarray<std::valarray<double>>&);

   static double lambda(double, double, double);
   double beta(std::valarray<double>, std::valarray<double> );
};
