#include <valarray>
#include <string>

#include <constants.hpp>
#include "ColorFull/Core/Tree_level_gluon_basis.h"

#include "CS_Dipole.hpp"

class CS_FF_Dipole : public CS_Dipole {

public:
    explicit CS_FF_Dipole(
            double (*)(std::vector<Vec4D<double>> const&),
            int, int, int,
            std::vector<Vec4D<double>> const&
    );
    virtual double get_dipoles_value(std::vector<Vec4D<double>> const&) override final;

    // @todo
    double unsu(const std::vector<Vec4D<double>>&);

    double Born(const std::vector<Vec4D<double>>&);
private:
   // to costruct a subtraction counter term we need ME with color infromation
   // and dipole indices i j and spectaror index k
   double (*born_me2_)(std::vector<Vec4D<double>> const&);
   double mi = 0.;
   double mj = 1500.;
   double mk = 1500.;
   double Alfa = 1/137.;
   double Alfa2 = pow(Alfa, 2);
   double Alfas = 0.1184;

   virtual std::vector<Vec4D<double>> mom_shuffle() override final;
   double ColorMatrix(const std::string&, const std::string&);
   double yijk();
   double z(int);
   double Vijk();


    Vec4D<double> q_;
};
