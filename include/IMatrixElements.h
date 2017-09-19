#ifndef RSYMSQCD_IMATRIXELEMENTS_H
#define RSYMSQCD_IMATRIXELEMENTS_H

#include <vector>
#include <boost/property_tree/ptree.hpp>
#include "Vec4D.hpp"

#include "LHAPDF/LHAPDF.h"
#include "ColorFull/Core/Tree_level_gluon_basis.h"

/**
 *    @param
 */
enum class EpsOrd {
   DoublePole,
   SinglePole,
   Eps0,
   Eps1,
   Eps2
};

enum class Particle {
   u, ubar, d, dbar, s, sbar, c, cbar, b, bbar,
   e, ebar,
   suL, suLdagger,
   suR, suRdagger
};

class IMatrixElements {

public:
   const LHAPDF::PDF* pdf = LHAPDF::mkPDF( "MMHT2014nlo68cl", 0);
   double mu_r, mu_f;
   virtual double BornME(std::vector<Particle>, double, double) const noexcept = 0;
   virtual double BornCCME(std::vector<Particle>, int, int, EpsOrd, std::vector<Vec4D<double>> const&) const noexcept = 0;
   virtual double VirtualME(std::vector<Particle>, EpsOrd, double, double) const noexcept = 0;
   virtual double RealME(std::vector<Particle>, std::vector<Vec4D<double>> const&) const noexcept = 0;

/**
 * 
 * @return
 */
   static IMatrixElements* create_process(boost::property_tree::ptree const&);

protected:
   /**
    *
    * @param emitter - emitting particle
    * @param spectator
    * @param col_str1
    * @param col_str2
    * @return
    */
   double ColorMatrix(
         int emitter,  int spectator,
         std::string const& col_str1, std::string const& col_str2
   ) const noexcept {

      ColorFull::Col_amp Ca1(col_str1);
      ColorFull::Col_amp Ca2(col_str2);

      ColorFull::Col_functions Col_fun;

      // FormCalc numbers particles from 1, hence i_+1
      if(emitter < 0 && spectator < 0) {
         return Col_fun.double_num(
              Col_fun.scalar_product(Ca1, Ca2)
         );
      } else {
      return Col_fun.double_num(
             Col_fun.scalar_product(
                      // emit gluon from emitter and attach it to the spectator
                      Col_fun.emit_gluon(Ca1, emitter+1, 77), Col_fun.emit_gluon(Ca2, spectator+1, 77)
              )
      );
      }
   }
};
#endif //RSYMSQCD_IMATRIXELEMENTS_H
