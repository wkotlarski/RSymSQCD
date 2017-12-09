#include "SM.h"
#include "constants.hpp"

#include "dilog.hpp"
#include "../../../include/IMatrixElements.h"

#include <complex>
#include "clooptools.h"

double delta = 1e-3;
std::vector<Vec4D<double>> mandelstam_to_p2 (double s, double t) {
   double mass = 5;
   double tp = t - pow(mass,2);
   double p = sqrt(s/4. - pow(mass,2));
   double costh = (s/2. + tp)/(p*sqrt(s));
   std::vector<Vec4D<double>> temp = {
           Vec4D<double> { sqrt(s)/2., 0, 0, sqrt(s)/2.},
           Vec4D<double> { sqrt(s)/2., 0, 0, -sqrt(s)/2.},
           Vec4D<double> { sqrt(s)/2., 0., p*sqrt(1-costh*costh), p*costh},
           Vec4D<double> { sqrt(s)/2., 0., -p*sqrt(1-costh*costh), -p*costh}
   };
   return temp;
}

SM::SM(boost::property_tree::ptree const& ptree) :
      MB_(ptree.get<double>("masses.b")),
      Alfa_(1./ptree.get<double>("couplings.inv_alpha_em")){
}

double SM::BornME(std::vector<Particle> part, EpsOrd ord, double S, double T) const noexcept {
   if (part[0] == Particle::e && part[1] == Particle::ebar && part[2] == Particle::b && part[3] == Particle::bbar) {
      return eebar_bbbar_B(ord, mandelstam_to_p2(S, T));
   }
   else {
      throw std::runtime_error("Requested unknown Born matrix element");
   }
}


double SM::BornCCME(std::vector<Particle> part, int emitter, int spectator, EpsOrd ord, std::vector<Vec4D<double>> const& p) const noexcept {
   if (part[0] == Particle::e && part[1] == Particle::ebar && part[2] == Particle::b && part[3] == Particle::bbar) {
      return eebar_bbbar_CCB(emitter, spectator, ord, 2.*p[0]*p[1], MB_*MB_ - 2.*p[1]*p[3]);
   }

}

double SM::VirtualME(std::vector<Particle> part, EpsOrd ord, double S, double T) const noexcept {
   if (part[0] == Particle::e && part[1] == Particle::ebar && part[2] == Particle::b && part[3] == Particle::bbar) {
      std::vector<Vec4D<double>> const& p = mandelstam_to_p2(S, T);
//      double os1 = eebar_bbbar_V_OS(ord, 2*(p[0]*p[1]), 25-2*(p[0]*p[2]));
      double os2 = eebar_bbbar_V_OS2(ord, mandelstam_to_p2(S, T));
//
      double ms = eebar_bbbar_V_MSbar(ord, mandelstam_to_p2(S, T));
      return ms;
   }
}

double SM::RealME(std::vector<Particle> part, std::vector<Vec4D<double>> const& p) const noexcept {
   if (part[0] == Particle::e && part[1] == Particle::ebar && part[2] == Particle::b && part[3] == Particle::bbar) {
      return eebar_bbbar_R (p);
   }
}

// --------------------------------- e+e- -> bbbar --------------------------------
double SM::eebar_bbbar_B(EpsOrd ord, std::vector<Vec4D<double>> const& p) const noexcept {
   double k12 = p[0]*p[1];
   double k14 = p[0]*p[3];
   double k34 = p[2]*p[3];
   double k44 = p[3]*p[3];
   double k23 = p[1]*p[2];
   double k24 = p[1]*p[3];
   double k33 = p[2]*p[2];
   double k13 = p[0]*p[2];
   const double Alfa2 = pow(Alfa_, 2);
   const double MB2OS = 25;
   switch (ord) {
      case EpsOrd::Eps0:
         return Alfa2*(16*(2*(k14*k23 + k13*k24) + 2*k12*MB2OS)*(pi*pi))/(3.*(k12*k12));
      case EpsOrd::Eps1:
         return Alfa2*(-32*(k34 + MB2OS)*(pi*pi))/(3.*k12);
   }
}

double SM::eebar_bbbar_CCB(int emitter, int spectator, EpsOrd ord, double S, double T) const noexcept {
   const double MB2 = pow(MB_, 2);
   const double U = 2*MB2 - T - S;
   const double Alfa2 = pow(Alfa_, 2);
   switch (ord) {
      case EpsOrd::Eps0:
         return 4*(8*Alfa2*ColorMatrix(emitter, spectator, "[{3,4}]","[{3,4}]")*pow(pi,2)*pow(S,-2)*(2*MB2*(S - T - U) + 2*pow(MB2,2) + pow(T,2) + pow(U,2)))/9.;
      case EpsOrd::Eps1:
         return 4*(-8*Alfa2*ColorMatrix(emitter, spectator, "[{3,4}]","[{3,4}]")*pow(pi,2))/9.;
   }
}

double SM::eebar_bbbar_V_OS(EpsOrd ord, double S, double T) const noexcept {
   using dilogarithm::dilog;
   double alphaS = 0.12;
   double Alfa2 = pow(Alfa_, 2);
   double MB2 = pow(MB_, 2);
   double U = 2 * MB2 - S - T;
   double b = sqrt(1 - 4 * MB2 / S);
   const double muR = 91.188;
   double CF = 4 / 3.;
   // eq. 3.24 of Harris & Owens '02
   const double A1 = -2 * CF * (1 - (1 + b * b) / (2 * b) * log((1 + b) / (1 - b)));
   switch (ord) {
      case EpsOrd::DoublePole:
         return 0.;
      case EpsOrd::SinglePole:
         return eebar_bbbar_B(EpsOrd::Eps0, mandelstam_to_p2(S, T)) * alphaS / (2. * pi) * A1;
      case EpsOrd::Eps0:
         return eebar_bbbar_B(EpsOrd::Eps0, mandelstam_to_p2(S, T)) * alphaS / (2. * pi) * (
               A1 * log(muR*muR / S) +
               CF * (-2 * (1 - (1 + b * b) / (2 * b) * log((1 + b) / (1 - b))) * log(S / MB2)
                     + 3 * b * log((1 + b) / (1 - b)) - 4 + (1 + b * b) / b * (-0.5 * pow(log((1 - b) / (1 + b)), 2)
                                                                               + 2 * log((1 - b) / (1 + b)) *
                                                                                 log(2 * b / (1 + b)) +
                                                                               2 * dilog((1 - b) / (1 + b)) +
                                                                               2 / 3. * pi * pi)
               )
         )
                // eq. 3.25
                + 8 * pi * pi * Alfa2 / S * pow(1 / 3., 2) * alphaS / (2 * pi)
                  *
                  (4 * 3 * CF * (b * b - 1) / b * log((1 - b) / (1 + b)) * (MB2 / S - (T - MB2) * (U - MB2) / (S * S)));
   }
}

/*
 *    There is a numerical difference between
 */
double SM::eebar_bbbar_V_OS2(EpsOrd ord, std::vector<Vec4D<double>> const& p) const noexcept {
   ltini();
   const double muR = 91.188;
   double Alfas = 0.12;
   double Alfa2 = pow(Alfa_, 2);
   double MB2 = pow(MB_, 2);
   double k12 = p[0]*p[1];
   double k14 = p[0]*p[3];
   double k34 = p[2]*p[3];
   double k44 = p[3]*p[3];
   double k23 = p[1]*p[2];
   double k24 = p[1]*p[3];
   double k33 = p[2]*p[2];
   double k13 = p[0]*p[2];
   if(ord == EpsOrd::DoublePole) {
      return 0.;
   }
   else if (ord  == EpsOrd::SinglePole) {
      setuvdiv(0);
      setlambda(-1);
      RealType temp1 = Alfa2 * Alfas *
            (32*pi*(4*(k14*k23 + k13*k24 + k12*MB2)*
                  (Re(B0i(bb0,2*k12,MB2,MB2)) -
                        2*(k34*Re(C0i(cc0,MB2,2*k12,MB2,0,MB2,MB2)) +
                              Re(C0i(cc00,MB2,2*k12,MB2,0,MB2,MB2)) +
                              (k34 - MB2)*(Re(C0i(cc1,MB2,2*k12,MB2,0,MB2,MB2)) +
                                    Re(C0i(cc2,MB2,2*k12,MB2,0,MB2,MB2))))) -
                  8*MB2*(-2*k13*k23 + k12*MB2)*
                        (Re(C0i(cc1,MB2,2*k12,MB2,0,MB2,MB2)) + Re(C0i(cc11,MB2,2*k12,MB2,0,MB2,MB2)) +
                              2*Re(C0i(cc12,MB2,2*k12,MB2,0,MB2,MB2)) +
                              Re(C0i(cc2,MB2,2*k12,MB2,0,MB2,MB2)) + Re(C0i(cc22,MB2,2*k12,MB2,0,MB2,MB2))))
            )/(9.*(k12*k12));

      // add IR part from the Zb to the vertex
      return (double(temp1) - 2 * 2 * eebar_bbbar_B(EpsOrd::Eps0, p) * Alfas / (3 * pi));
   }
   else if (ord == EpsOrd::Eps0) {
      // check UV-finiteness
      setuvdiv(1);
      setdelta(0);

      // finite part
      setmudim(pow(91.188,2));
      setlambda(0);
      //    - triangle
      RealType temp1 = Alfa2 * Alfas *
            (32*pi*(4*(k14*k23 + k13*k24 + k12*MB2)*
                  (Re(B0i(bb0,2*k12,MB2,MB2)) -
                        2*(k34*Re(C0i(cc0,MB2,2*k12,MB2,0,MB2,MB2)) +
                              Re(C0i(cc00,MB2,2*k12,MB2,0,MB2,MB2)) +
                              (k34 - MB2)*(Re(C0i(cc1,MB2,2*k12,MB2,0,MB2,MB2)) +
                                    Re(C0i(cc2,MB2,2*k12,MB2,0,MB2,MB2))))) -
                  8*MB2*(-2*k13*k23 + k12*MB2)*
                        (Re(C0i(cc1,MB2,2*k12,MB2,0,MB2,MB2)) + Re(C0i(cc11,MB2,2*k12,MB2,0,MB2,MB2)) +
                              2*Re(C0i(cc12,MB2,2*k12,MB2,0,MB2,MB2)) +
                              Re(C0i(cc2,MB2,2*k12,MB2,0,MB2,MB2)) + Re(C0i(cc22,MB2,2*k12,MB2,0,MB2,MB2))))
            )/(9.*(k12*k12))
      ;
      //    - CT apart from Dminus4
      RealType ren_const1 =  +2.*eebar_bbbar_B(EpsOrd::Eps0, p) *
            (-2*Alfas*(Re(B0i(bb0,MB2,0,MB2)) + Re(B0i(bb1,MB2,0,MB2)) +
                  2*MB2*(-Re(B0i(dbb0,MB2,0,MB2)) + Re(B0i(dbb1,MB2,0,MB2)))))/(3.*pi)
            ;
      // -----------------------------------------------------

      // eps time poles
      setlambda(-1);

      // temp3 + temp4
      // Dminus4 coeff. times triangle pole
      RealType temp2 = Alfa2 * Alfas * (
            (-128*pi*((k14*k23 + k13*k24 + k12*(k34 + 2*MB2))*Re(B0i(bb0,2*k12,MB2,MB2)) -
                  2*(k12*k34*(k34 + MB2)*Re(C0i(cc0,MB2,2*k12,MB2,0,MB2,MB2)) +
                        (k14*k23 + k13*k24 + k12*k34 + 2*k12*MB2)*
                              Re(C0i(cc00,MB2,2*k12,MB2,0,MB2,MB2)) +
                        k12*(k34 - MB2)*(k34 + MB2)*Re(C0i(cc1,MB2,2*k12,MB2,0,MB2,MB2)) +
                        MB2*(-2*k13*k23 + k12*MB2)*Re(C0i(cc11,MB2,2*k12,MB2,0,MB2,MB2)) -
                        4*k13*k23*MB2*Re(C0i(cc12,MB2,2*k12,MB2,0,MB2,MB2)) +
                        2*k12*(MB2*MB2)*Re(C0i(cc12,MB2,2*k12,MB2,0,MB2,MB2)) +
                        k12*(k34*k34)*Re(C0i(cc2,MB2,2*k12,MB2,0,MB2,MB2)) -
                        k12*(MB2*MB2)*Re(C0i(cc2,MB2,2*k12,MB2,0,MB2,MB2)) -
                        2*k13*k23*MB2*Re(C0i(cc22,MB2,2*k12,MB2,0,MB2,MB2)) +
                        k12*(MB2*MB2)*Re(C0i(cc22,MB2,2*k12,MB2,0,MB2,MB2)))))/(9.*(k12*k12))
      );
      // Dminus4 x the ren. const. pole time 4d born
      RealType ren_const2 =  +2*eebar_bbbar_B(EpsOrd::Eps0, p) *
            (2*Alfas*(Re(B0i(bb0,MB2,0,MB2)) + Re(B0i(bb1,MB2,0,MB2)) +
                  2*MB2*Re(B0i(dbb1,MB2,0,MB2))))/(3.*pi)
      ;
      // Dd born times  ren const pole
      RealType temp5 =  +2*eebar_bbbar_CCB(-1,-1,EpsOrd::Eps1, 2*k12, MB2 - 2*k13) *
            (-2*Alfas*(Re(B0i(bb0,MB2,0,MB2)) + Re(B0i(bb1,MB2,0,MB2)) +
                  2*MB2*(-Re(B0i(dbb0,MB2,0,MB2)) + Re(B0i(dbb1,MB2,0,MB2)))))/(3.*pi)
      ;

      // temp3 + temp4 have very simple analytic form
      assert(
         abs(
       1 - (ren_const1 + ren_const2)/(2*eebar_bbbar_B(EpsOrd::Eps0, p) * (-Alfas/pi) * (4/3. + log(muR*muR/MB2)))
      ) < 1e-15
      );
      return ren_const1 + ren_const2 + temp1 + temp2 + temp5;
   }
}

double SM::eebar_bbbar_V_MSbar(EpsOrd ord, std::vector<Vec4D<double>> const& p) const noexcept {
   ltini();
   const double muR = 91.188;
   double Alfas = 0.12;
   double Alfa2 = pow(Alfa_, 2);
   const double MB2OS {25.};
   const double MBOS {sqrt(MB2OS)};
   double MB2 = pow(MB_-delta, 2);
   double MB = sqrt(MB2);
   double k12 = p[0]*p[1];
   double k14 = p[0]*p[3];
   double k34 = p[2]*p[3];
   double k44 = p[3]*p[3];
   double k23 = p[1]*p[2];
   double k24 = p[1]*p[3];
   double k33 = p[2]*p[2];
   double k13 = p[0]*p[2];
   if(ord == EpsOrd::DoublePole) {
      return 0.;
   }
   else if (ord  == EpsOrd::SinglePole) {
      setuvdiv(0);
      setlambda(-1);
      RealType temp1 = Alfa2 * Alfas *
            (32*pi*(4*(k14*k23 + k13*k24 + k12*MB2)*
                  (Re(B0i(bb0,2*k12,MB2,MB2)) -
                        2*(k34*Re(C0i(cc0,MB2,2*k12,MB2,0,MB2,MB2)) +
                              Re(C0i(cc00,MB2,2*k12,MB2,0,MB2,MB2)) +
                              (k34 - MB2)*(Re(C0i(cc1,MB2,2*k12,MB2,0,MB2,MB2)) +
                                    Re(C0i(cc2,MB2,2*k12,MB2,0,MB2,MB2))))) -
                  8*MB2*(-2*k13*k23 + k12*MB2)*
                        (Re(C0i(cc1,MB2,2*k12,MB2,0,MB2,MB2)) + Re(C0i(cc11,MB2,2*k12,MB2,0,MB2,MB2)) +
                              2*Re(C0i(cc12,MB2,2*k12,MB2,0,MB2,MB2)) +
                              Re(C0i(cc2,MB2,2*k12,MB2,0,MB2,MB2)) + Re(C0i(cc22,MB2,2*k12,MB2,0,MB2,MB2))))
            )/(9.*(k12*k12));

      // add IR part from the Zb to the vertex
      return (double(temp1) - 2 * 2 * eebar_bbbar_B(EpsOrd::Eps0, p) * Alfas / (3 * pi));
   }
   else if (ord == EpsOrd::Eps0) {
      setuvdiv(1);
      // finite part
      setmudim(pow(91.188,2));
      setlambda(0);
      //    - triangle
      RealType triangle_finite = Alfa2 * Alfas *
            (32*pi*(4*MB2*(-4*k13*k23 + 2*k12*MB2OS)*
                  (Re(C0i(cc0,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                        2*Re(C0i(cc1,MB2OS,2*k12,MB2OS,0,MB2,MB2))) -
                  4*MB2OS*(k14*k23 + k13*(-2*k23 + k24) + k12*(-k34 + MB2OS))*
                        (Re(C0i(cc1,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              Re(C0i(cc11,MB2OS,2*k12,MB2OS,0,MB2,MB2))) -
                  4*(k14*k23*k34 + k13*k24*k34 - k13*k23*MB2OS - k14*k24*MB2OS)*
                        (Re(C0i(cc0,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              Re(C0i(cc1,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2))) -
                  2*MB2OS*(k14*(k23 - 2*k24) + k13*k24 + k12*(-k34 + MB2OS))*
                        (Re(C0i(cc0,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              Re(C0i(cc1,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              2*Re(C0i(cc12,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2))) -
                  2*MB2OS*(k14*k23 + k13*(-2*k23 + k24) + k12*(-k34 + MB2OS))*
                        (Re(C0i(cc0,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              Re(C0i(cc1,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              2*Re(C0i(cc12,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2))) +
                  4*MB2*(-4*k14*k24 + 2*k12*MB2OS)*
                        (Re(C0i(cc0,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              2*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2))) -
                  2*(k14*k23 + k13*k24 + k12*MB2OS)*
                        (-2*Re(B0i(bb0,2*k12,MB2,MB2)) +
                              2*(k34 + MB2)*Re(C0i(cc0,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              4*Re(C0i(cc00,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                              2*(k13 + k23 - 2*k34)*Re(C0i(cc1,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              2*(k14 + k24 - 2*MB2OS)*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2))) -
                  4*MB2OS*(k14*(k23 - 2*k24) + k13*k24 + k12*(-k34 + MB2OS))*
                        (Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                              Re(C0i(cc22,MB2OS,2*k12,MB2OS,0,MB2,MB2)))))/(9.*(k12*k12))
      ;

      // eps time poles
      setlambda(-1);
      // Dminus4 coeff. times triangle pole
      RealType triangle_pole_times_epsilon = Alfa2 * Alfas * (
            (-128*pi*((k14*k23 + k13*k24 + k12*k34 + 2*k12*MB2OS)*Re(B0i(bb0,2*k12,MB2,MB2)) -
                  (2*k12*(k34*k34) + 2*k13*k23*MB2 + 2*k14*k24*MB2 -
                        k12*(k13 + k14 + k23 + k24)*MB2 + k14*k23*(MB2 + MB2OS) +
                        k13*k24*(MB2 + MB2OS) + k12*k34*(MB2 + MB2OS))*
                        Re(C0i(cc0,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  2*(k14*k23 + k13*k24 + k12*k34 + 2*k12*MB2OS)*
                        Re(C0i(cc00,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  (-(k13*k13*k24) + k14*k23*(-k23 + k34) +
                        k13*(-(k14*k23) + (k12 + k24)*k34 - k23*(k24 + 4*MB2 - 2*MB2OS)) +
                        k12*((k23 - 3*k34)*k34 + 2*MB2*MB2OS + MB2OS*MB2OS))*
                        Re(C0i(cc1,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  2*k13*k23*MB2OS*Re(C0i(cc11,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  k14*k23*MB2OS*Re(C0i(cc11,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  k13*k24*MB2OS*Re(C0i(cc11,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  k12*k34*MB2OS*Re(C0i(cc11,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  k12*(MB2OS*MB2OS)*Re(C0i(cc11,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  2*k13*k23*MB2OS*Re(C0i(cc12,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  2*k14*k23*MB2OS*Re(C0i(cc12,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  2*k13*k24*MB2OS*Re(C0i(cc12,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  2*k14*k24*MB2OS*Re(C0i(cc12,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  2*k12*k34*MB2OS*Re(C0i(cc12,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  2*k12*(MB2OS*MB2OS)*Re(C0i(cc12,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  k14*k14*k23*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  k13*k14*k24*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  k14*k23*k24*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  k13*(k24*k24)*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  k12*k14*k34*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  k14*k23*k34*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  k12*k24*k34*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  k13*k24*k34*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  k12*(k34*k34)*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  4*k14*k24*MB2*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  2*k14*k23*MB2OS*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  2*k13*k24*MB2OS*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  2*k14*k24*MB2OS*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  2*k12*k34*MB2OS*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  2*k12*MB2*MB2OS*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) +
                  k12*(MB2OS*MB2OS)*Re(C0i(cc2,MB2OS,2*k12,MB2OS,0,MB2,MB2)) -
                  MB2OS*(k14*k23 + k13*k24 - 2*k14*k24 - k12*k34 + k12*MB2OS)*
                        Re(C0i(cc22,MB2OS,2*k12,MB2OS,0,MB2,MB2))))/(9.*(k12*k12))
      );
      RealType temp6 = 2*eebar_bbbar_B(EpsOrd::Eps1, p) *
            (-Alfas/(3.*pi));
      setlambda(0);
      RealType temp3 = -2*eebar_bbbar_B(EpsOrd::Eps0, p) *
            (Alfas*(2*Re(B0i(bb0,MB2OS,0,MB2)) + 2*Re(B0i(bb1,MB2OS,0,MB2)) +
                  2*(2*MB2OS - 4*MB*MBOS)*Re(B0i(dbb0,MB2OS,0,MB2)) +
                  4*MB2OS*Re(B0i(dbb1,MB2OS,0,MB2))))/(3.*pi);
      setlambda(-1);
      RealType temp4 = -2*eebar_bbbar_B(EpsOrd::Eps0, p) *
            (-2*Alfas*(Re(B0i(bb0,MB2OS,0,MB2)) + Re(B0i(bb1,MB2OS,0,MB2)) +
                  2*MB2OS*Re(B0i(dbb0,MB2OS,0,MB2)) - 2*MB*MBOS*Re(B0i(dbb0,MB2OS,0,MB2)) +
                  2*MB2OS*Re(B0i(dbb1,MB2OS,0,MB2))))/(3.*pi);
      RealType temp5 = -2*eebar_bbbar_B(EpsOrd::Eps1, p) *
            (Alfas*(2*Re(B0i(bb0,MB2OS,0,MB2)) + 2*Re(B0i(bb1,MB2OS,0,MB2)) +
                  2*(2*MB2OS - 4*MB*MBOS)*Re(B0i(dbb0,MB2OS,0,MB2)) +
                  4*MB2OS*Re(B0i(dbb1,MB2OS,0,MB2))))/(3.*pi);
//      std::cout << temp5/temp6 << std::endl;
      return triangle_finite + triangle_pole_times_epsilon + temp5 + temp3 + temp4;
   }
}

double SM::eebar_bbbar_R (std::vector<Vec4D<double>> const& p) const noexcept {
   double Alfas = 0.12;
   double Alfa2 = pow(Alfa_, 2);
   double k12 = p[0]*p[1];
   double k35 = p[2]*p[4];
   double k25 = p[1]*p[4];
   double k24 = p[1]*p[3];
   double k15 = p[0]*p[4];
   double k45 = p[3]*p[4];
   double k13 = p[0]*p[2];
   double k23 = p[1]*p[2];
   double k14 = p[0]*p[3];
   double k34 = p[2]*p[3];
   double k44 = p[3]*p[3];
   double k33 = p[2]*p[2];
   double MB2OS {25.};
   double MB2 = pow(5 - delta, 2);
   return
         Alfa2*Alfas*pow(k12,-2)*((1024*((2*(-2*k15*k23 - 2*k15*k25 - 2*k13*(k23 + k25) +
               k12*(k13 + k15 + k23 + k25))*k35 +
               k12*((-2*k35 - MB2)*MB2OS + (MB2 - MB2OS)*MB2OS) +
               MB2*(4*k14*k24 + k12*(-4*k14 - 4*k24 - MB2OS) + 2*pow(k12,2)))*
               pow(2*k35 - MB2 + MB2OS,-2) +
               ((k34 + k35 + MB2OS)*(-2*k14*k45 + k15*MB2OS) +
                     k12*(-2*k45*MB2OS + (MB2 - MB2OS)*MB2OS) -
                     MB2*(-4*k13*k23 + k15*k23 + k13*k25 + 2*k12*MB2OS + 2*pow(k12,2)) +
                     k13*(4*k14*k45 - 2*k34*k45 - 2*k15*MB2OS + k35*MB2OS - k45*MB2OS -
                           2*pow(k45,2)))*pow(2*k45 - MB2 + MB2OS,-2)))/9. -
               (512*((k14*k23 + k12*(k13 + k23) + k13*(-2*k23 + k24))*(MB2 - 2*MB2OS) -
                     4*(-2*k15*k25*(MB2 + MB2OS) + k12*(k34*(MB2 - MB2OS) + MB2*MB2OS)) -
                     MB2*(6*k13*k23 - 5*k14*k23 - 5*k13*k24 - 4*k14*k24 +
                           k12*(3*k13 + 3*k23 + 4*k34 - 8*MB2OS) + 4*pow(k12,2)) +
                     2*(2*k15*k23*k25 + 4*k15*k23*k35 - 2*k15*k23*MB2OS - 4*k15*k25*MB2OS -
                           2*k15*k25*(MB2 + 2*MB2OS) +
                           MB2OS*(-2*k15*k23 + 2*k13*(-2*k23 - k25) +
                                 k12*(k15 + k25 + 2*(k35 + MB2OS))) + 2*k25*pow(k13,2) -
                           2*k23*pow(k15,2) + 2*k15*pow(k23,2) -
                           k13*(2*k15*(k23 - k25) + 2*k23*k25 - 8*k23*k35 - 4*k25*k35 + 2*k25*MB2OS +
                                 2*pow(k25,2)) + k12*
                           (2*k15*k35 + 2*k25*k35 - 2*k13*(k15 - 2*MB2OS) - 2*k23*(k25 - 2*MB2OS) +
                                 3*k15*MB2OS + 3*k25*MB2OS - 4*k35*MB2OS - 4*pow(k13,2) - 4*pow(k23,2) -
                                 4*pow(k35,2) - 2*pow(MB2OS,2))))*pow(2*k35 - MB2 + MB2OS,-1)*
                     pow(2*k45 - MB2 + MB2OS,-1))/9.)*pow(pi,3);
//         -4*Alfa2*Alfas*pow(1/(2.*k12),2)*((1024*pow(1/(k33 + 2*k35 - MB2),2)*
//               (-2*(-2*k15*k23 - 2*k15*k25 - 2*k13*(k23 + k25) + k12*(k13 + k15 + k23 + k25))*
//                     k35 + MB2*(-4*k14*k24 + k12*(-2*k13 - 2*k23 + 2*k33 + 2*k34 + k44) +
//                     2*pow(k12,2)) + k12*pow(MB2,2)))/9. -
//               (512*((-2*(2*k15*k23*k25 - 2*k15*k23*k33 - 4*k15*k25*k33 + 4*k15*k23*k35 +
//                     (6*k14*k23 + 6*k13*k24 + 4*k14*k24 -
//                           k12*(3*k13 + k14 + 3*k23 + k24 - 3*k33 + 6*k34 + k44))*MB2 +
//                     2*k25*pow(k13,2) - 2*k23*pow(k15,2) + 2*k15*pow(k23,2) -
//                     k13*(2*k15*(k23 - k25) + 2*k23*k25 + 2*k25*k33 - 8*k23*k35 -
//                           4*k25*k35 + 2*pow(k25,2)) +
//                     k12*(-2*k13*(k15 - 2*k33) - 2*k23*(k25 - 2*k33) + 3*k15*k33 +
//                           3*k25*k33 + 2*k15*k35 + 2*k25*k35 - 4*k33*k35 - 4*pow(k13,2) -
//                           4*pow(k23,2) - 2*pow(k33,2) - 4*pow(k35,2))) + 4*k12*pow(MB2,2))/
//                     ((k33 + 2*k35 - MB2)*(k44 + 2*k45 - MB2)) +
//                     pow(1/(k44 + 2*k45 - MB2),2)*
//                           (2*((k33 + k34 + k35)*(k15*k44 - 2*k14*k45) +
//                                 k13*(-2*k15*k44 + k35*k44 + 4*k14*k45 - 2*k34*k45 - k44*k45 -
//                                       2*pow(k45,2))) -
//                                 2*(MB2*(-6*k13*k23 - k14*k23 - k13*k24 - k12*(k13 + k23 - 2*k33 + k44) +
//                                       4*pow(k12,2)) + k12*pow(MB2,2)))))/9.)*pow(pi,3);
}

// --------------------------------- uubar -> e+e- --------------------------------
