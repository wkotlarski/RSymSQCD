#include "SM.h"
#include "constants.hpp"

#include "dilog.hpp"
#include "../../../include/IMatrixElements.h"

#include <complex>
#include "clooptools.h"

double delta = 0E-4;
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
//      return eebar_bbbar_V_OS(ord, S, T);
      //std::cout << eebar_bbbar_V_MSbar(ord, mandelstam_to_p2(S, T))/eebar_bbbar_V_OS(ord, S, T) << std::endl;
      double os = eebar_bbbar_V_OS(ord, S, T);
      double ms = eebar_bbbar_V_MSbar(ord, mandelstam_to_p2(S, T));
//      std::cout << os << ' ' << os/ms << ' ' << ms - os << std::endl;
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
   const double MB2 = pow(5-delta, 2);
   switch (ord) {
      case EpsOrd::Eps0:
         //return 4. * (8 * Alfa2 * (pi * pi) * (S * S + 2 * (MB2 * MB2 - 2 * MB2 * T + T * (S + T)))) / (3. * (S * S));
         return 42.666666666666664*Alfa2*(k14*k23 + k13*k24 + k12*MB2)*pow(0.5/k12,2.)*pow(pi,2.);
      case EpsOrd::Eps1:
         return Alfa2*(-105.27578027828648*(k34 + MB2))/k12;
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
         return 4*(8*Alfa2*ColorMatrix(emitter, spectator, "[{3,4}]","[{3,4}]")*pow(pi,2))*MB2/(9.*S) * (-1);
   }
}

double SM::eebar_bbbar_V_OS(EpsOrd ord, double S, double T) const noexcept {
   using dilogarithm::dilog;
   double alphaS = 0.12;
   double Alfa2 = pow(Alfa_, 2);
   double MB2 = pow(MB_, 2);
   double U = 2 * MB2 - S - T;
   double b = sqrt(1 - 4 * MB2 / S);
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
               A1 * log(91.188 * 91.188 / S) +
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

double SM::eebar_bbbar_V_MSbar(EpsOrd ord, std::vector<Vec4D<double>> const& p) const noexcept {
   ltini();
   setmudim (91.188*91.188);
   double Alfas = 0.12;
   double Alfa2 = pow(Alfa_, 2);
   double MB2 = pow(5, 2);
   double k12 = p[0]*p[1];
   double k14 = p[0]*p[3];
   double k34 = p[2]*p[3];
   double k44 = p[3]*p[3];
   double k23 = p[1]*p[2];
   double k24 = p[1]*p[3];
   double k33 = p[2]*p[2];
   double k13 = p[0]*p[2];
   if(ord == EpsOrd::DoublePole)
         return 0.;
   if (ord  == EpsOrd::SinglePole) {
      setlambda(-1);
      setuvdiv(0);
      std::complex<double> temp1 = Alfa2 * Alfas * (22.340214425527417 * ((-2. * (-2. * k13 * k23 + k12 * MB2) *
                                                                           (-1. * (k34 + MB2) * A0i(aa0, MB2) +
                                                                            MB2 * ((k12 - 2. * (k34 + MB2)) *
                                                                                   B0i(bb0, 2. * k12, MB2, MB2) +
                                                                                   2. * (k34 + MB2) *
                                                                                   B0i(bb0, MB2, 0., MB2)))) /
                                                                          (-1. * (k34 * k34) + MB2 * MB2) +
                                                                          (k14 * k23 + k13 * k24 + k12 * MB2) *
                                                                          (-3. * B0i(bb0, 2. * k12, MB2, MB2) +
                                                                           4. * B0i(bb0, MB2, 0., MB2) -
                                                                           4. * k34 *
                                                                           C0i(cc0, MB2, 2. * k12, MB2, 0., MB2,
                                                                               MB2)))) / (k12 * k12);
      // add IR part from the Zb to the vertex
      return (temp1.real() - 2 * 2 * eebar_bbbar_B(EpsOrd::Eps0, p) * Alfas / (3 * pi));
   }
   if (ord == EpsOrd::Eps0) {
      // check UV-finiteness
      setuvdiv(0);
      setdelta(0e+6);
      setmudim(91.188*91.188);

      // finite part
      setlambda(0);
      //    - triangle
      std::complex<double> temp1 = Alfa2 * Alfas * (22.340214425527417 * ((-2. * (-2. * k13 * k23 + k12 * MB2) *
                                                                           (-1. * (k34 + MB2) * A0i(aa0, MB2) +
                                                                            MB2 * ((k12 - 2. * (k34 + MB2)) *
                                                                                   B0i(bb0, 2. * k12, MB2, MB2) +
                                                                                   2. * (k34 + MB2) *
                                                                                   B0i(bb0, MB2, 0., MB2)))) /
                                                                          (-1. * (k34 * k34) + MB2 * MB2) +
                                                                          (k14 * k23 + k13 * k24 + k12 * MB2) *
                                                                          (-3. * B0i(bb0, 2. * k12, MB2, MB2) +
                                                                           4. * B0i(bb0, MB2, 0., MB2) -
                                                                           4. * k34 *
                                                                           C0i(cc0, MB2, 2. * k12, MB2, 0., MB2,
                                                                               MB2)))) / (k12 * k12);
      //    - CT apart from Dminus4
      std::complex<double> temp3 =  +2.*eebar_bbbar_B(EpsOrd::Eps0, p) *
            (-0.2122065907891938*Alfas*(Re(B0i(bb0,MB2,0.,MB2)) + Re(B0i(bb1,MB2,0.,MB2)) +
                                                                     2.*MB2*(-1.*Re(B0i(dbb0,MB2,0.,MB2)) + Re(B0i(dbb1,MB2,0.,MB2))))
            );
      // -----------------------------------------------------

      // eps time poles
      setlambda(-1);

      // Dminus4 coeff. times triangle pole
      std::complex<double> temp2 = Alfa2 * Alfas *
            (
                  (-22.340214425527417*((k34 - 1.*MB2)*(2.*k14*k23 + 2.*k13*k24 - 1.*k12*(3.*k34 + MB2))*
                                        B0i(bb0,2.*k12,MB2,MB2) - 4.*
                                                                  (k12*(-1.*(k34*k34) + MB2*MB2)*B0i(bb0,MB2,0.,MB2) +
                                                                   k12*k34*(k34*k34 - 1.*(MB2*MB2))*C0i(cc0,MB2,2.*k12,MB2,0.,MB2,MB2) +
                                                                   MB2*(2.*k13*k23 - 1.*k12*MB2)*
                                                                   (2.*C0i(cc00,MB2,2.*k12,MB2,0.,MB2,MB2) -
                                                                    1.*(k34 - 1.*MB2)*(C0i(cc11,MB2,2.*k12,MB2,0.,MB2,MB2) +
                                                                                       2.*C0i(cc12,MB2,2.*k12,MB2,0.,MB2,MB2) +
                                                                                       C0i(cc22,MB2,2.*k12,MB2,0.,MB2,MB2))))))/(k12*k12*(k34 - 1.*MB2))
            );
      // Dminus4 x the ren. const. pole time 4d born
      std::complex<double> temp4 =  +2*eebar_bbbar_B(EpsOrd::Eps0, p) *
            (0.2122065907891938*Alfas*(Re(B0i(bb0,MB2,0.,MB2)) + Re(B0i(bb1,MB2,0.,MB2)) +
                                       2.*MB2*Re(B0i(dbb1,MB2,0.,MB2)))
            );
      // Dd born times  ren const pole
      std::complex<double> temp5 =  +2*eebar_bbbar_B(EpsOrd::Eps1, p) *
            (-0.2122065907891938*Alfas*(Re(B0i(bb0,MB2,0.,MB2)) + Re(B0i(bb1,MB2,0.,MB2)) +
                                        2.*MB2*(-1.*Re(B0i(dbb0,MB2,0.,MB2)) + Re(B0i(dbb1,MB2,0.,MB2))))
            );
//      std::cout << temp1.real() << ' ' << temp2.real() << ' ' << temp3.real() << ' ' << temp4.real() <<
//          ' ' << temp5.real() << " sum: " << (temp1+temp2+temp3+temp4+temp5).real() << std::endl;
      return (1.*(temp1+temp3) + 1.*(temp2+temp4+temp5)).real();
   }
//   std::cout << "A " << temp1.real()  << ' ' <<  3*eebar_bbbar_B(EpsOrd::Eps0, p) * Alfas/(2*pi) << std::endl;
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
   double MB2 = pow(5 - delta, 2);
   return -4*Alfa2*Alfas*pow(1/(2.*k12),2)*((1024*pow(1/(k33 + 2*k35 - MB2),2)*
                                             (-2*(-2*k15*k23 - 2*k15*k25 - 2*k13*(k23 + k25) + k12*(k13 + k15 + k23 + k25))*
                                              k35 + MB2*(-4*k14*k24 + k12*(-2*k13 - 2*k23 + 2*k33 + 2*k34 + k44) +
                                                         2*pow(k12,2)) + k12*pow(MB2,2)))/9. -
                                            (512*((-2*(2*k15*k23*k25 - 2*k15*k23*k33 - 4*k15*k25*k33 + 4*k15*k23*k35 +
                                                       (6*k14*k23 + 6*k13*k24 + 4*k14*k24 -
                                                        k12*(3*k13 + k14 + 3*k23 + k24 - 3*k33 + 6*k34 + k44))*MB2 +
                                                       2*k25*pow(k13,2) - 2*k23*pow(k15,2) + 2*k15*pow(k23,2) -
                                                       k13*(2*k15*(k23 - k25) + 2*k23*k25 + 2*k25*k33 - 8*k23*k35 -
                                                            4*k25*k35 + 2*pow(k25,2)) +
                                                       k12*(-2*k13*(k15 - 2*k33) - 2*k23*(k25 - 2*k33) + 3*k15*k33 +
                                                            3*k25*k33 + 2*k15*k35 + 2*k25*k35 - 4*k33*k35 - 4*pow(k13,2) -
                                                            4*pow(k23,2) - 2*pow(k33,2) - 4*pow(k35,2))) + 4*k12*pow(MB2,2))/
                                                  ((k33 + 2*k35 - MB2)*(k44 + 2*k45 - MB2)) +
                                                  pow(1/(k44 + 2*k45 - MB2),2)*
                                                  (2*((k33 + k34 + k35)*(k15*k44 - 2*k14*k45) +
                                                      k13*(-2*k15*k44 + k35*k44 + 4*k14*k45 - 2*k34*k45 - k44*k45 -
                                                           2*pow(k45,2))) -
                                                   2*(MB2*(-6*k13*k23 - k14*k23 - k13*k24 - k12*(k13 + k23 - 2*k33 + k44) +
                                                           4*pow(k12,2)) + k12*pow(MB2,2)))))/9.)*pow(pi,3);
}

// --------------------------------- uubar -> e+e- --------------------------------
