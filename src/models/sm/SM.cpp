#include "SM.h"
#include "constants.hpp"

#include "dilog.hpp"
#include "../../../include/IModel.h"

SM::SM(boost::property_tree::ptree const& ptree) :
      MB_(ptree.get<double>("masses.b")),
      Alfa_(1./ptree.get<double>("couplings.inv_alpha_em")){
}

double SM::BornME(std::vector<Particle> part, double S, double T) const noexcept {
   if (part[0] == Particle::e && part[1] == Particle::ebar && part[2] == Particle::b && part[3] == Particle::bbar) {
      return eebar_bbbar_B(S, T);
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
      return eebar_bbbar_V(ord, S, T);
   }
}

double SM::RealME(std::vector<Particle> part, std::vector<Vec4D<double>> const& p) const noexcept {
   if (part[0] == Particle::e && part[1] == Particle::ebar && part[2] == Particle::b && part[3] == Particle::bbar) {
      return eebar_bbbar_R (p);
   }
}

// --------------------------------- e+e- -> bbbar --------------------------------
double SM::eebar_bbbar_B(double S, double T) const noexcept {
   const double Alfa2 = pow(Alfa_, 2);
   const double MB2 = pow(MB_, 2);
   return 4.*(8*Alfa2*(pi*pi)*(S*S + 2*(MB2*MB2 - 2*MB2*T + T*(S + T))))/(3.*(S*S));
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

double SM::eebar_bbbar_V(EpsOrd ord, double S, double T) const noexcept {
   using dilogarithm::dilog;
   double alphaS = 0.12;
   double Alfa2 = pow(Alfa_, 2);
   double MB2 = pow(MB_, 2);
   double U = 2 * MB2 - S - T;
   double b = sqrt(1 - 4 * MB2 / S);
   double CF = 4 / 3.;
   const double born =
         4. * (8 * Alfa2 * (pi * pi) * (S * S + 2 * (MB2 * MB2 - 2 * MB2 * T + T * (S + T)))) / (3. * (S * S));
// eq. 3.24 of Harris & Owens '02
   const double A1 = -2 * CF * (1 - (1 + b * b) / (2 * b) * log((1 + b) / (1 - b)));
   switch (ord) {
      case EpsOrd::DoublePole:
         return 0.;
      case EpsOrd::SinglePole:
         return born * alphaS / (2. * pi) * A1;
      case EpsOrd::Eps0:
         return born * alphaS / (2. * pi) * (
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

double SM::eebar_bbbar_R (std::vector<Vec4D<double>> const& p) const noexcept {
   double Alfas = 0.12;
   double Alfa2 = pow(137., -2);
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
   double S12 = 2. * k12;
   double m1 = 5;
   double m2 = 5;
   double S35 = m1 * m1 + 2 * k35;
   double S45 = m2 * m2 + 2 * k45;
   double T = m1 * m1 - 2. * k13;
   double T14 = m2 * m2 - 2. * k14;
   double MB2 = 5*5;
   return 4*(-512*Alfa2*Alfas*((-2*k14*k24*MB2 - k12*k35*MB2 + k12*MB2*(k12 + MB2) +
                                (-2*k14*k24 + k12*(k14 + k24))*(-k13 - k23 + k34 + MB2))/(4.*(k35*k35)) -
                               (k15*k15*k23 - k15*(k23*k23) - 2*(k12*k12)*(k13 + k23) - k13*k13*k25 - k15*k23*k25 - 2*k15*k23*k35 - 2*k14*k23*MB2 +
                                k13*(k15*(k23 - k25) + k25*(k23 + k25) + 2*(-2*k23 - k25)*k35 - 2*k24*MB2) +
                                k12*(4*(k13*k13) + 4*(k23*k23) + k23*(2*k15 + 3*k25 - 2*k35) - (k15 + k25 - 2*k35)*k35 +
                                     k13*(3*k15 - 2*(-2*k23 - k25 + k35)) + (-k15 - k25 + 2*(k34 + k35))*MB2))/(4.*k35*k45) +
                               (k15*(-(k13*k23) + k23*k23) + k13*k13*k25 - k13*k23*k25 + 2*k13*k23*k35 + k12*k12*MB2 - 2*k13*k23*MB2 +
                                k12*(MB2*MB2 + k15*(k13 + MB2) + k25*(k23 + MB2) - k35*(k13 + k23 + MB2)))/(4.*(k45*k45)))*pow(pi,3))/
          (9.*(k12*k12));
}
