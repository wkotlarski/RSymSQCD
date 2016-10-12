#ifndef PROCESS_H_
#define PROCESS_H_

#include "MRSSM_1L_uu_su1su4.h"
#include "MSSM_1L_uu_su1su4.h"
#include "MRSSM_1L_ud_su1sd4.h"
#include "MSSM_1L_ud_su1sd4.h"
#include "MRSSM_1L_uubar_su1su1dagger.h"
#include "MRSSM_1L_GG_su1su1dagger.h"

#include <string>
#include <array>

class Process {
   private:
      // particle masses
      double MassGlu = 1.0e+3;
      double MassSq = 1.5e+3;
      double mass_OA  = 1e+3;
      double mass_OS  = sqrt( pow(mass_OA,2) + 4.0 * pow(MassGlu, 2) );
      std::array< std::array<double, 2>, 6 > squark_mass;
      double mass_top;
       
      // tree-level matrix elements
      double matrixMSSMTree_uu_suLsuR(double, double, double, double);
      double matrixMRSSMTree_uubar_suLsuLdagger(double, double, double, double);
      double matrixMSSMTree_ud_suLsdR(double, double, double, double);
      double matrixMRSSMTree_ud_suLsdR(double, double, double, double);
      double matrixMRSSMTree_GG_suLsuLdagger(double, double, double, double);
      
      // loop-level matrix elements
      
   public: 
      Process(std::string); 
      double (Process::* matrixelementTree)(double, double, double, double); // matrix element squared 
      double (*matrixelementVirt)(double, double, double, double, double, double,
                                    double, double, double, double, double, double); // Re[M^1L M^(B star)]
      double m1, m2;
      double f1,f2; // flavours of initial partons
      double k;     // 1/k = average over initial state colors and helicities
      double h;     // h = sum over initial and final state helicities of fermions (_hel = 0 in FormCalc)
};
#endif /* PROCESS_H_ */
