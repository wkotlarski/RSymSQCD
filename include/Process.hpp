#ifndef PROCESS_H_
#define PROCESS_H_

#include <math.h>
#include <complex>
#include "clooptools.h"
#include <string>
#include <array>
#include <iostream>
#include <boost/property_tree/ptree.hpp>

class Process {
   private:
      // particle masses
      double MasssigmaO, MassphiO, MassGlu, MassTop, MassSq,
         MassSuL, MassSuR, MassSdL, MassSdR, MassSsL, MassSsR,
         MassScL, MassScR, MassSbL, MassSbR, MassStL, MassStR;
       
      // tree-level matrix elements
      double matrixMSSMTree_uu_suLsuR(double, double, double, double);
      double matrixMRSSMTree_uubar_suLsuLdagger(double, double, double, double);
      double matrixMSSMTree_ud_suLsdR(double, double, double, double);
      double matrixMRSSMTree_ud_suLsdR(double, double, double, double);
      double matrixMRSSMTree_GG_suLsuLdagger(double, double, double, double);
      
      // loop-level matrix elements
      double matrixMSSMVirt_uu_suLsuR(double, double, double, double, 
                 double, const double, const double, int);
   public: 
      Process(std::string,  boost::property_tree::ptree); 
      double (Process::* matrixelementTree)(double, double, double, double); // matrix element squared 
      double (Process::* matrixelementVirt)(double, double, double, double, double, double,
                                    double, int); // Re[M^1L M^(B star)]
      double m1, m2;
      double f1,f2; // flavours of initial partons
      double k;     // 1/k = average over initial state colors and helicities
      double h;     // h = sum over initial and final state helicities of fermions (_hel = 0 in FormCalc)
};
#endif /* PROCESS_H_ */
