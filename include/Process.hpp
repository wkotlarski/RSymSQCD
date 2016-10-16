#ifndef PROCESS_H_
#define PROCESS_H_

#include <cmath>
#include <complex>
#include "clooptools.h"
#include <string>
#include <array>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <mathematica_wrapper.hpp>
#include <constants.hpp>
#include "LHAPDF/LHAPDF.h"

class Process {
   private:
      // particle masses
      double MasssigmaO, MassphiO, MassGlu, MassTop, MassSq,
         MassSuL, MassSuR, MassSdL, MassSdR, MassSsL, MassSsR,
         MassScL, MassScR, MassSbL, MassSbR, MassStL, MassStR, 
         mu_r, mu_f, dS;
      const LHAPDF::PDF* pdf;
       
      // tree-level matrix elements
      double matrixMSSMTree_uu_suLsuR( double );
      double matrixMRSSMTree_uubar_suLsuLdagger(double, double, double, double);
      double matrixMSSMTree_ud_suLsdR(double, double, double, double);
      double matrixMRSSMTree_ud_suLsdR(double, double, double, double);
      double matrixMRSSMTree_GG_suLsuLdagger(double, double, double, double);
      inline double matrixSgluonTree_qqbar_OO(double);
      inline double matrixSgluonTree_gg_OO(double);
      
      // loop-level matrix elements
      double matrixMSSMVirt_uu_suLsuR(double, double, double, double, int);
      
      // soft matrix elements
      double matrixMRSSMSoft_uu_suLsuRg(double, double);
      double f( double, double, double, double, int );
      double g( double, double );
      
   public: 
      Process(std::string,  boost::property_tree::ptree); 
      double (Process::* matrixelementTree)(double); // matrix element squared 
      double (Process::* matrixelementVirt)(double, double, double, //double, double, double,
                                    double, int); // Re[M^1L M^(B star)]
      double (Process::* matrixelementReal_SC)(double, double);
      double (Process::* matrixelementReal_HC1)(double, double, double);
      double (Process::* matrixelementReal_HC2)(double, double, double);
      
      double m1, m2;
      double f1,f2; // flavours of initial partons
      double k;     // 1/k = average over initial state colors and helicities
      double h;     // h = sum over initial and final state helicities of fermions (_hel = 0 in FormCalc)
};
#endif /* PROCESS_H_ */
