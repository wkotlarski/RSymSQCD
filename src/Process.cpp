#include "Process.h"

Process::Process(std::string processID) {	        
    /* squark production, MRSSM */
	if(processID == "MRSSM,uu_suLsuR") {
      // tree-level is the same as in the MSSM
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      matrixelementVirt = &matrixMRSSMVirt_uu_suLsuR;
      f1 = 2.;
      f2 = 2.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,ud_suLsdR") {
      matrixelementTree = &Process::matrixMRSSMTree_ud_suLsdR;
      matrixelementVirt = &matrixMRSSMVirt_ud_suLsdR;
      f1 = 2.;
      f2 = 1.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
    /* squark production, MSSM */
   else if(processID == "MSSM,uu_suLsuR") {
      matrixelementTree = &Process::matrixMSSMTree_uu_suLsuR;
      matrixelementVirt = &matrixMSSMVirt_uu_suLsuR;
      f1 = 2.;
      f2 = 2.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MSSM,ud_suLsdR") {
      matrixelementTree = &Process::matrixMSSMTree_ud_suLsdR;
      matrixelementVirt = &matrixMSSMVirt_ud_suLsdR;
      f1 = 2.;
      f2 = 1.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
//    /* squark anti-squark production, MRSSM */
   else if(processID == "MRSSM,uubar_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_uubar_suLsuLdagger;
      matrixelementVirt = &matrixMRSSMVirt_uubar_suLsuLdagger;
      f1 = 2.;
      f2 = -2.;
      k = 2.*2*3*3;
      h = 2.*2;
   }
   else if(processID == "MRSSM,GG_suLsuLdagger") {
      matrixelementTree = &Process::matrixMRSSMTree_GG_suLsuLdagger;
      matrixelementVirt = &matrixMRSSMVirt_GG_suLsuLdagger;
      f1 = 0.;
      f2 = 0.;
      k = 2.*2*8*8;
      h = 1.;
   }
}

double Process::matrixMSSMTree_uu_suLsuR(double alphaS, double T, double U, double S) {
   double msquaredReal = 315.82734083485946*(alphaS*alphaS)*(1/pow(MassGlu*MassGlu - 1.*T,2) + 
      1/pow(MassGlu*MassGlu - 1.*U,2))*(-1.*pow(MassSq,4) + T*U); /*left-right*/
                    //+ 105.27578027828648*(alphaS*alphaS)*(MassGlu*MassGlu)*S*(3./pow(MassGlu*MassGlu - 1.*T,2) + 3./pow(MassGlu*MassGlu - 1.*U,2) - 2./((MassGlu*MassGlu - 1.*T)*(MassGlu*MassGlu - 1.*U))); /* left-left + right-right or 1/2 of left-left */
   return msquaredReal;
}

double Process::matrixMRSSMTree_uubar_suLsuLdagger(double alphaS, double T, double U, double S) {

   double msquaredtree =  (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/(S*S) + 355.3057584392169*(alphaS*alphaS)*pow(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U) - 236.8705056261446*(alphaS*alphaS)*(0.3333333333333333/S - 1./(-1.*(MassGlu*MassGlu) + T))*(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T))*(-1.*pow(MassSq,4) + T*U) + 355.3057584392169*(alphaS*alphaS)*pow(1/S - 0.3333333333333333/(-1.*(MassGlu*MassGlu) + T),2)*(-1.*pow(MassSq,4) + T*U);

   return msquaredtree;
}

double Process::matrixMSSMTree_ud_suLsdR(double alphaS, double T, double U, double S) {
double msquaredReal = (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(MassGlu*MassGlu - 1.*T,2); /*  ONLY HAVE OF THE VALUE WHICH IS COMMENTED: left-right + right-left or 2 of left-right*/
                   // + 2.*(315.82734083485946*(alphaS*alphaS)*(MassGlu*MassGlu)*S)/pow(MassGlu*MassGlu - 1.*T,2); /* left-left + right-right or 2 of left-left */

return msquaredReal;
}

double Process::matrixMRSSMTree_ud_suLsdR(double alphaS, double T, double U, double S) {
double msquaredTree = (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2);
return msquaredTree;
}

double Process::matrixMRSSMTree_GG_suLsuLdagger(double alphaS, double T, double U, double S) {
    double msquaredtree =  -157.91367041742973*(alphaS*alphaS)*((96.*(pow(MassSq,4) + T*U - 1.*(MassSq*MassSq)*(T + U)))/(S*S) + MassSq*MassSq*(37.333333333333336/(-1.*(MassSq*MassSq) + U) - 85.33333333333333*(T/pow(-1.*(MassSq*MassSq) + T,2) + U/pow(-1.*(MassSq*MassSq) + U,2))) + (37.333333333333336*(MassSq*MassSq) + (5.333333333333333*(9.*pow(MassSq,4) - 3.*(MassSq*MassSq)*(2.*(MassSq*MassSq) + S) + (S + T)*(S + U)))/(-1.*(MassSq*MassSq) + U))/(-1.*(MassSq*MassSq) + T) - (48.*((5.*pow(MassSq,4) + MassSq*MassSq*U - 1.*T*(-1.*(MassSq*MassSq) + U) - 1.*(MassSq*MassSq)*(3.*S + 4.*T + 2.*U))/(-1.*(MassSq*MassSq) + U) + (5.*pow(MassSq,4) + MassSq*MassSq*U - 1.*T*(-1.*(MassSq*MassSq) + U) - 1.*(MassSq*MassSq)*(3.*S + 2.*T + 4.*U))/(-1.*(MassSq*MassSq) + T)))/S);

    return msquaredtree;
}