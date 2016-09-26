#include "MRSSM_Tree_ud_su1sd4.h"

double matrixMRSSMTree_ud_suLsdR(double alphaS, double MassSq, double MassGlu, 
                     double T, double U, double S)
{
double msquaredTree = (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2);
return msquaredTree;
}
