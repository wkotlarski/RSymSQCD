#include "MRSSM_Tree_uu_su1su4.h"

double MsquaredMRSSMTree_uu_suLsuR(double alphaS, double MassSq, double MassGlu, 
                    double T, double U, double S)
{
double MsquaredReal = (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + T,2) + (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(-1.*(MassGlu*MassGlu) + U,2);
return MsquaredReal;
}
