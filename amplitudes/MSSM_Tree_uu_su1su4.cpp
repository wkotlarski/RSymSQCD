/* MSSM-tree level */
#include "MSSM_Tree_uu_su1su4.h"

double matrixMSSMTree_uu_suLsuR(double alphaS, double MassSq, double MassGlu, 
                    double S, double T, double U)
{
double msquaredReal = 315.82734083485946*(alphaS*alphaS)*(1/pow(MassGlu*MassGlu - 1.*T,2) + 1/pow(MassGlu*MassGlu - 1.*U,2))*(-1.*pow(MassSq,4) + T*U); /*left-right*/
                    //+ 105.27578027828648*(alphaS*alphaS)*(MassGlu*MassGlu)*S*(3./pow(MassGlu*MassGlu - 1.*T,2) + 3./pow(MassGlu*MassGlu - 1.*U,2) - 2./((MassGlu*MassGlu - 1.*T)*(MassGlu*MassGlu - 1.*U))); /* left-left + right-right or 1/2 of left-left */

return msquaredReal;
}
