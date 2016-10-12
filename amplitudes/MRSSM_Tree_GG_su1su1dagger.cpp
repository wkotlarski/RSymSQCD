#include "MRSSM_Tree_GG_su1su1dagger.h"

double matrixMRSSMTree_GG_suLsuLdagger(double alphaS, double MassSq, double MassGlu, 
                double T, double U, double S)
{

    double msquaredtree =  -157.91367041742973*(alphaS*alphaS)*((96.*(pow(MassSq,4) + T*U - 1.*(MassSq*MassSq)*(T + U)))/(S*S) + MassSq*MassSq*(37.333333333333336/(-1.*(MassSq*MassSq) + U) - 85.33333333333333*(T/pow(-1.*(MassSq*MassSq) + T,2) + U/pow(-1.*(MassSq*MassSq) + U,2))) + (37.333333333333336*(MassSq*MassSq) + (5.333333333333333*(9.*pow(MassSq,4) - 3.*(MassSq*MassSq)*(2.*(MassSq*MassSq) + S) + (S + T)*(S + U)))/(-1.*(MassSq*MassSq) + U))/(-1.*(MassSq*MassSq) + T) - (48.*((5.*pow(MassSq,4) + MassSq*MassSq*U - 1.*T*(-1.*(MassSq*MassSq) + U) - 1.*(MassSq*MassSq)*(3.*S + 4.*T + 2.*U))/(-1.*(MassSq*MassSq) + U) + (5.*pow(MassSq,4) + MassSq*MassSq*U - 1.*T*(-1.*(MassSq*MassSq) + U) - 1.*(MassSq*MassSq)*(3.*S + 2.*T + 4.*U))/(-1.*(MassSq*MassSq) + T)))/S);

    return msquaredtree;
}
