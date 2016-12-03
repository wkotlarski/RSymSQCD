inline double Process::matrixMSSMTree_ud_suLsdR(double alphaS, double T, double U, double S) {
      return (315.82734083485946*(alphaS*alphaS)*(-1.*pow(MassSq,4) + T*U))/pow(MassGlu*MassGlu - 1.*T,2); /*  ONLY HAVE OF THE VALUE WHICH IS COMMENTED: left-right + right-left or 2 of left-right*/
                         // + 2.*(315.82734083485946*(alphaS*alphaS)*(MassGlu*MassGlu)*S)/pow(MassGlu*MassGlu - 1.*T,2); /* left-left + right-right or 2 of left-left */
                         //
                         }
