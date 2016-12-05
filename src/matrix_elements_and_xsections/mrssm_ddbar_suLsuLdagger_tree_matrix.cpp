double Process::matrixMRSSMTree_ddbar_suLsuLdagger( double S, double T ) { // agrees with Philip
	double alphaS = pdf->alphasQ( mu_r );
	double U = 2*MassSq*MassSq - S - T;
    return  (315.82734083485946*(alphaS*alphaS)*(MassGlu*MassGlu)*S)/pow(-1.*(MassGlu*MassGlu) + T,2);
} 
