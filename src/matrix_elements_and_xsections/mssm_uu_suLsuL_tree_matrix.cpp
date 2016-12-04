/*
 *    checked with MadGraph    
 */

double Process::matrixMSSMTree_uu_suLsuL( double S, double T ) {
   double alphaS = pdf->alphasQ( mu_r );
	double U = 2*MassSq*MassSq - S - T;
	/* factor of 0.5 introduced by hand, as mathematica output for left-left is actually twice as real left-left */
	/* because of double counting of final state particles */
   return 0.5*(105.27578027828648*(alphaS*alphaS)*(MassGlu*MassGlu)*S*
      (3./pow(MassGlu*MassGlu - 1.*T,2) + 3./pow(MassGlu*MassGlu - 1.*U,2) -
      2./((MassGlu*MassGlu - 1.*T)*(MassGlu*MassGlu - 1.*U))));
}  