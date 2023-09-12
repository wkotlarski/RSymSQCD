double MSSM::matrixHard_gu_suLsuLdaggeru(const double alphas, std::array<std::array<double, 4>, 5> const &p) const {
   const double k12 = p[0][0]*p[1][0]-p[0][3]*p[1][3];
   const double k13 = p[0][0]*p[2][0]-p[0][3]*p[2][3];
   const double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   const double T15 = -2*(p[0][0]*p[4][0]-p[0][3]*p[4][3]);
   const double T25 = -2*(p[1][0]*p[4][0]-p[1][3]*p[4][3]);
   const double S   = 2*k12;
   const double T   = pow<2>(MassSq) - 2*k13;
   const double S35 = pow<2>(MassSq) + 2*k35;
   const double MassGlu2 = pow<2>(MassGlu);
   const double MassSq2 = pow<2>(MassSq);
   const double res = -128/(9.*S) + 2176/(9.*(MassSq2 - S - T - T15)) - 224/T15 + (8*(12*MassSq2 - 24*S35 - 12*T))/(27.*S*T15) + (8*(24*S35 + 12*T15))/(27.*(MassSq2 - T)*(MassSq2 - S - T - T15)) + (8*(168*S35 + 84*T15))/(27.*S*(MassSq2 - T)) + (8*(-168*S35 + 144*T + 48*T15 - 48*T25))/(27.*S*(MassSq2 - S - T - T15)) + (8*(-432*MassSq2 - 432*S + 864*S35 - 864*T - 432*T15))/(27.*(T25*T25)) + (8*(-48*(S35*S35)*T + 96*S35*(T*T) - 48*pow<3>(T) - 48*S35*T*T15 + 48*(T*T)*T15))/(27.*(MassSq2 - T)*(MassSq2 - S - T - T15)*(T25*T25)) + (8*(-624*(S*S) + 864*S*S35 - 240*(S35*S35) - 1296*S*T + 960*S35*T - 720*(T*T) - 1632*S*T15 + 1248*S35*T15 - 2064*T*T15 - 1008*(T15*T15)))/(27.*(MassSq2 - S - T - T15)*(T25*T25)) + (8*(-192*(S35*S35) + 768*S35*T - 576*(T*T) - 384*S35*T15 + 768*T*T15 - 192*(T15*T15)))/(27.*(MassSq2 - T)*(T25*T25)) + (8*(-192*(S35*S35)*T + 384*S35*(T*T) - 192*pow<3>(T) - 384*S35*T*T15 + 384*(T*T)*T15 - 192*T*(T15*T15)))/(27.*pow<2>(MassSq2 - T)*(T25*T25)) + (8*(-192*pow<3>(S) + 384*(S*S)*S35 - 192*S*(S35*S35) - 576*(S*S)*T + 768*S*S35*T - 192*(S35*S35)*T - 576*S*(T*T) + 384*S35*(T*T) - 192*pow<3>(T) - 576*(S*S)*T15 + 768*S*S35*T15 - 192*(S35*S35)*T15 - 1152*S*T*T15 + 768*S35*T*T15 - 576*(T*T)*T15 - 576*S*(T15*T15) + 384*S35*(T15*T15) - 576*T*(T15*T15) - 192*pow<3>(T15)))/(27.*pow<2>(MassSq2 - S - T - T15)*(T25*T25)) - 1856/(9.*T25) + (8*(1164*MassSq2 - 1320*S35 - 60*T - 576*T15))/(27.*S*T25) + (8*(-144*S35 + 120*T - 120*T15))/(27.*(MassSq2 - T)*T25) + (8*(312*MassSq2 + 336*S - 708*S35 + 504*T))/(27.*T15*T25) + (8*(-24*(S35*S35) + 24*S35*T))/(27.*(MassSq2 - T)*T15*T25) + (8*(540*pow<2>(MassSq2) - 1188*MassSq2*S35 + 864*(S35*S35) + 108*MassSq2*T - 540*S35*T + 216*(T*T)))/(27.*S*T15*T25) + (8*(-48*pow<3>(S35) + 144*(S35*S35)*T - 144*S35*(T*T) + 48*pow<3>(T)))/(27.*S*(MassSq2 - T)*T15*T25) + (8*(336*(S*S) - 756*S*S35 + 588*(S35*S35) + 840*S*T - 1260*S35*T + 672*(T*T)))/(27.*(MassSq2 - S - T - T15)*T15*T25) + (8*(-168*pow<3>(S35) + 504*(S35*S35)*T - 504*S35*(T*T) + 168*pow<3>(T)))/(27.*S*(MassSq2 - S - T - T15)*T15*T25) + (8*(-408*S + 372*S35 - 1080*T + 24*T15))/(27.*(MassSq2 - S - T - T15)*T25) + (8*(-192*S*S35 - 192*S35*T - 192*S35*T15))/(27.*pow<2>(MassSq2 - S - T - T15)*T25) + (8*(-192*S35*T - 192*T*T15))/(27.*pow<2>(MassSq2 - T)*T25) + (8*(24*(S35*S35) - 96*S35*T + 24*(T*T) + 24*S35*T15 - 48*T*T15))/(27.*(MassSq2 - T)*(MassSq2 - S - T - T15)*T25) + (8*(-120*(S35*S35) + 264*S35*T - 144*(T*T) + 120*S35*T15 - 144*T*T15 - 48*(T15*T15)))/(27.*S*(MassSq2 - S - T - T15)*T25) + (8*(420*(S35*S35) - 756*S35*T + 336*(T*T) + 336*S35*T15 - 252*T*T15 + 84*(T15*T15)))/(27.*S*(MassSq2 - T)*T25) - (1088*MassSq2)/(9.*pow<2>(S + T15 + T25)) - (128*MassSq2*(T15*T15))/(T25*T25*pow<2>(S + T15 + T25)) - (128*MassSq2*T15)/(T25*pow<2>(S + T15 + T25)) - (256*MassSq2*T25)/(9.*S*pow<2>(S + T15 + T25)) - (256*MassSq2*T25)/(9.*T15*pow<2>(S + T15 + T25)) + 4448/(9.*(S + T15 + T25)) + (128*T)/(9.*(MassSq2 - T)*(S + T15 + T25)) + (8*(-84*(S35*S35) + 252*S35*T - 168*(T*T)))/(27.*S*(MassSq2 - T)*(S + T15 + T25)) + (8*(-24*(S35*S35) + 24*S35*T))/(27.*(MassSq2 - T)*T15*(S + T15 + T25)) + (8*(-600*pow<2>(MassSq2) + 1260*MassSq2*S35 - 888*(S35*S35) - 60*MassSq2*T + 516*S35*T - 228*(T*T)))/(27.*S*T15*(S + T15 + T25)) + (8*(48*pow<3>(S35) - 144*(S35*S35)*T + 144*S35*(T*T) - 48*pow<3>(T)))/(27.*S*(MassSq2 - T)*T15*(S + T15 + T25)) + (8*(168*pow<3>(S35) - 504*(S35*S35)*T + 504*S35*(T*T) - 168*pow<3>(T)))/(27.*S*(MassSq2 - S - T - T15)*T15*(S + T15 + T25)) + (8*(-1800*S35 + 2196*T - 972*T15 - 1620*T25))/(27.*(MassSq2 - S - T - T15)*(S + T15 + T25)) + (8*(-492*MassSq2 + 264*S35 + 24*T - 48*T25))/(27.*S*(S + T15 + T25)) + (8*(432*S35*T*T15 - 432*(T*T)*T15))/(27.*(MassSq2 - S - T - T15)*(T25*T25)*(S + T15 + T25)) + (8*(-864*MassSq2*T15 + 864*S35*T15 - 864*T*T15 + 432*(T15*T15)))/(27.*(T25*T25)*(S + T15 + T25)) + (8*(432*S35*T*T15 - 432*(T*T)*T15 + 432*T*(T15*T15)))/(27.*(MassSq2 - T)*(T25*T25)*(S + T15 + T25)) + (8*(-1620*pow<2>(MassSq2) + 3024*MassSq2*S35 - 1944*(S35*S35) + 216*MassSq2*T + 864*S35*T - 540*(T*T)))/(27.*S*T25*(S + T15 + T25)) + (8*(216*pow<3>(S35) - 648*(S35*S35)*T + 648*S35*(T*T) - 216*pow<3>(T)))/(27.*S*(MassSq2 - T)*T25*(S + T15 + T25)) + (8*(216*pow<3>(S35) - 648*(S35*S35)*T + 648*S35*(T*T) - 216*pow<3>(T)))/(27.*S*(MassSq2 - S - T - T15)*T25*(S + T15 + T25)) + (8*(-2160*MassSq2 + 2376*S35 - 864*T + 1284*T15))/(27.*T25*(S + T15 + T25)) + (8*(-216*(S35*S35) + 528*S35*T - 312*(T*T) + 312*T*T15))/(27.*(MassSq2 - T)*T25*(S + T15 + T25)) + (8*(-432*(S35*S35) + 1200*S35*T - 768*(T*T) - 648*S35*T15 + 1296*T*T15))/(27.*(MassSq2 - S - T - T15)*T25*(S + T15 + T25)) + (8*(-312*MassSq2 - 372*S + 444*S35 - 312*T - 348*T15 - 6*T25))/(27.*pow<2>(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) - 232/(9.*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-44*S35 + 42*T - 13*T15))/(27.*(MassSq2 - T)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-40*(S35*S35) + 40*S35*T))/(27.*(MassSq2 - T)*T15*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-64*S*S35 - 64*S35*T - 64*S35*T15))/(27.*pow<2>(MassSq2 - S - T - T15)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-344*MassSq2 - 163*S + 502*S35 - 266*T - 124*T15))/(27.*T25*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(90*pow<2>(MassSq2) - 28*MassSq2*S - 8*(S*S) - 198*MassSq2*S35 + 58*S*S35 + 144*(S35*S35) + 18*MassSq2*T - 12*S*T - 90*S35*T + 36*(T*T)))/(27.*T15*T25*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-20*S*(S35*S35) - 40*pow<3>(S35) + 20*S*S35*T + 120*(S35*S35)*T - 120*S35*(T*T) + 40*pow<3>(T)))/(27.*(MassSq2 - T)*T15*T25*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-8*pow<3>(S) + 18*(S*S)*S35 - 14*S*(S35*S35) + 4*pow<3>(S35) - 20*(S*S)*T + 30*S*S35*T - 12*(S35*S35)*T - 16*S*(T*T) + 12*S35*(T*T) - 4*pow<3>(T)))/(27.*(MassSq2 - S - T - T15)*T15*T25*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-8*(S35*S35)*T + 16*S35*(T*T) - 8*pow<3>(T) - 8*S35*T*T15 + 8*(T*T)*T15))/(27.*(MassSq2 - T)*(MassSq2 - S - T - T15)*T25*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-266*(S*S) + 364*S*S35 - 100*(S35*S35) - 532*S*T + 374*S35*T - 274*(T*T) - 572*S*T15 + 410*S35*T15 - 640*T*T15 - 314*(T15*T15)))/(27.*(MassSq2 - S - T - T15)*T25*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-31*S*S35 - 82*(S35*S35) + 11*S*T + 188*S35*T - 106*(T*T) - 11*S*T15 - 53*S35*T15 + 77*T*T15 - 11*(T15*T15)))/(27.*(MassSq2 - T)*T25*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-64*pow<3>(S) + 128*(S*S)*S35 - 64*S*(S35*S35) - 192*(S*S)*T + 256*S*S35*T - 64*(S35*S35)*T - 192*S*(T*T) + 128*S35*(T*T) - 64*pow<3>(T) - 192*(S*S)*T15 + 256*S*S35*T15 - 64*(S35*S35)*T15 - 384*S*T*T15 + 256*S35*T*T15 - 192*(T*T)*T15 - 192*S*(T15*T15) + 128*S35*(T15*T15) - 192*T*(T15*T15) - 64*pow<3>(T15)))/(27.*pow<2>(MassSq2 - S - T - T15)*T25*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(90*pow<2>(MassSq2) - 108*MassSq2*S35 + 36*(S35*S35) - 72*MassSq2*T + 36*S35*T + 18*(T*T) - 198*MassSq2*T15 + 108*S35*T15 - 54*T*T15))/(27.*T25*(S + T15 + T25)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(72*S35*T*T15 - 72*(T*T)*T15))/(27.*(MassSq2 - S - T - T15)*T25*(S + T15 + T25)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-12*MassSq2 + 2*S + 56*S35 + 4*T + 16*T25))/(27.*T15*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-30*MassSq2 + 74*S35 - 84*T - 42*T15 + 38*T25))/(27.*S*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-194*MassSq2 + 204*S35 - 48*T + 162*T15 + 126*T25))/(27.*(S + T15 + T25)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-152*S + 162*S35 - 302*T - 58*T15 + 184*T25))/(27.*(MassSq2 - S - T - T15)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-756*S + 1428*S35 - 1428*T + 840*T25))/(27.*(MassSq2 - S - T - T15)*T15) + (8*(-840*MassSq2 + 1548*S35 - 588*T + 840*T25))/(27.*T15*(S + T15 + T25)) + (8*(-10*pow<2>(MassSq2) + 12*MassSq2*S35 - 4*(S35*S35) + 8*MassSq2*T - 4*S35*T - 2*(T*T) + 2*MassSq2*T25 - 4*S35*T25 - 2*T*T25))/(27.*S*T15*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(4*(S35*S35) - 16*S35*T + 4*(T*T) + 4*S35*T15 - 8*T*T15 + 4*S35*T25 + 2*T15*T25))/(27.*(MassSq2 - T)*(MassSq2 - S - T - T15)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-486*(S*S) + 660*S*S35 - 180*(S35*S35) - 948*S*T + 642*S35*T - 462*(T*T) - 900*S*T15 + 606*S35*T15 - 888*T*T15 - 438*(T15*T15) + 36*S*T25 - 180*S35*T25 + 18*T*T25 + 6*T15*T25))/(27.*(MassSq2 - S - T - T15)*pow<2>(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-96*pow<3>(S) + 192*(S*S)*S35 - 96*S*(S35*S35) - 288*(S*S)*T + 384*S*S35*T - 96*(S35*S35)*T - 288*S*(T*T) + 192*S35*(T*T) - 96*pow<3>(T) - 288*(S*S)*T15 + 384*S*S35*T15 - 96*(S35*S35)*T15 - 576*S*T*T15 + 384*S35*T*T15 - 288*(T*T)*T15 - 288*S*(T15*T15) + 192*S35*(T15*T15) - 288*T*(T15*T15) - 96*pow<3>(T15) - 96*S*S35*T25 - 96*S35*T*T25 - 96*S35*T15*T25))/(27.*pow<2>(MassSq2 - S - T - T15)*pow<2>(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-924*(S35*S35) + 1764*S35*T - 840*(T*T) - 1596*S35*T25 + 1428*T*T25 - 840*(T25*T25)))/(27.*(MassSq2 - S - T - T15)*T15*(S + T15 + T25)) + (8*(-72*(S35*S35) + 136*S35*T - 64*(T*T) - 108*S35*T15 + 216*T*T15 - 204*S35*T25 + 174*T*T25 - 162*T15*T25 - 126*(T25*T25)))/(27.*(MassSq2 - S - T - T15)*(S + T15 + T25)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-82*pow<2>(MassSq2) + 202*MassSq2*S35 - 160*(S35*S35) - 38*MassSq2*T + 118*S35*T - 40*(T*T) + 40*MassSq2*T25 - 140*S35*T25 + 80*T*T25 - 40*(T25*T25)))/(27.*S*(S + T15 + T25)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(-100*(S35*S35) + 220*S35*T - 120*(T*T) + 100*S35*T15 - 120*T*T15 - 40*(T15*T15) - 140*S35*T25 + 120*T*T25 + 40*T15*T25 - 40*(T25*T25)))/(27.*S*(MassSq2 - S - T - T15)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(16*pow<2>(MassSq2) + 16*MassSq2*S35 - 20*(S35*S35) - 48*MassSq2*T + 24*S35*T + 12*(T*T) - 12*MassSq2*T25 - 38*S35*T25 + 14*T*T25 - 20*(T25*T25)))/(27.*T15*(S + T15 + T25)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(2*(S*S) + 2*S*S35 - 6*(S35*S35) - 6*S*T + 18*S35*T - 12*(T*T) + 16*S*T25 - 30*S35*T25 + 34*T*T25 - 20*(T25*T25)))/(27.*(MassSq2 - S - T - T15)*T15*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(48*pow<2>(MassSq2) - 18*MassSq2*S - 24*(S*S) + 48*MassSq2*S35 + 48*S*S35 - 60*(S35*S35) - 144*MassSq2*T + 12*S*T + 72*S35*T + 36*(T*T) - 6*MassSq2*T25 + 30*S*T25 - 66*S35*T25 + 54*T*T25 - 6*(T25*T25)))/(27.*T15*pow<2>(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(168*(S35*S35) - 312*S35*T + 144*(T*T) + 168*S35*T25 - 144*T*T25 + 48*(T25*T25)))/(27.*S*(MassSq2 - S - T - T15)*(S + T15 + T25)) + (8*(-24*pow<3>(S) + 54*(S*S)*S35 - 42*S*(S35*S35) + 12*pow<3>(S35) - 60*(S*S)*T + 90*S*S35*T - 36*(S35*S35)*T - 48*S*(T*T) + 36*S35*(T*T) - 12*pow<3>(T) + 30*(S*S)*T25 - 48*S*S35*T25 + 24*(S35*S35)*T25 + 42*S*T*T25 - 36*S35*T*T25 + 12*(T*T)*T25 - 6*S*(T25*T25) + 12*S35*(T25*T25)))/(27.*(MassSq2 - S - T - T15)*T15*pow<2>(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(4*pow<3>(S35) - 12*(S35*S35)*T + 12*S35*(T*T) - 4*pow<3>(T) + 22*(S35*S35)*T25 - 42*S35*T*T25 + 20*(T*T)*T25 + 38*S35*(T25*T25) - 34*T*(T25*T25) + 20*pow<3>(T25)))/(27.*(MassSq2 - S - T - T15)*T15*(S + T15 + T25)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25)) + (8*(40*pow<3>(S35) - 120*(S35*S35)*T + 120*S35*(T*T) - 40*pow<3>(T) + 140*(S35*S35)*T25 - 260*S35*T*T25 + 120*(T*T)*T25 + 140*S35*(T25*T25) - 120*T*(T25*T25) + 40*pow<3>(T25)))/(27.*S*(MassSq2 - S - T - T15)*(S + T15 + T25)*(MassGlu2 - 3*MassSq2 + S + S35 + T + T15 + T25));
   return pow<3>(alphas*pi)*res;
}