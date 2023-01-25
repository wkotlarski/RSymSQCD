(* ----------------------------------------------------------------------------- *) 
(* This model file was automatically created by SARAH version4.5.8  *) 
(* SARAH References: arXiv:0806.0538, 0909.2863, 1002.0840, 1207.0906, 1309.7223 *) 
(* (c) Florian Staub, 2013  *) 
(* ----------------------------------------------------------------------------- *) 
(* File created at 18:25 on 2.12.2015  *) 
(* ---------------------------------------------------------------------- *) 
 
 
IndexRange[  Index[Colour]  ] =NoUnfold[Range[3]]; 
IndexStyle[  Index[Colour, i_Integer ] ] := Greek[i];  
IndexRange[  Index[I3Gen]  ] =Range[3]; 
IndexStyle[  Index[I3Gen, i_Integer ] ] := Alph[ 8+i];  
IndexRange[  Index[I6Gen]  ] =Range[6]; 
IndexStyle[  Index[I6Gen, i_Integer ] ] := Alph[ 8+i];  
IndexRange[  Index[Gluon]  ] =NoUnfold[Range[8]]; 
IndexStyle[  Index[Gluon, i_Integer ] ] := Alph[ 8+i];  

 
 
(* Definitions for trigonometric functions  
Sin[ThetaW]: STW
Sin[2*ThetaW]: S2TW
Cos[ThetaW]: CTW
Cos[2*ThetaW]: C2TW
*) 
 
Conjugate[STW] ^= STW
Conjugate[S2TW] ^= S2TW
Conjugate[CTW] ^= CTW
Conjugate[C2TW] ^= C2TW
 
 
Lam[a_,b_,c_]:=2*SUNT[a,b,c]; 
fSU3[a_,b_,c_]:=SUNF[a,b,c]; 
LambdaProd[a_,b_][c_,d_]:=4*SUNT[a,b,c,d]; 
 
 
M$ClassesDescription= {
S[26] == {SelfConjugate -> True,
Indices -> {Index[Gluon]},
Mass -> MasssigmaO,
PropagatorLabel->ComposedChar["O","P"],
PropagatorType -> ScalarDash,
PropagatorArrow -> None},

 
S[25] == {SelfConjugate -> True,
Indices -> {Index[Gluon]},
Mass -> MassphiO,
PropagatorLabel->ComposedChar["O","S"],
PropagatorType -> ScalarDash,
PropagatorArrow -> None},

 
S[14] == {SelfConjugate -> False,
Indices -> {Index[I6Gen], Index[Colour]},
Mass -> MassSd,
PropagatorLabel->ComposedChar["d",Index[I6Gen],Null,"\\tilde"],
PropagatorType -> ScalarDash,
PropagatorArrow -> Forward},

 
S[13] == {SelfConjugate -> False,
Indices -> {Index[I6Gen], Index[Colour]},
Mass -> MassSu,
PropagatorLabel->ComposedChar["u",Index[I6Gen],Null,"\\tilde"],
PropagatorType -> ScalarDash,
PropagatorArrow -> Forward},

 
F[15] == {SelfConjugate -> False,
Indices -> {Index[Gluon]},
Mass -> MassGlu,
PropagatorLabel->ComposedChar["g", Null, Null ,"\\tilde"],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
F[4] == {SelfConjugate -> False,
Indices -> {Index[I3Gen], Index[Colour]},
Mass -> MassFd,
PropagatorLabel->ComposedChar["d",Index[I3Gen]],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
F[3] == {SelfConjugate -> False,
Indices -> {Index[I3Gen], Index[Colour]},
Mass -> MassFu,
PropagatorLabel->ComposedChar["u",Index[I3Gen]],
PropagatorType -> Straight,
PropagatorArrow -> Forward},

 
V[5] == {SelfConjugate -> True,
Indices -> {Index[Gluon]},
Mass -> 0,
PropagatorLabel->ComposedChar["g"],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[1] == {SelfConjugate -> True,
Indices -> {},
Mass -> 0,
PropagatorLabel->ComposedChar["\\gamma"],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[2] == {SelfConjugate -> True,
Indices -> {},
Mass -> MassVZ,
PropagatorLabel->ComposedChar["Z"],
PropagatorType -> Sine,
PropagatorArrow -> None},

 
V[3] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassVWm,
PropagatorLabel->ComposedChar["W","-"],
PropagatorType -> Sine,
PropagatorArrow -> Forward},

 
U[5] == {SelfConjugate -> False,
Indices -> {Index[Gluon]},
Mass -> 0,
PropagatorLabel->ComposedChar["\\eta",Index[Gluon],"G"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[1] == {SelfConjugate -> False,
Indices -> {},
Mass -> 0,
PropagatorLabel->ComposedChar["\\eta","\\gamma"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[2] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassVZ,
PropagatorLabel->ComposedChar["\\eta","Z"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[3] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassVWm,
PropagatorLabel->ComposedChar["\\eta","-"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward},

 
U[4] == {SelfConjugate -> False,
Indices -> {},
Mass -> MassVWm,
PropagatorLabel->ComposedChar["\\eta","+"],
PropagatorType -> GhostDash,
PropagatorArrow -> Forward}
 
}

 
MasssigmaO[gen_, y_] = MasssigmaO
MassphiO[gen_, y_] = MassphiO
MassSd[gen_, y_] = MassSd[gen]
MassSu[gen_, y_] = MassSu[gen]
MassGlu[gen_, y_] = MassGlu
MassFd[gen_, y_] = MassFd[gen]
MassFu[gen_, y_] = MassFu[gen]


GaugeXi[S[26,___]] = 1 
GaugeXi[S[25,___]] = 1 
GaugeXi[S[14,___]] = 1 
GaugeXi[S[13,___]] = 1 


GaugeXi[V[5,___]] = GaugeXi[G]
GaugeXi[V[1,___]] = GaugeXi[P]
GaugeXi[V[2,___]] = GaugeXi[Z]
GaugeXi[V[3,___]] = GaugeXi[Wm]


M$CouplingMatrices= {

 C[V[5, {ct1}], V[5, {ct2}], V[5, {ct3}], V[5, {ct4}]] == 
  {{(-I)*g3^2*(SUNF[ct1, ct3, ct2, ct4] - SUNF[ct1, ct4, ct3, ct2]), 
    (-2*I)*(dZGG1 + dZgs1)*g3^2*(SUNF[ct1, ct3, ct2, ct4] - 
      SUNF[ct1, ct4, ct3, ct2])}, 
   {(-I)*g3^2*(SUNF[ct1, ct2, ct3, ct4] + SUNF[ct1, ct4, ct3, ct2]), 
    (-2*I)*(dZGG1 + dZgs1)*g3^2*(SUNF[ct1, ct2, ct3, ct4] + 
      SUNF[ct1, ct4, ct3, ct2])}, 
   {I*g3^2*(SUNF[ct1, ct2, ct3, ct4] + SUNF[ct1, ct3, ct2, ct4]), 
    (2*I)*(dZGG1 + dZgs1)*g3^2*(SUNF[ct1, ct2, ct3, ct4] + 
      SUNF[ct1, ct3, ct2, ct4])}}, 
C[V[5, {ct1}], V[5, {ct2}], V[5, {ct3}]] == 
  {{g3*SUNF[ct1, ct2, ct3], ((3*dZGG1 + 2*dZgs1)*g3*SUNF[ct1, ct2, ct3])/2}}, 
(*taking gluon self interaction just from mssmct*)
(* C[V[5, {ct1}], V[5, {ct2}], V[5, {ct3}], V[5, {ct4}]] == { *)
(* {I*g3^2*(-(fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3]) - fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] - fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] - fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] - fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] - fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] - fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] - fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] - fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] - fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] - fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] - fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] - fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] - fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] - fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3] - fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4]),  *)
(* 2*I*(dZGG1 + dZgs1) g3^2 *(-(fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3]) - fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] - fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] - fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] - fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] - fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] - fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] - fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] - fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] - fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] - fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] - fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] - fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] - fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] - fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3] - fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4])},  *)
(* {I*g3^2*(fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3] - fSU3[1, ct1, ct2]*fSU3[1, ct3, ct4] + fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] - fSU3[2, ct1, ct2]*fSU3[2, ct3, ct4] + fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] - fSU3[3, ct1, ct2]*fSU3[3, ct3, ct4] + fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] - fSU3[4, ct1, ct2]*fSU3[4, ct3, ct4] + fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] - fSU3[5, ct1, ct2]*fSU3[5, ct3, ct4] + fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] - fSU3[6, ct1, ct2]*fSU3[6, ct3, ct4] + fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] - fSU3[7, ct1, ct2]*fSU3[7, ct3, ct4] + fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3] - fSU3[8, ct1, ct2]*fSU3[8, ct3, ct4]),  *)
(* 2 I (dZGG1 + dZgs1) g3^2 (fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3] - fSU3[1, ct1, ct2]*fSU3[1, ct3, ct4] + fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] - fSU3[2, ct1, ct2]*fSU3[2, ct3, ct4] + fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] - fSU3[3, ct1, ct2]*fSU3[3, ct3, ct4] + fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] - fSU3[4, ct1, ct2]*fSU3[4, ct3, ct4] + fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] - fSU3[5, ct1, ct2]*fSU3[5, ct3, ct4] + fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] - fSU3[6, ct1, ct2]*fSU3[6, ct3, ct4] + fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] - fSU3[7, ct1, ct2]*fSU3[7, ct3, ct4] + fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3] - fSU3[8, ct1, ct2]*fSU3[8, ct3, ct4])},  *)
(* {I*g3^2*(fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] + fSU3[1, ct1, ct2]*fSU3[1, ct3, ct4] + fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] + fSU3[2, ct1, ct2]*fSU3[2, ct3, ct4] + fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] + fSU3[3, ct1, ct2]*fSU3[3, ct3, ct4] + fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] + fSU3[4, ct1, ct2]*fSU3[4, ct3, ct4] + fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] + fSU3[5, ct1, ct2]*fSU3[5, ct3, ct4] + fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] + fSU3[6, ct1, ct2]*fSU3[6, ct3, ct4] + fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] + fSU3[7, ct1, ct2]*fSU3[7, ct3, ct4] + fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4] + fSU3[8, ct1, ct2]*fSU3[8, ct3, ct4]),  *)
(* 2 I (dZGG1 + dZgs1) g3^2*(fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] + fSU3[1, ct1, ct2]*fSU3[1, ct3, ct4] + fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] + fSU3[2, ct1, ct2]*fSU3[2, ct3, ct4] + fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] + fSU3[3, ct1, ct2]*fSU3[3, ct3, ct4] + fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] + fSU3[4, ct1, ct2]*fSU3[4, ct3, ct4] + fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] + fSU3[5, ct1, ct2]*fSU3[5, ct3, ct4] + fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] + fSU3[6, ct1, ct2]*fSU3[6, ct3, ct4] + fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] + fSU3[7, ct1, ct2]*fSU3[7, ct3, ct4] + fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4] + fSU3[8, ct1, ct2]*fSU3[8, ct3, ct4]) *)
(* }}, *)


(*  C[V[5, {ct1}], V[5, {ct2}], V[5, {ct3}]] == {{g3*fSU3[ct1, ct2, ct3],g3*(dZgs1+3*dZGG1/2)*fSU3[ct1, ct2, ct3]}}, *)

 (* C[-V[3], V[1], V[1], V[3]] == {{I*g2^2*STW^2}, {I*g2^2*STW^2}, {(-2*I)*g2^2*STW^2}}, *)
 (* C[-V[3], V[1], V[3], V[2]] == {{(I/2)*g2^2*S2TW}, {(-I)*g2^2*S2TW}, {(I/2)*g2^2*S2TW}}, *)
 (* C[-V[3], -V[3], V[3], V[3]] == {{(2*I)*g2^2}, {(-I)*g2^2}, {(-I)*g2^2}}, *)
 (* C[-V[3], V[3], V[2], V[2]] == {{(-2*I)*CTW^2*g2^2}, {I*CTW^2*g2^2}, {I*CTW^2*g2^2}}, *)
 C[S[25, {ct1}], S[25, {ct2}], V[5, {ct3}], V[5, {ct4}]] == {{
I*g3^2*(SUNF[ct1, ct3, ct2, ct4] + SUNF[ct1, ct4, ct2, ct3])
(*(fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3] + fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] + fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] + fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] + fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] + fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] + fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] + fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3]) + I*g3^2*(fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] + fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] + fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] + fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] + fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] + fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] + fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] + fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4])*),
I*(dZPOc1+dZGG1 +2*dZgs1)*g3^2*(SUNF[ct1, ct3, ct2, ct4] + SUNF[ct1, ct4, ct2, ct3])
(*(fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3] + fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] + fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] + fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] + fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] + fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] + fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] + fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3]) + I*g3^2*(fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] + fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] + fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] + fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] + fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] + fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] + fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] + fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4])*)}},
 
C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[5, {ct3}], V[5, {ct4}]] == {{
(I/4)*g3^2*IndexDelta[gt1, gt2]*LambdaProd[ct3,ct4][ct2,ct1] (*Lam[ct3, ct2, 1]*Lam[ct4, 1, ct1] + Lam[ct3, ct2, 2]*Lam[ct4, 2, ct1] + Lam[ct3, ct2, 3]*Lam[ct4, 3, ct1])*) 
+ (I/4)*g3^2*IndexDelta[gt1, gt2]*LambdaProd[ct4,ct3][ct2,ct1] (*(Lam[ct3, 1, ct1]*Lam[ct4, ct2, 1] + Lam[ct3, 2, ct1]*Lam[ct4, ct2, 2] + Lam[ct3, 3, ct1]*Lam[ct4, ct2, 3])*),
(I/4)*g3^2*(dZGG1+2*dZgs1+dZSd[gt1]/2+dZSd[gt2]/2)*IndexDelta[gt1, gt2]*LambdaProd[ct3,ct4][ct2,ct1] + 
(I/4)*g3^2*(dZGG1+2*dZgs1+dZSd[gt1]/2+dZSd[gt2]/2)*IndexDelta[gt1, gt2]*LambdaProd[ct4,ct3][ct2,ct1] }},

 (* C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[5, {ct3}], V[1]] == {{(I/6)*CTW*g1*g3*Lam[ct3, ct2, ct1]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) - (I/2)*g2*g3*STW*Lam[ct3, ct2, ct1]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) - (I/3)*CTW*g1*g3*Lam[ct3, ct2, ct1]*(Conjugate[ZD[gt1, 4]]*ZD[gt2, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt2, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt2, 6])}}, *)
 (* C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[5, {ct3}], V[2]] == {{(-I/2)*CTW*g2*g3*Lam[ct3, ct2, ct1]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) - (I/6)*g1*g3*STW*Lam[ct3, ct2, ct1]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) + (I/3)*g1*g3*STW*Lam[ct3, ct2, ct1]*(Conjugate[ZD[gt1, 4]]*ZD[gt2, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt2, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt2, 6])}}, *)
 (* C[S[14, {gt1, ct1}], -S[13, {gt2, ct2}], -V[3], V[5, {ct4}]] == {{(I*g2*g3*Lam[ct4, ct2, ct1]*(Conjugate[ZD[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZU[gt2, 3]))/Sqrt[2]}}, *)

 (* C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[1], V[1]] == {{(I/18)*CTW^2*g1^2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) - (I/3)*CTW*g1*g2*STW*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) + (I/2)*g2^2*STW^2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) + ((2*I)/9)*CTW^2*g1^2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 4]]*ZD[gt2, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt2, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt2, 6])}}, *)
 (* C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[1], V[2]] == {{(-I/6)*C2TW*g1*g2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) - (I/36)*g1^2*S2TW*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) + (I/4)*g2^2*S2TW*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) - (I/9)*g1^2*S2TW*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 4]]*ZD[gt2, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt2, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt2, 6])}}, *)
 (* C[S[14, {gt1, ct1}], -S[13, {gt2, ct2}], -V[3], V[1]] == {{((I/3)*CTW*g1*g2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZU[gt2, 3]))/Sqrt[2]}}, *)
 (* C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], -V[3], V[3]] == {{(I/2)*g2^2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3])}}, *)
 (* C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[2], V[2]] == {{(I/2)*CTW^2*g2^2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) + (I/3)*CTW*g1*g2*STW*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) + (I/18)*g1^2*STW^2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) + ((2*I)/9)*g1^2*STW^2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 4]]*ZD[gt2, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt2, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt2, 6])}}, *)
 (* C[S[14, {gt1, ct1}], -S[13, {gt2, ct2}], -V[3], V[2]] == {{((-I/3)*g1*g2*STW*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZU[gt2, 3]))/Sqrt[2]}}, *)

 C[S[26, {ct1}], S[26, {ct2}], V[5, {ct3}], V[5, {ct4}]] == {{
I*g3^2*(SUNF[ct1, ct3, ct2, ct4] + SUNF[ct1, ct4, ct2, ct3])
(**(fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3] + fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] + fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] + fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] + fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] + fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] + fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] + fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3]) + I*g3^2*(fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] + fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] + fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] + fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] + fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] + fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] + fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] + fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4])*),
I*(dZSOc1+dZGG1 +2*dZgs1)*g3^2*(SUNF[ct1, ct3, ct2, ct4] + SUNF[ct1, ct4, ct2, ct3])}},

 C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[5, {ct3}], V[5, {ct4}]] == {{
(I/4)*g3^2*IndexDelta[gt1, gt2]*LambdaProd[ct3,ct4][ct2,ct1] (**(Lam[ct3, ct2, 1]*Lam[ct4, 1, ct1] + Lam[ct3, ct2, 2]*Lam[ct4, 2, ct1] + Lam[ct3, ct2, 3]*Lam[ct4, 3, ct1])*) + 
(I/4)*g3^2*IndexDelta[gt1, gt2]*LambdaProd[ct4,ct3][ct2,ct1] (**(Lam[ct3, 1, ct1]*Lam[ct4, ct2, 1] + Lam[ct3, 2, ct1]*Lam[ct4, ct2, 2] + Lam[ct3, 3, ct1]*Lam[ct4, ct2, 3])*),
(I/4)*g3^2*(dZGG1+2*dZgs1+dZSu1[gt1]/2+dZSu1[gt2]/2)*IndexDelta[gt1, gt2]*LambdaProd[ct3,ct4][ct2,ct1] + (I/4)*g3^2*(dZGG1+2*dZgs1+dZSu1[gt1]/2+dZSu1[gt2]/2)*IndexDelta[gt1, gt2]*LambdaProd[ct4,ct3][ct2,ct1]}},

 (* C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[5, {ct3}], V[1]] == {{(I/6)*CTW*g1*g3*Lam[ct3, ct2, ct1]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + (I/2)*g2*g3*STW*Lam[ct3, ct2, ct1]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + ((2*I)/3)*CTW*g1*g3*Lam[ct3, ct2, ct1]*(Conjugate[ZU[gt1, 4]]*ZU[gt2, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt2, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt2, 6])}}, *)
 (* C[S[13, {gt1, ct1}], -S[14, {gt2, ct2}], V[5, {ct3}], V[3]] == {{(I*g2*g3*Lam[ct3, ct2, ct1]*(Conjugate[ZU[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZD[gt2, 3]))/Sqrt[2]}}, *)
 (* C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[5, {ct3}], V[2]] == {{(I/2)*CTW*g2*g3*Lam[ct3, ct2, ct1]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) - (I/6)*g1*g3*STW*Lam[ct3, ct2, ct1]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) - ((2*I)/3)*g1*g3*STW*Lam[ct3, ct2, ct1]*(Conjugate[ZU[gt1, 4]]*ZU[gt2, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt2, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt2, 6])}}, *)
 (* C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[1], V[1]] == {{(I/18)*CTW^2*g1^2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + (I/3)*CTW*g1*g2*STW*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + (I/2)*g2^2*STW^2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + ((8*I)/9)*CTW^2*g1^2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 4]]*ZU[gt2, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt2, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt2, 6])}}, *)
 (* C[S[13, {gt1, ct1}], -S[14, {gt2, ct2}], V[1], V[3]] == {{((I/3)*CTW*g1*g2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZD[gt2, 3]))/Sqrt[2]}}, *)
 (* C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[1], V[2]] == {{(I/6)*C2TW*g1*g2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) - (I/36)*g1^2*S2TW*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + (I/4)*g2^2*S2TW*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) - ((4*I)/9)*g1^2*S2TW*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 4]]*ZU[gt2, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt2, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt2, 6])}}, *)
 (* C[S[13, {gt1, ct1}], -S[14, {gt2, ct2}], V[3], V[2]] == {{((-I/3)*g1*g2*STW*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZD[gt2, 3]))/Sqrt[2]}}, *)
 (* C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], -V[3], V[3]] == {{(I/2)*g2^2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3])}}, *)
 (* C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[2], V[2]] == {{(I/2)*CTW^2*g2^2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) - (I/3)*CTW*g1*g2*STW*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + (I/18)*g1^2*STW^2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + ((8*I)/9)*g1^2*STW^2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 4]]*ZU[gt2, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt2, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt2, 6])}}, *)

 C[F[4, {gt1, ct1}], F[15, {ct2}], -S[14, {gt3, ct3}]] == {{
((-I)*g3*Lam[ct2, ct3, ct1]*(Conjugate[ZDL[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZDL[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZDL[gt1, 3]]*ZD[gt3, 3]))/Sqrt[2],
((-I)*(dZFd1[gt1]/2+dZSd1[gt3]/2+dZGlL1/2+dZgs1+dZgs1restore)*g3*Lam[ct2, ct3, ct1]*(Conjugate[ZDL[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZDL[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZDL[gt1, 3]]*ZD[gt3, 3]))/Sqrt[2]
}, {0,0}},

 C[-F[15, {ct1}], F[4, {gt2, ct2}], -S[14, {gt3, ct3}]] == {{0,0}, {
(I*g3*Lam[ct1, ct3, ct2]*(ZD[gt3, 4]*ZDR[gt2, 1] + ZD[gt3, 5]*ZDR[gt2, 2] + ZD[gt3, 6]*ZDR[gt2, 3]))/Sqrt[2],
(I*(dZFd1[gt2]/2+dZSd1[gt3]/2+dZGlL1/2+dZgs1+dZgs1restore)*g3*Lam[ct1, ct3, ct2]*(ZD[gt3, 4]*ZDR[gt2, 1] + ZD[gt3, 5]*ZDR[gt2, 2] + ZD[gt3, 6]*ZDR[gt2, 3]))/Sqrt[2]
}},

 C[F[3, {gt1, ct1}], F[15, {ct2}], -S[13, {gt3, ct3}]] == {{
((-I)*g3*Lam[ct2, ct3, ct1]*(Conjugate[ZUL[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZUL[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZUL[gt1, 3]]*ZU[gt3, 3]))/Sqrt[2],
((-I)*(dZFu1[gt1]/2+dZSu1[gt3]/2+dZGlL1/2+dZgs1+dZgs1restore)*g3*Lam[ct2, ct3, ct1]*(Conjugate[ZUL[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZUL[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZUL[gt1, 3]]*ZU[gt3, 3]))/Sqrt[2]
}, {0,0}},

 C[-F[15, {ct1}], F[3, {gt2, ct2}], -S[13, {gt3, ct3}]] == {{0,0}, {
(I*g3*Lam[ct1, ct3, ct2]*(ZU[gt3, 4]*ZUR[gt2, 1] + ZU[gt3, 5]*ZUR[gt2, 2] + ZU[gt3, 6]*ZUR[gt2, 3]))/Sqrt[2],
(I*(dZFu1[gt2]/2+dZSu1[gt3]/2+dZGlL1/2+dZgs1+dZgs1restore)*g3*Lam[ct1, ct3, ct2]*(ZU[gt3, 4]*ZUR[gt2, 1] + ZU[gt3, 5]*ZUR[gt2, 2] + ZU[gt3, 6]*ZUR[gt2, 3]))/Sqrt[2]}},
 
C[-F[15, {ct1}], F[15, {ct2}], S[25, {ct3}]] == {{-(g3*fSU3[ct1, ct2, ct3]),-(dZPOc1+dZGlL1/2+dZGlR1/2+dZgs1)*(g3*fSU3[ct1, ct2, ct3])}, {-(g3*fSU3[ct1, ct2, ct3]),-(dZPOc1+dZGlL1/2+dZGlR1/2+dZgs1)*(g3*fSU3[ct1, ct2, ct3])}},
 
C[-F[4, {gt1, ct1}], F[15, {ct2}], S[14, {gt3, ct3}]] == {{
(I*g3*(Conjugate[ZD[gt3, 4]]*Conjugate[ZDR[gt1, 1]] + Conjugate[ZD[gt3, 5]]*Conjugate[ZDR[gt1, 2]] + Conjugate[ZD[gt3, 6]]*Conjugate[ZDR[gt1, 3]])*Lam[ct2, ct1, ct3])/Sqrt[2],
(I*(dZFd1[gt1]/2+dZSd1[gt3]/2+dZGlL1/2+dZgs1+dZgs1restore)*g3*(Conjugate[ZD[gt3, 4]]*Conjugate[ZDR[gt1, 1]] + Conjugate[ZD[gt3, 5]]*Conjugate[ZDR[gt1, 2]] + Conjugate[ZD[gt3, 6]]*Conjugate[ZDR[gt1, 3]])*Lam[ct2, ct1, ct3])/Sqrt[2]}, {0,0}},

 C[-F[15, {ct1}], F[15, {ct2}], S[26, {ct3}]] == {{I*g3*fSU3[ct1, ct2, ct3],I*(dZSOc1+dZGlL1/2+dZGlR1/2+dZgs1)*g3*fSU3[ct1, ct2, ct3]}, {(-I)*g3*fSU3[ct1, ct2, ct3],(-I)*(dZSOc1+dZGlL1/2+dZGlR1/2+dZgs1)*g3*fSU3[ct1, ct2, ct3]}},

 C[-F[3, {gt1, ct1}], F[15, {ct2}], S[13, {gt3, ct3}]] == {{
(I*g3*(Conjugate[ZU[gt3, 4]]*Conjugate[ZUR[gt1, 1]] + Conjugate[ZU[gt3, 5]]*Conjugate[ZUR[gt1, 2]] + Conjugate[ZU[gt3, 6]]*Conjugate[ZUR[gt1, 3]])*Lam[ct2, ct1, ct3])/Sqrt[2],
(I*(dZFu1[gt1]/2+dZSu1[gt3]/2+dZGlL1/2+dZgs1+dZgs1restore)*g3*(Conjugate[ZU[gt3, 4]]*Conjugate[ZUR[gt1, 1]] + Conjugate[ZU[gt3, 5]]*Conjugate[ZUR[gt1, 2]] + Conjugate[ZU[gt3, 6]]*Conjugate[ZUR[gt1, 3]])*Lam[ct2, ct1, ct3])/Sqrt[2]}, {0,0}},

 C[-F[4, {gt1, ct1}], -F[15, {ct2}], S[14, {gt3, ct3}]] == {{0,0}, {
((-I)*g3*Lam[ct2, ct1, ct3]*(Conjugate[ZD[gt3, 1]]*ZDL[gt1, 1] + Conjugate[ZD[gt3, 2]]*ZDL[gt1, 2] + Conjugate[ZD[gt3, 3]]*ZDL[gt1, 3]))/Sqrt[2],
((-I)*(dZFd1[gt1]/2+dZSd1[gt3]/2+dZGlL1/2+dZgs1+dZgs1restore)*g3*Lam[ct2, ct1, ct3]*(Conjugate[ZD[gt3, 1]]*ZDL[gt1, 1] + Conjugate[ZD[gt3, 2]]*ZDL[gt1, 2] + Conjugate[ZD[gt3, 3]]*ZDL[gt1, 3]))/Sqrt[2]}},

 C[-F[3, {gt1, ct1}], -F[15, {ct2}], S[13, {gt3, ct3}]] == {{0,0}, {
((-I)*g3*Lam[ct2, ct1, ct3]*(Conjugate[ZU[gt3, 1]]*ZUL[gt1, 1] + Conjugate[ZU[gt3, 2]]*ZUL[gt1, 2] + Conjugate[ZU[gt3, 3]]*ZUL[gt1, 3]))/Sqrt[2],
((-I)*(dZFu1[gt1]/2+dZSu1[gt3]/2+dZGlL1/2+dZgs1+dZgs1restore)*g3*Lam[ct2, ct1, ct3]*(Conjugate[ZU[gt3, 1]]*ZUL[gt1, 1] + Conjugate[ZU[gt3, 2]]*ZUL[gt1, 2] + Conjugate[ZU[gt3, 3]]*ZUL[gt1, 3]))/Sqrt[2]}},

 C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}], V[5, {ct3}]] == {{(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2],(-I/2)*(dZFd1[gt1]/2+dZFd1[gt2]/2+dZGG1/2+dZgs1)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2]}, 
{(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2],(-I/2)*(dZFd1[gt1]/2+dZFd1[gt2]/2+dZGG1/2+dZgs1)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2]}},

 (* C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}], V[1]] == {{(-I/6)*(CTW*g1 - 3*g2*STW)*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}, {(I/3)*CTW*g1*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}}, *)
 (* C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}], V[2]] == {{(I/6)*(3*CTW*g2 + g1*STW)*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}, {(-I/3)*g1*STW*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}}, *)
 (* C[-F[3, {gt1, ct1}], F[4, {gt2, ct2}], -V[3]] == {{((-I)*g2*IndexDelta[ct1, ct2]*(Conjugate[ZDL[gt2, 1]]*ZUL[gt1, 1] + Conjugate[ZDL[gt2, 2]]*ZUL[gt1, 2] + Conjugate[ZDL[gt2, 3]]*ZUL[gt1, 3]))/Sqrt[2]}, {0}}, *)

 C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}], V[5, {ct3}]] == {{(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2],(-I/2)*(dZFu1[gt1]/2+dZFu1[gt2]/2+dZGG1/2+dZgs1)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2]}, 
{(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2],(-I/2)*(dZFu1[gt1]/2+dZFu1[gt2]/2+dZGG1/2+dZgs1)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct1, ct2]}},

 (* C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}], V[1]] == {{(-I/6)*(CTW*g1 + 3*g2*STW)*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}, {((-2*I)/3)*CTW*g1*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}}, *)
 (* C[-F[4, {gt1, ct1}], F[3, {gt2, ct2}], V[3]] == {{((-I)*g2*IndexDelta[ct1, ct2]*(Conjugate[ZUL[gt2, 1]]*ZDL[gt1, 1] + Conjugate[ZUL[gt2, 2]]*ZDL[gt1, 2] + Conjugate[ZUL[gt2, 3]]*ZDL[gt1, 3]))/Sqrt[2]}, {0}}, *)
 (* C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}], V[2]] == {{(-I/6)*(3*CTW*g2 - g1*STW)*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}, {((2*I)/3)*g1*STW*IndexDelta[ct1, ct2]*IndexDelta[gt1, gt2]}}, *)

 C[-F[15, {ct1}], F[15, {ct2}], V[5, {ct3}]] == {{-(g3*fSU3[ct1, ct2, ct3]),-(g3*(dZGlL1+dZGG1/2+dZgs1)*fSU3[ct1, ct2, ct3])}, 
{-(g3*fSU3[ct1, ct2, ct3]),-(g3*(dZGlR1+dZGG1/2+dZgs1)*fSU3[ct1, ct2, ct3])}},

 C[S[25, {ct1}], S[14, {gt2, ct2}], -S[14, {gt3, ct3}]] == {{
(-I/2)*g3*(MDO + Conjugate[MDO])*Lam[ct1, ct3, ct2]*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3] - Conjugate[ZD[gt2, 4]]*ZD[gt3, 4] - Conjugate[ZD[gt2, 5]]*ZD[gt3, 5] - Conjugate[ZD[gt2, 6]]*ZD[gt3, 6]),
(-I/2)*g3*((MDO + Conjugate[MDO])*(dZgs1+dZSd1[gt2]/2+dZSd1[gt3]/2+dZPOc1/2)+2*dMGl1)*Lam[ct1, ct3, ct2]*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3] - Conjugate[ZD[gt2, 4]]*ZD[gt3, 4] - Conjugate[ZD[gt2, 5]]*ZD[gt3, 5] - Conjugate[ZD[gt2, 6]]*ZD[gt3, 6])
}},

 C[S[25, {ct1}], S[13, {gt2, ct2}], -S[13, {gt3, ct3}]] == {{
(-I/2)*g3*(MDO + Conjugate[MDO])*Lam[ct1, ct3, ct2]*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3] - Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] - Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] - Conjugate[ZU[gt2, 6]]*ZU[gt3, 6]),
(-I/2)*g3*((MDO + Conjugate[MDO])*(dZgs1+dZSu1[gt2]/2+dZSu1[gt3]/2+dZPOc1/2)+2*dMGl1)*Lam[ct1, ct3, ct2]*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3] - Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] - Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] - Conjugate[ZU[gt2, 6]]*ZU[gt3, 6])
}},

 C[S[14, {gt1, ct1}], S[26, {ct2}], -S[14, {gt3, ct3}]] == {{
(g3*(MDO - Conjugate[MDO])*Lam[ct2, ct3, ct1]*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3] - Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] - Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] - Conjugate[ZD[gt1, 6]]*ZD[gt3, 6]))/2,
(g3*((MDO - Conjugate[MDO])*(dZgs1+dZSd1[gt2]/2+dZSd1[gt3]/2+dZSOc1/2))*Lam[ct2, ct3, ct1]*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3] - Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] - Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] - Conjugate[ZD[gt1, 6]]*ZD[gt3, 6]))/2}},

 C[S[26, {ct1}], S[13, {gt2, ct2}], -S[13, {gt3, ct3}]] == {{
(g3*(MDO - Conjugate[MDO])*Lam[ct1, ct3, ct2]*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3] - Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] - Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] - Conjugate[ZU[gt2, 6]]*ZU[gt3, 6]))/2,
(g3*(MDO - Conjugate[MDO])*(dZgs1+dZSu1[gt2]/2+dZSu1[gt3]/2+dZSOc1/2)*Lam[ct1, ct3, ct2]*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3] - Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] - Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] - Conjugate[ZU[gt2, 6]]*ZU[gt3, 6]))/2}},

 C[S[25, {ct1}], S[25, {ct2}], V[5, {ct3}]] == {{g3*fSU3[ct1, ct2, ct3],(dZgs1+dZGG1/2+dZPOc1)*g3*fSU3[ct1, ct2, ct3]}},

 C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[5, {ct3}]] == {{(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct2, ct1],(-I/2)*(dZgs1+dZGG1/2+dZSd1[gt1]/2+dZSd1[gt2]/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct2, ct1]}},

 (* C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[1]] == {{(-I/6)*IndexDelta[ct1, ct2]*((CTW*g1 - 3*g2*STW)*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) - 2*CTW*g1*(Conjugate[ZD[gt1, 4]]*ZD[gt2, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt2, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt2, 6]))}}, *)
 (* C[S[14, {gt1, ct1}], -S[14, {gt2, ct2}], V[2]] == {{(I/6)*IndexDelta[ct1, ct2]*((3*CTW*g2 + g1*STW)*(Conjugate[ZD[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt2, 3]) - 2*g1*STW*(Conjugate[ZD[gt1, 4]]*ZD[gt2, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt2, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt2, 6]))}}, *)
 (* C[S[14, {gt1, ct1}], -S[13, {gt2, ct2}], -V[3]] == {{((-I)*g2*IndexDelta[ct1, ct2]*(Conjugate[ZD[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZD[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZD[gt1, 3]]*ZU[gt2, 3]))/Sqrt[2]}}, *)

 C[S[26, {ct1}], S[26, {ct2}], V[5, {ct3}]] == {g3*fSU3[ct1, ct2, ct3]*{1,dZgs1+dZGG1/2+dZSOc1}},

 C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[5, {ct3}]] == {{(-I/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct2, ct1],(-I/2)*(dZgs1+dZGG1/2+dZSu1[gt1]/2+dZSu1[gt2]/2)*g3*IndexDelta[gt1, gt2]*Lam[ct3, ct2, ct1]}},

 (* C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[1]] == {{(-I/6)*IndexDelta[ct1, ct2]*((CTW*g1 + 3*g2*STW)*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) + 4*CTW*g1*(Conjugate[ZU[gt1, 4]]*ZU[gt2, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt2, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt2, 6]))}}, *)
 (* C[S[13, {gt1, ct1}], -S[14, {gt2, ct2}], V[3]] == {{((-I)*g2*IndexDelta[ct1, ct2]*(Conjugate[ZU[gt1, 1]]*ZD[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZD[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZD[gt2, 3]))/Sqrt[2]}}, *)
 (* C[S[13, {gt1, ct1}], -S[13, {gt2, ct2}], V[2]] == {{(-I/6)*IndexDelta[ct1, ct2]*((3*CTW*g2 - g1*STW)*(Conjugate[ZU[gt1, 1]]*ZU[gt2, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt2, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt2, 3]) - 4*g1*STW*(Conjugate[ZU[gt1, 4]]*ZU[gt2, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt2, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt2, 6]))}}, *)

 C[S[25, {ct1}], S[25, {ct2}], S[26, {ct3}], S[26, {ct4}]] == {{(-I)*g3^2*(SUNF[ct1, ct3, ct2, ct4] + SUNF[ct1, ct4, ct2, ct3])
(* (fSU3[1, ct1, ct4]*fSU3[1, ct2, ct3] + fSU3[1, ct1, ct3]*fSU3[1, ct2, ct4] + fSU3[2, ct1, ct4]*fSU3[2, ct2, ct3] + fSU3[2, ct1, ct3]*fSU3[2, ct2, ct4] + fSU3[3, ct1, ct4]*fSU3[3, ct2, ct3] + fSU3[3, ct1, ct3]*fSU3[3, ct2, ct4] + fSU3[4, ct1, ct4]*fSU3[4, ct2, ct3] + fSU3[4, ct1, ct3]*fSU3[4, ct2, ct4] + fSU3[5, ct1, ct4]*fSU3[5, ct2, ct3] + fSU3[5, ct1, ct3]*fSU3[5, ct2, ct4] + fSU3[6, ct1, ct4]*fSU3[6, ct2, ct3] + fSU3[6, ct1, ct3]*fSU3[6, ct2, ct4] + fSU3[7, ct1, ct4]*fSU3[7, ct2, ct3] + fSU3[7, ct1, ct3]*fSU3[7, ct2, ct4] + fSU3[8, ct1, ct4]*fSU3[8, ct2, ct3] + fSU3[8, ct1, ct3]*fSU3[8, ct2, ct4]) *)}},

 C[S[25, {ct1}], S[14, {gt2, ct2}], S[26, {ct3}], -S[14, {gt4, ct4}]] == {{(*(I/2)*g3^2*(fSU3[1, ct1, ct3]*Lam[1, ct4, ct2] + fSU3[2, ct1, ct3]*Lam[2, ct4, ct2] + fSU3[3, ct1, ct3]*Lam[3, ct4, ct2] + fSU3[4, ct1, ct3]*Lam[4, ct4, ct2] + fSU3[5, ct1, ct3]*Lam[5, ct4, ct2] + fSU3[6, ct1, ct3]*Lam[6, ct4, ct2] + fSU3[7, ct1, ct3]*Lam[7, ct4, ct2] + fSU3[8, ct1, ct3]*Lam[8, ct4, ct2])*)
g3^2*(SUNT[ct1,ct3,ct4,ct2]-SUNT[ct3,ct1,ct4,ct2])*(-(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1]) - Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] - Conjugate[ZD[gt2, 3]]*ZD[gt4, 3] + Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6])}},

C[S[25, {ct1}], S[26, {ct2}], S[13, {gt3, ct3}], -S[13, {gt4, ct4}]] == {{(*(I/2)*g3^2*(fSU3[1, ct1, ct2]*Lam[1, ct4, ct3] + fSU3[2, ct1, ct2]*Lam[2, ct4, ct3] + fSU3[3, ct1, ct2]*Lam[3, ct4, ct3] + fSU3[4, ct1, ct2]*Lam[4, ct4, ct3] + fSU3[5, ct1, ct2]*Lam[5, ct4, ct3] + fSU3[6, ct1, ct2]*Lam[6, ct4, ct3] + fSU3[7, ct1, ct2]*Lam[7, ct4, ct3] + fSU3[8, ct1, ct2]*Lam[8, ct4, ct3])*)
g3^2*(SUNT[ct1,ct2,ct4,ct3]-SUNT[ct2,ct1,ct4,ct3])*(-(Conjugate[ZU[gt3, 1]]*ZU[gt4, 1]) - Conjugate[ZU[gt3, 2]]*ZU[gt4, 2] - Conjugate[ZU[gt3, 3]]*ZU[gt4, 3] + Conjugate[ZU[gt3, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt3, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt3, 6]]*ZU[gt4, 6])}},

 C[S[14, {gt1, ct1}], S[14, {gt2, ct2}], -S[14, {gt3, ct3}], -S[14, {gt4, ct4}]] == {{(I/72)*(-(IndexDelta[ct1, ct4]*IndexDelta[ct2, ct3]*(2*g1^2*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt1, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt4, 3]) + 18*g2^2*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt1, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt4, 3]) - 12*g3^2*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt1, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt4, 3]) + 4*g1^2*(Conjugate[ZD[gt2, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt1, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt4, 3]) + 12*g3^2*(Conjugate[ZD[gt2, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt1, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt4, 3]) + 72*(Conjugate[ZD[gt2, 1]]*(Yd[1, 1]*ZD[gt3, 4] + Yd[2, 1]*ZD[gt3, 5] + Yd[3, 1]*ZD[gt3, 6]) + Conjugate[ZD[gt2, 2]]*(Yd[1, 2]*ZD[gt3, 4] + Yd[2, 2]*ZD[gt3, 5] + Yd[3, 2]*ZD[gt3, 6]) + Conjugate[ZD[gt2, 3]]*(Yd[1, 3]*ZD[gt3, 4] + Yd[2, 3]*ZD[gt3, 5] + Yd[3, 3]*ZD[gt3, 6]))*((Conjugate[Yd[1, 1]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 1]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 1]]*Conjugate[ZD[gt1, 6]])*ZD[gt4, 1] + (Conjugate[Yd[1, 2]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 2]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 2]]*Conjugate[ZD[gt1, 6]])*ZD[gt4, 2] + (Conjugate[Yd[1, 3]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 3]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 3]]*Conjugate[ZD[gt1, 6]])*ZD[gt4, 3]) + 18*g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) + 18*g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3] - Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] - Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] - Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) - 18*g3^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) + 4*g1^2*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt1, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt4, 6]) + 12*g3^2*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt1, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt4, 6]) + 8*g1^2*(Conjugate[ZD[gt2, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt1, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt4, 6]) - 12*g3^2*(Conjugate[ZD[gt2, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt1, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt4, 6]) - 18*g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6]) - 18*g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3] - Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] - Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] - Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6]) + 18*g3^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6]) + 72*((Conjugate[Yd[1, 1]]*Conjugate[ZD[gt2, 4]] + Conjugate[Yd[2, 1]]*Conjugate[ZD[gt2, 5]] + Conjugate[Yd[3, 1]]*Conjugate[ZD[gt2, 6]])*ZD[gt3, 1] + (Conjugate[Yd[1, 2]]*Conjugate[ZD[gt2, 4]] + Conjugate[Yd[2, 2]]*Conjugate[ZD[gt2, 5]] + Conjugate[Yd[3, 2]]*Conjugate[ZD[gt2, 6]])*ZD[gt3, 2] + (Conjugate[Yd[1, 3]]*Conjugate[ZD[gt2, 4]] + Conjugate[Yd[2, 3]]*Conjugate[ZD[gt2, 5]] + Conjugate[Yd[3, 3]]*Conjugate[ZD[gt2, 6]])*ZD[gt3, 3])*(Conjugate[ZD[gt1, 1]]*(Yd[1, 1]*ZD[gt4, 4] + Yd[2, 1]*ZD[gt4, 5] + Yd[3, 1]*ZD[gt4, 6]) + Conjugate[ZD[gt1, 2]]*(Yd[1, 2]*ZD[gt4, 4] + Yd[2, 2]*ZD[gt4, 5] + Yd[3, 2]*ZD[gt4, 6]) + Conjugate[ZD[gt1, 3]]*(Yd[1, 3]*ZD[gt4, 4] + Yd[2, 3]*ZD[gt4, 5] + Yd[3, 3]*ZD[gt4, 6])))) - IndexDelta[ct1, ct3]*IndexDelta[ct2, ct4]*(36*g3^2*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt1, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt4, 3]) - 36*g3^2*(Conjugate[ZD[gt2, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt1, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt4, 3]) + g1^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) + 9*g2^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) - 6*g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) + 2*g1^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) + 6*g3^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) + ((g1^2 + 9*g2^2 - 6*g3^2)*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3]) + 2*(g1^2 + 3*g3^2)*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6]))*(Conjugate[ZD[gt2, 1]]*ZD[gt4, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt4, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt4, 3]) + 72*(Conjugate[ZD[gt1, 1]]*(Yd[1, 1]*ZD[gt3, 4] + Yd[2, 1]*ZD[gt3, 5] + Yd[3, 1]*ZD[gt3, 6]) + Conjugate[ZD[gt1, 2]]*(Yd[1, 2]*ZD[gt3, 4] + Yd[2, 2]*ZD[gt3, 5] + Yd[3, 2]*ZD[gt3, 6]) + Conjugate[ZD[gt1, 3]]*(Yd[1, 3]*ZD[gt3, 4] + Yd[2, 3]*ZD[gt3, 5] + Yd[3, 3]*ZD[gt3, 6]))*((Conjugate[Yd[1, 1]]*Conjugate[ZD[gt2, 4]] + Conjugate[Yd[2, 1]]*Conjugate[ZD[gt2, 5]] + Conjugate[Yd[3, 1]]*Conjugate[ZD[gt2, 6]])*ZD[gt4, 1] + (Conjugate[Yd[1, 2]]*Conjugate[ZD[gt2, 4]] + Conjugate[Yd[2, 2]]*Conjugate[ZD[gt2, 5]] + Conjugate[Yd[3, 2]]*Conjugate[ZD[gt2, 6]])*ZD[gt4, 2] + (Conjugate[Yd[1, 3]]*Conjugate[ZD[gt2, 4]] + Conjugate[Yd[2, 3]]*Conjugate[ZD[gt2, 5]] + Conjugate[Yd[3, 3]]*Conjugate[ZD[gt2, 6]])*ZD[gt4, 3]) - 36*g3^2*(Conjugate[ZD[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt2, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt1, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt4, 6]) + 36*g3^2*(Conjugate[ZD[gt2, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt1, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt4, 6]) + 2*g1^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6]) + 6*g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6]) + 4*g1^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6]) - 6*g3^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6]) + 2*((g1^2 + 3*g3^2)*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3]) + (2*g1^2 - 3*g3^2)*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6]))*(Conjugate[ZD[gt2, 4]]*ZD[gt4, 4] + Conjugate[ZD[gt2, 5]]*ZD[gt4, 5] + Conjugate[ZD[gt2, 6]]*ZD[gt4, 6]) + 72*((Conjugate[Yd[1, 1]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 1]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 1]]*Conjugate[ZD[gt1, 6]])*ZD[gt3, 1] + (Conjugate[Yd[1, 2]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 2]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 2]]*Conjugate[ZD[gt1, 6]])*ZD[gt3, 2] + (Conjugate[Yd[1, 3]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 3]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 3]]*Conjugate[ZD[gt1, 6]])*ZD[gt3, 3])*(Conjugate[ZD[gt2, 1]]*(Yd[1, 1]*ZD[gt4, 4] + Yd[2, 1]*ZD[gt4, 5] + Yd[3, 1]*ZD[gt4, 6]) + Conjugate[ZD[gt2, 2]]*(Yd[1, 2]*ZD[gt4, 4] + Yd[2, 2]*ZD[gt4, 5] + Yd[3, 2]*ZD[gt4, 6]) + Conjugate[ZD[gt2, 3]]*(Yd[1, 3]*ZD[gt4, 4] + Yd[2, 3]*ZD[gt4, 5] + Yd[3, 3]*ZD[gt4, 6]))))}},
 C[S[14, {gt1, ct1}], S[13, {gt2, ct2}], -S[14, {gt3, ct3}], -S[13, {gt4, ct4}]] == {{(I/72)*(IndexDelta[ct1, ct3]*IndexDelta[ct2, ct4]*(-(g1^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3])) + 9*g2^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) + 6*g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - 2*g1^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - 6*g3^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - ((g1^2 - 9*g2^2 - 6*g3^2)*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3]) + 2*(g1^2 + 3*g3^2)*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6]))*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) + 4*g1^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) - 6*g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + 8*g1^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + 6*g3^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + ((4*g1^2 - 6*g3^2)*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3]) + 2*(4*g1^2 + 3*g3^2)*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6]))*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6])) - 18*IndexDelta[ct1, ct4]*IndexDelta[ct2, ct3]*(2*g2^2*(Conjugate[ZU[gt2, 1]]*ZD[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZD[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZD[gt3, 3])*(Conjugate[ZD[gt1, 1]]*ZU[gt4, 1] + Conjugate[ZD[gt1, 2]]*ZU[gt4, 2] + Conjugate[ZD[gt1, 3]]*ZU[gt4, 3]) + 4*(Conjugate[ZU[gt2, 1]]*(Yd[1, 1]*ZD[gt3, 4] + Yd[2, 1]*ZD[gt3, 5] + Yd[3, 1]*ZD[gt3, 6]) + Conjugate[ZU[gt2, 2]]*(Yd[1, 2]*ZD[gt3, 4] + Yd[2, 2]*ZD[gt3, 5] + Yd[3, 2]*ZD[gt3, 6]) + Conjugate[ZU[gt2, 3]]*(Yd[1, 3]*ZD[gt3, 4] + Yd[2, 3]*ZD[gt3, 5] + Yd[3, 3]*ZD[gt3, 6]))*((Conjugate[Yd[1, 1]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 1]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 1]]*Conjugate[ZD[gt1, 6]])*ZU[gt4, 1] + (Conjugate[Yd[1, 2]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 2]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 2]]*Conjugate[ZD[gt1, 6]])*ZU[gt4, 2] + (Conjugate[Yd[1, 3]]*Conjugate[ZD[gt1, 4]] + Conjugate[Yd[2, 3]]*Conjugate[ZD[gt1, 5]] + Conjugate[Yd[3, 3]]*Conjugate[ZD[gt1, 6]])*ZU[gt4, 3]) + g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) + g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3] - Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] - Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] - Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - g3^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - g3^2*(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1] + Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] + Conjugate[ZD[gt1, 3]]*ZD[gt3, 3])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + g3^2*(Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + g3^2*(-(Conjugate[ZD[gt1, 1]]*ZD[gt3, 1]) - Conjugate[ZD[gt1, 2]]*ZD[gt3, 2] - Conjugate[ZD[gt1, 3]]*ZD[gt3, 3] + Conjugate[ZD[gt1, 4]]*ZD[gt3, 4] + Conjugate[ZD[gt1, 5]]*ZD[gt3, 5] + Conjugate[ZD[gt1, 6]]*ZD[gt3, 6])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + 4*((Conjugate[Yu[1, 1]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 1]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 1]]*Conjugate[ZU[gt2, 6]])*ZD[gt3, 1] + (Conjugate[Yu[1, 2]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 2]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 2]]*Conjugate[ZU[gt2, 6]])*ZD[gt3, 2] + (Conjugate[Yu[1, 3]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 3]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 3]]*Conjugate[ZU[gt2, 6]])*ZD[gt3, 3])*(Conjugate[ZD[gt1, 1]]*(Yu[1, 1]*ZU[gt4, 4] + Yu[2, 1]*ZU[gt4, 5] + Yu[3, 1]*ZU[gt4, 6]) + Conjugate[ZD[gt1, 2]]*(Yu[1, 2]*ZU[gt4, 4] + Yu[2, 2]*ZU[gt4, 5] + Yu[3, 2]*ZU[gt4, 6]) + Conjugate[ZD[gt1, 3]]*(Yu[1, 3]*ZU[gt4, 4] + Yu[2, 3]*ZU[gt4, 5] + Yu[3, 3]*ZU[gt4, 6]))))}},
 C[S[13, {gt1, ct1}], S[13, {gt2, ct2}], -S[13, {gt3, ct3}], -S[13, {gt4, ct4}]] == {{(I/72)*(-(IndexDelta[ct1, ct4]*IndexDelta[ct2, ct3]*(2*g1^2*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt1, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt4, 3]) + 18*g2^2*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt1, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt4, 3]) - 12*g3^2*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt1, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt4, 3]) - 8*g1^2*(Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt1, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt4, 3]) + 12*g3^2*(Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt1, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt4, 3]) + 72*(Conjugate[ZU[gt2, 1]]*(Yu[1, 1]*ZU[gt3, 4] + Yu[2, 1]*ZU[gt3, 5] + Yu[3, 1]*ZU[gt3, 6]) + Conjugate[ZU[gt2, 2]]*(Yu[1, 2]*ZU[gt3, 4] + Yu[2, 2]*ZU[gt3, 5] + Yu[3, 2]*ZU[gt3, 6]) + Conjugate[ZU[gt2, 3]]*(Yu[1, 3]*ZU[gt3, 4] + Yu[2, 3]*ZU[gt3, 5] + Yu[3, 3]*ZU[gt3, 6]))*((Conjugate[Yu[1, 1]]*Conjugate[ZU[gt1, 4]] + Conjugate[Yu[2, 1]]*Conjugate[ZU[gt1, 5]] + Conjugate[Yu[3, 1]]*Conjugate[ZU[gt1, 6]])*ZU[gt4, 1] + (Conjugate[Yu[1, 2]]*Conjugate[ZU[gt1, 4]] + Conjugate[Yu[2, 2]]*Conjugate[ZU[gt1, 5]] + Conjugate[Yu[3, 2]]*Conjugate[ZU[gt1, 6]])*ZU[gt4, 2] + (Conjugate[Yu[1, 3]]*Conjugate[ZU[gt1, 4]] + Conjugate[Yu[2, 3]]*Conjugate[ZU[gt1, 5]] + Conjugate[Yu[3, 3]]*Conjugate[ZU[gt1, 6]])*ZU[gt4, 3]) + 18*g3^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) + 18*g3^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3] - Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] - Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] - Conjugate[ZU[gt1, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - 18*g3^2*(Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - 8*g1^2*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt1, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt4, 6]) + 12*g3^2*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt1, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt4, 6]) + 32*g1^2*(Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt1, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt4, 6]) - 12*g3^2*(Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt1, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt4, 6]) - 18*g3^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) - 18*g3^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3] - Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] - Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] - Conjugate[ZU[gt1, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + 18*g3^2*(Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + 72*((Conjugate[Yu[1, 1]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 1]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 1]]*Conjugate[ZU[gt2, 6]])*ZU[gt3, 1] + (Conjugate[Yu[1, 2]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 2]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 2]]*Conjugate[ZU[gt2, 6]])*ZU[gt3, 2] + (Conjugate[Yu[1, 3]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 3]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 3]]*Conjugate[ZU[gt2, 6]])*ZU[gt3, 3])*(Conjugate[ZU[gt1, 1]]*(Yu[1, 1]*ZU[gt4, 4] + Yu[2, 1]*ZU[gt4, 5] + Yu[3, 1]*ZU[gt4, 6]) + Conjugate[ZU[gt1, 2]]*(Yu[1, 2]*ZU[gt4, 4] + Yu[2, 2]*ZU[gt4, 5] + Yu[3, 2]*ZU[gt4, 6]) + Conjugate[ZU[gt1, 3]]*(Yu[1, 3]*ZU[gt4, 4] + Yu[2, 3]*ZU[gt4, 5] + Yu[3, 3]*ZU[gt4, 6])))) - IndexDelta[ct1, ct3]*IndexDelta[ct2, ct4]*(36*g3^2*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt1, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt4, 3]) - 36*g3^2*(Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt1, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt4, 3]) + g1^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) + 9*g2^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - 6*g3^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) - 4*g1^2*(Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) + 6*g3^2*(Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) + ((g1^2 + 9*g2^2 - 6*g3^2)*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3]) + 2*(-2*g1^2 + 3*g3^2)*(Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt3, 6]))*(Conjugate[ZU[gt2, 1]]*ZU[gt4, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt4, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt4, 3]) + 72*(Conjugate[ZU[gt1, 1]]*(Yu[1, 1]*ZU[gt3, 4] + Yu[2, 1]*ZU[gt3, 5] + Yu[3, 1]*ZU[gt3, 6]) + Conjugate[ZU[gt1, 2]]*(Yu[1, 2]*ZU[gt3, 4] + Yu[2, 2]*ZU[gt3, 5] + Yu[3, 2]*ZU[gt3, 6]) + Conjugate[ZU[gt1, 3]]*(Yu[1, 3]*ZU[gt3, 4] + Yu[2, 3]*ZU[gt3, 5] + Yu[3, 3]*ZU[gt3, 6]))*((Conjugate[Yu[1, 1]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 1]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 1]]*Conjugate[ZU[gt2, 6]])*ZU[gt4, 1] + (Conjugate[Yu[1, 2]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 2]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 2]]*Conjugate[ZU[gt2, 6]])*ZU[gt4, 2] + (Conjugate[Yu[1, 3]]*Conjugate[ZU[gt2, 4]] + Conjugate[Yu[2, 3]]*Conjugate[ZU[gt2, 5]] + Conjugate[Yu[3, 3]]*Conjugate[ZU[gt2, 6]])*ZU[gt4, 3]) - 36*g3^2*(Conjugate[ZU[gt2, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt2, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt2, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt1, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt4, 6]) + 36*g3^2*(Conjugate[ZU[gt2, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt1, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt4, 6]) - 4*g1^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + 6*g3^2*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + 16*g1^2*(Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) - 6*g3^2*(Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt3, 6])*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + ((-4*g1^2 + 6*g3^2)*(Conjugate[ZU[gt1, 1]]*ZU[gt3, 1] + Conjugate[ZU[gt1, 2]]*ZU[gt3, 2] + Conjugate[ZU[gt1, 3]]*ZU[gt3, 3]) + 2*(8*g1^2 - 3*g3^2)*(Conjugate[ZU[gt1, 4]]*ZU[gt3, 4] + Conjugate[ZU[gt1, 5]]*ZU[gt3, 5] + Conjugate[ZU[gt1, 6]]*ZU[gt3, 6]))*(Conjugate[ZU[gt2, 4]]*ZU[gt4, 4] + Conjugate[ZU[gt2, 5]]*ZU[gt4, 5] + Conjugate[ZU[gt2, 6]]*ZU[gt4, 6]) + 72*((Conjugate[Yu[1, 1]]*Conjugate[ZU[gt1, 4]] + Conjugate[Yu[2, 1]]*Conjugate[ZU[gt1, 5]] + Conjugate[Yu[3, 1]]*Conjugate[ZU[gt1, 6]])*ZU[gt3, 1] + (Conjugate[Yu[1, 2]]*Conjugate[ZU[gt1, 4]] + Conjugate[Yu[2, 2]]*Conjugate[ZU[gt1, 5]] + Conjugate[Yu[3, 2]]*Conjugate[ZU[gt1, 6]])*ZU[gt3, 2] + (Conjugate[Yu[1, 3]]*Conjugate[ZU[gt1, 4]] + Conjugate[Yu[2, 3]]*Conjugate[ZU[gt1, 5]] + Conjugate[Yu[3, 3]]*Conjugate[ZU[gt1, 6]])*ZU[gt3, 3])*(Conjugate[ZU[gt2, 1]]*(Yu[1, 1]*ZU[gt4, 4] + Yu[2, 1]*ZU[gt4, 5] + Yu[3, 1]*ZU[gt4, 6]) + Conjugate[ZU[gt2, 2]]*(Yu[1, 2]*ZU[gt4, 4] + Yu[2, 2]*ZU[gt4, 5] + Yu[3, 2]*ZU[gt4, 6]) + Conjugate[ZU[gt2, 3]]*(Yu[1, 3]*ZU[gt4, 4] + Yu[2, 3]*ZU[gt4, 5] + Yu[3, 3]*ZU[gt4, 6]))))}},


 (* C[-V[3], V[1], V[3]] == {{I*g2*STW}}, *)
 (* C[-V[3], V[3], V[2]] == {{(-I)*CTW*g2}}, *)

 C[-U[5, {ct1}], U[5, {ct2}], V[5, {ct3}]] == {{g3*fSU3[ct1, ct2, ct3]}, {0}},

 (* C[-U[3], U[1], V[3]] == {{I*g2*STW}, {0}}, *)
 (* C[-U[4], U[1], -V[3]] == {{(-I)*g2*STW}, {0}}, *)
 (* C[-U[3], U[3], V[1]] == {{(-I)*g2*STW}, {0}}, *)
 (* C[-U[3], U[3], V[2]] == {{(-I)*CTW*g2}, {0}}, *)
 (* C[-U[1], U[3], -V[3]] == {{I*g2*STW}, {0}}, *)
 (* C[-U[2], U[3], -V[3]] == {{I*CTW*g2}, {0}}, *)
 (* C[-U[4], U[4], V[1]] == {{I*g2*STW}, {0}}, *)
 (* C[-U[1], U[4], V[3]] == {{(-I)*g2*STW}, {0}}, *)
 (* C[-U[2], U[4], V[3]] == {{(-I)*CTW*g2}, {0}}, *)
 (* C[-U[4], U[4], V[2]] == {{I*CTW*g2}, {0}}, *)
 (* C[-U[3], U[2], V[3]] == {{I*CTW*g2}, {0}}, *)
 (* C[-U[4], U[2], -V[3]] == {{(-I)*CTW*g2}, {0}}, *)
 (* C[-F[3, {j1, o1}], F[3, {j2, o2}]] ==  *)
 (*  {{0, (-I/2)*(dZbarfL1[3, j2, j1] + dZfL1[3, j1, j2])*IndexDelta[o1, o2]},  *)
 (*   {0, (I/2)*(dZbarfR1[3, j2, j1] + dZfR1[3, j1, j2])*IndexDelta[o1, o2]},  *)
 (*   {0, (-I/2)*IndexDelta[o1, o2]*(2*dMf1[3, j1]*IndexDelta[j1, j2] +  *)
 (*      dZfL1[3, j1, j2]*Mass[F[3, {j1}]] + dZbarfR1[3, j1, j2]* *)
 (*       Mass[F[3, {j2}]])}, {0, (-I/2)*IndexDelta[o1, o2]* *)
 (*     (2*dMf1[3, j1]*IndexDelta[j1, j2] + dZfR1[3, j1, j2]*Mass[F[3, {j1}]] +  *)
 (*      dZbarfL1[3, j1, j2]*Mass[F[3, {j2}]])}},  *)
 (* C[-F[4, {j1, o1}], F[4, {j2, o2}]] ==  *)
 (*  {{0, (-I/2)*(dZbarfL1[4, j2, j1] + dZfL1[4, j1, j2])*IndexDelta[o1, o2]},  *)
 (*   {0, (I/2)*(dZbarfR1[4, j2, j1] + dZfR1[4, j1, j2])*IndexDelta[o1, o2]},  *)
 (*   {0, (-I/2)*IndexDelta[o1, o2]*(2*dMf1[4, j1]*IndexDelta[j1, j2] +  *)
 (*      dZfL1[4, j1, j2]*Mass[F[4, {j1}]] + dZbarfR1[4, j1, j2]* *)
 (*       Mass[F[4, {j2}]])}, {0, (-I/2)*IndexDelta[o1, o2]* *)
 (*     (2*dMf1[4, j1]*IndexDelta[j1, j2] + dZfR1[4, j1, j2]*Mass[F[4, {j1}]] +  *)
 (*      dZbarfL1[4, j1, j2]*Mass[F[4, {j2}]])}}, *)
C[-F[3, {gt1, ct1}], F[3, {gt2, ct2}]] == {
{0, -(1/2) I (dZFu1[gt1] + dZFu1[gt2]) IndexDelta[ct1,ct2] IndexDelta[gt1,gt2]}, 
{0, 1/2 I (dZFu1[gt1] + dZFu1[gt2]) IndexDelta[ct1, ct2] IndexDelta[gt1,gt2]}, {0, 0}, {0, 0}}, 

C[-F[4, {gt1, ct1}], F[4, {gt2, ct2}]] == {
{0, -(1/2) I (dZFd1[gt1] + dZFd1[gt2]) IndexDelta[ct1,ct2] IndexDelta[gt1,gt2]}, 
{0, 1/2 I (dZFd1[gt1] + dZFd1[gt2]) IndexDelta[ct1, ct2] IndexDelta[gt1,gt2]}, {0, 0}, {0, 0}},

C[-S[13, {gt1,ct1}], S[13, {gt2, ct2}]] == 
  {{0, (-I/2)*(dZSu1[gt1] + dZSu1[gt2])*IndexDelta[gt1, gt2]*IndexDelta[ct1, ct2]}, 
   {0, (-I/2)*IndexDelta[gt1, gt2]*IndexDelta[ct1, ct2]*
     (2*dMSu1[gt1] + dZSu1[gt1]*
       TheMass[S[13, {gt1, ct1}]]^2 + dZSu1[gt2]*
       TheMass[S[13, {gt2, ct2}]]^2)}},
 
 C[-S[14, {gt1,ct1}], S[14, {gt2, ct2}]] == 
  {{0, (-I/2)*(dZSd1[gt1] + dZSd1[gt2])*IndexDelta[gt1, gt2]*IndexDelta[ct1, ct2]}, 
   {0, (-I/2)*IndexDelta[gt1, gt2]*IndexDelta[ct1, ct2]*
     (2*dMSd1[gt1] + dZSd1[gt1]*
       TheMass[S[14, {gt1, ct1}]]^2 + dZSd1[gt2]*
       TheMass[S[14, {gt2, ct2}]]^2)}},
  (* {{0, (-I/2)*(dZbarSf1[s2, s1, 4, j2] + dZSf1[s1, s2, 4, j1])* *)
  (*    IndexDelta[j1, j2]*IndexDelta[o1, o2]},  *)
  (*  {0, (-I/2)*IndexDelta[j1, j2]*IndexDelta[o1, o2]* *)
  (*    (2*dMSfsq1[s1, s2, 4, j1] + dZSf1[s1, s2, 4, j1]* *)
  (*      TheMass[S[14, {s1, j1}]]^2 + dZbarSf1[s2, s1, 4, j2]* *)
  (*      TheMass[S[14, {s2, j2}]]^2)}},  *)

 C[S[25, {ct1}], S[25,{ct2}]] == 
  {{0, (-I)*(dZPOc1)*
     IndexDelta[ct1, ct2]}, 
   {0, (-I/2)*IndexDelta[ct1, ct2]*
     (2*dMO21+8*dMGl1 + dZPOc1*
       TheMass[S[25, {ct1}]]^2 + dZPOc1*
       TheMass[S[25, {ct2}]]^2)}}, 

 C[S[26, {ct1}], S[26,{ct2}]] == 
  {{0, (-I)*(dZSOc1)*
     IndexDelta[ct1, ct2]}, 
   {0, (-I/2)*IndexDelta[ct1, ct2]*
     (2*dMO21 + dZSOc1*
       TheMass[S[26, {ct1}]]^2 + dZSOc1*
       TheMass[S[26, {ct2}]]^2)}}, 

C[V[5, {ct1}], V[5, {ct2}]] == 
  {{0, I*dZGG1*IndexDelta[ct1, ct2]}, {0, 0}, 
   {0, (-I)*dZGG1*IndexDelta[ct1, ct2]}}, 

C[-F[15, {ct1}], F[15, {ct2}]] == 
  {{0, (-I/2)*(dZGlL1 + dZGlL1)*IndexDelta[ct1, ct2]}, 
   {0, (I/2)*(dZGlR1 + dZGlR1)*IndexDelta[ct1, ct2]}, 
   {0, (-I/2)*IndexDelta[ct1, ct2]*(2*dMGl1 + (dZGlR1 + dZGlL1)*
       TheMass[F[15, {ct1}]])}, {0, (-I/2)*IndexDelta[ct1, ct2]*
     (2*Conjugate[dMGl1] + (dZGlL1 + dZGlR1)*TheMass[F[15, {ct1}]])}}
 

}

 
Conjugate[g1] ^= g1; 
Conjugate[g2] ^= g2; 
Conjugate[g3] ^= g3; 
Conjugate[mHd2] ^= mHd2; 
Conjugate[mHu2] ^= mHu2; 
Conjugate[mRd2] ^= mRd2; 
Conjugate[mRu2] ^= mRu2; 
Conjugate[vd] ^= vd; 
Conjugate[vu] ^= vu; 
Conjugate[vT] ^= vT; 
Conjugate[vS] ^= vS; 
Conjugate[ZZ[a___]] ^= ZZ[a]; 
Conjugate[aEWinv] ^= aEWinv; 
Conjugate[v] ^= v; 
Conjugate[Gf] ^= Gf; 

Conjugate[MDO] ^= MDO;
Conjugate[MassGlu] ^= MassGlu;
Conjugate[MasssigmaO] ^= MasssigmaO;
Conjugate[MassphiO] ^= MassphiO;
Conjugate[MassFu[i_]] ^= MassFu[i];
Conjugate[MassFd[i_]] ^= MassFd[i];
Conjugate[MassSu[i_]] ^= MassSu[i];
Conjugate[MassSd[i_]] ^= MassSd[i];


(* ------------------------ Renormalization constants ---------------------- *)
(*taken from MSSMCT.mod*)


Clear[RenConst, xHC, xZf1, xZfL1, xZfR1, xZbarfL1, xZbarfR1, xZSf1]

ReDiag = Identity

ReOffDiag = Identity


xHC[Conjugate[mat_[j_, i_]]] := mat[i, j]

xHC[mat_[i_, j_]] := Conjugate[mat[j, i]]

xHC[x_] := Conjugate[x]


xZf1[f_, X_, dm_, s_] :=
Block[ {m = TheMass[f],
        se = ReDiag[SelfEnergy[f]],
        dse = ReDiag[DSelfEnergy[f]]},
  -X[se] - m^2 (LVectorCoeff[dse] + RVectorCoeff[dse]) -
             m (LScalarCoeff[dse] + RScalarCoeff[dse]) +
  s ((LScalarCoeff[se] - dm) - (RScalarCoeff[se] - xHC[dm]))/(2 m)
]

xZfL1[f_, dm_] := xZf1[f, LVectorCoeff, dm, +1]

xZfR1[f_, dm_] := xZf1[f, RVectorCoeff, dm, -1]

xZbarfL1[f_, dm_] := xZf1[f, LVectorCoeff, dm, -1]

xZbarfR1[f_, dm_] := xZf1[f, RVectorCoeff, dm, +1]


(*xZf1[proc_, m1_, m2_, X_, Y_, SX_, SY_, dm_] :=
Block[ {se = ReOffDiag[SelfEnergy[proc, m2]]},
  2/(m1^2 - m2^2) (m2 (m2 X[se] + m1 Y[se]) +
      m2 (SX[se] - xHC[dm]) + m1 (SY[se] - dm))
]

xZfL1[f1_, f2_, dm_] := xZf1[f2 -> f1, TheMass[f1], TheMass[f2],
  LVectorCoeff, RVectorCoeff, RScalarCoeff, LScalarCoeff, dm]

xZfR1[f1_, f2_, dm_] := xZf1[f2 -> f1, TheMass[f1], TheMass[f2],
  RVectorCoeff, LVectorCoeff, LScalarCoeff, RScalarCoeff, xHC[dm]]

xZbarfL1[f1_, f2_, dm_] := xZf1[f2 -> f1, TheMass[f2], TheMass[f1],
  LVectorCoeff, RVectorCoeff, LScalarCoeff, RScalarCoeff, xHC[dm]]

xZbarfR1[f1_, f2_, dm_] := xZf1[f2 -> f1, TheMass[f2], TheMass[f1],
  RVectorCoeff, LVectorCoeff, RScalarCoeff, LScalarCoeff, dm]*)



(* ---------------------- SM renormalization constants --------------------- *)

(*_UVMf1 = Identity

UVMf1[4, 3] = UVDivergentPart

RenConst[dMf1[t_, j1_]] := UVMf1[t, j1][MassRC[F[t, {j1}]]]


RenConst[dZfL1[t_, j1_, j1_]] := xZfL1[F[t, {j1}], 0]

RenConst[dZfR1[t_, j1_, j1_]] := xZfR1[F[t, {j1}], 0]

RenConst[dZbarfL1[t_, j1_, j1_]] := xZbarfL1[F[t, {j1}], 0]

RenConst[dZbarfR1[t_, j1_, j1_]] := xZbarfR1[F[t, {j1}], 0]


RenConst[dZfL1[t_, j1_, j2_]] := xZfL1[F[t, {j1}], F[t, {j2}], 0]

RenConst[dZfR1[t_, j1_, j2_]] := xZfR1[F[t, {j1}], F[t, {j2}], 0]

RenConst[dZbarfL1[t_, j1_, j2_]] := xZbarfL1[F[t, {j1}], F[t, {j2}], 0]*)


(* --------------------------------- SQCD ---------------------------------- *)
RenConst[dZFu1[gen_]] :=FieldRC[F[3,{gen}]][[1]]
RenConst[dZFd1[gen_]] :=FieldRC[F[4,{gen}]][[1]]  

RenConst[dZSu1[gen_]] := FieldRC[S[13,{gen}]]	
RenConst[dMSu1[gen_]] := MassRC[S[13,{gen}]]
RenConst[dZSd1[gen_]] := FieldRC[S[14,{gen}]]	
RenConst[dMSd1[gen_]] := MassRC[S[14,{gen}]]

RenConst[dZGlL1] := FieldRC[F[15]][[1]]  (*xZfL1[F[15], 0]*)
RenConst[dZGlR1] := FieldRC[F[15]][[2]]  (*xZfR1[F[15], 0]*)
RenConst[dZbarGlL1] := dZGlR1
RenConst[dZbarGlR1] := dZGlL1
RenConst[dMGl1] := MassRC[F[15]]

RenConst[dZgs1] := 0 - FiniteGs * g3*g3/(16*Pi*Pi)*(2*Log[(MassGlu/mu)^2]
                     +Log[(MassSq/mu)^2]+1/3*Log[(MassFu[3]/mu)^2]
                     +1/4*Log[(MasssigmaO/mu)^2]
                     +1/4*Log[(MassphiO/mu)^2])
dZgs1restore := FiniteGs * g3*g3/(12*Pi*Pi)

RenConst[dZGG1] := FieldRC[V[5]]
(*-------------------------------------------------------------------*)

(*Set mixing matrices to identities*) 
ZUR := KroneckerDelta 
ZUL := KroneckerDelta 
ZU := KroneckerDelta
ZDR := KroneckerDelta
ZDL := KroneckerDelta
ZD := KroneckerDelta

MassGlu[_] := MassGlu
MassphiO[_] := MassphiO
MasssigmaO[_] := MasssigmaO
