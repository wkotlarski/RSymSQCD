<< "FeynArts.m";
<< "FormCalc.m";

(* default is 1 for one-loop which is fine for UV divergences but not for IR
ones *)
$Dminus4MaxPower = 2;

model = $ScriptCommandLine[[2]];
processNr = ToExpression[$ScriptCommandLine[[3]]];

Switch[model,
   (* exclude W, Z, gamma *)
   "MSSM"|"MRSSM", SetOptions[InsertFields, Model -> model <> "EWSB_myRenConst", InsertionLevel -> {Particles}, ExcludeParticles -> {V[n_ /; n != 5]}],
   (* exclude W, Z *)
   "SMQEDQCD", SetOptions[InsertFields, Model -> "SMQCD", InsertionLevel -> {Particles}, ExcludeParticles -> {V[2], V[3]}],
   "SMQCD", SetOptions[InsertFields, Model -> "SMQCD", InsertionLevel -> {Particles}, ExcludeParticles -> {V[n_ /; n != 5]}]
];

Switch[processNr,
   1,
      process = {F[3, {1}], -F[3, {1}]} -> {S[26], S[26], V[5]};
      procName = "qqbar_OO",
   2,
      process = {V[5], V[5]} -> {S[26], S[26], V[5]};
      procName = "gg_OO",
   3,
      process = {F[3, {1}], F[3, {1}]} -> {S[13, {1}], S[13, {4}], V[5]};
      procName = "uu_suLsuRg";,
   4,
      process = {F[3, {1}], F[4, {1}]} -> {S[13, {1}], S[14, {4}], V[5]};,
   5,
      process = {F[3, {1}],-F[3, {1}]} -> {S[13, {1}],-S[13, {1}],V[5]};
      procName = "uubar_suLsuLdaggerg";,
   6,
      process = {F[4, {1}],-F[4, {1}]} -> {S[13, {1}],-S[13, {1}],V[5]};
      procName = "ddbar_suLsuLdaggerg";,
   7,
      process = {V[5], V[5]} -> {S[13, {1}], -S[13, {1}], V[5]};
      procName = "gg_suLsuLdaggerg";,
   8,
      process = {V[5], V[5]} -> {F[3, {3}], -F[3, {3}], V[5]};
      procName = "gg_ttbar";,
   9,
      process = {V[5], -F[3, {1}]} -> {S[13, {1}], -S[13, {1}], -F[3, {1}]};
      procName = "gubar_suLsuLdaggerubar";,
   10,
      process = {V[5], F[3, {1}]} -> {S[13, {1}], -S[13, {1}], F[3, {1}]};
      procName = "gu_suLsuLdaggeru";,
   11,
      process = {V[5], -F[4, {1}]} -> {S[13, {1}], -S[13, {1}], -F[4, {1}]};
      procName = "gdbar_suLsuLdaggerdbar";,
   13,
      process = {V[5], F[3, {1}]} -> {S[13, {1}], S[13, {4}], -F[3, {1}]};
      procName = "gu_suLsuRubar";,
   14,
      process = {V[5], F[4, {1}]} -> {S[13, {1}], -S[13, {1}], F[4, {1}]};
      procName = "gd_suLsuLdaggerd";,
   _, Print["Unknown process ", processNr]; Quit[1]
];

top = CreateTopologies[0, 2 -> 3];
ins = InsertFields[top, process];
(* diagram removal *)
(*ins = DiagramSelect[ins, (Not[SChannelQ[F[15]][##]] && Not[SChannelQ[-F[15]][##]])&];*)

Neglect[MassFu[1]]=0;
Neglect[MassFu[2]]=0;
Neglect[MassFd[_]]=0;

mOfLeg[i_Integer] := TheMass[Flatten[List@@process][[i]]];
Square[mOfLeg[3]];
Square[mOfLeg[4]];

(* For gg->suLsuL*g with Transverse->False FormCalc crashes because the
   expression is larger than a hardcoded maximal size. One needs to edit file
   ReadForm.tm in FormCalc directory and change TERMBUF to lets say 6000000 *)
SetOptions[CalcFeynAmp, FermionChains->Chiral];
SetOptions[PolarizationSum, GaugeTerms->Off];

amp = CalcFeynAmp[CreateFeynAmp[ins], Invariants->False, Normalized->False, Transverse->False];
result = PolarizationSum[amp] //. Abbr[] //. Subexpr[];

If[Count[Flatten[List@@process], V[5]] >= 2,
   If[processNr == 7,
      f[proc_] := Module[{ins, amp},
         ins = InsertFields[top, proc];
         amp = CalcFeynAmp[CreateFeynAmp[ins], Invariants->False, Normalized->False, Transverse->False];
         PolarizationSum[amp] //. Abbr[] //. Subexpr[]
      ];

      ghostProcesses = {
         { U[5], -U[5]} -> {S[13, {1}], -S[13, {1}],  V[5]},
         {-U[5],  U[5]} -> {S[13, {1}], -S[13, {1}],  V[5]},
         { U[5],  V[5]} -> {S[13, {1}], -S[13, {1}],  U[5]},
         {-U[5],  V[5]} -> {S[13, {1}], -S[13, {1}], -U[5]},
         { V[5],  U[5]} -> {S[13, {1}], -S[13, {1}],  U[5]},
         { V[5], -U[5]} -> {S[13, {1}], -S[13, {1}], -U[5]}
      };
      result -= Plus@@(f /@ ghostProcesses)
      ,
      Print["Ghost sutraction is only handled for process 7"];
      Quit[1]
   ]
];

momReplaceRules = {
   Pair[k[i_Integer], k[j_Integer]] /; j<i :> Pair[k[j], k[i]],
   Pair[k[i_Integer], k[4]] /; i<4 :> Pair[k[i], k[1]] + Pair[k[i], k[2]] - Pair[k[i], k[3]] - Pair[k[i], k[5]],
   Pair[k[4], k[5]] -> Pair[k[1], k[5]] + Pair[k[2], k[5]] - Pair[k[3], k[5]] - Pair[k[5], k[5]],
   Pair[eta[i_Integer], k[4]] :> Pair[eta[i], k[1]] + Pair[eta[i], k[2]] - Pair[eta[i], k[3]] - Pair[eta[i], k[5]],
   Eps[k[i1_], k[i2_], k[i3_], k[4]] :> Eps[k[i1], k[i2], k[i3], k[1]] + Eps[k[i1], k[i2], k[i3], k[2]] - Eps[k[i1], k[i2], k[i3], k[3]] + Eps[k[i1], k[i2], k[i3], k[5]],
   Eps[k[i1_], k[i2_], k[i3_], k[i4_]] :> Signature[{i1, i2, i3, i4}] * Eps@@(k /@ Sort[{i1, i2, i3, i4}]),
   Pair[k[i_Integer], k[i_Integer]] :> mOfLeg[i]^2
};
result = result //. momReplaceRules;

(* We have eliminated one momenta (k[4] in this case).
   Out of 5-1 momenta one can construct only 6 pair.
   We chose: S, T, U, T15, T25 *)
toMandelstam = {
   Pair[k[1], k[5]] -> -(T15 - mOfLeg[1]^2 - mOfLeg[5]^2)/2,
   Pair[k[2], k[5]] -> -(T25 - mOfLeg[2]^2 - mOfLeg[5]^2)/2,
   Pair[k[1], k[2]] ->  (S   - mOfLeg[1]^2 - mOfLeg[2]^2)/2,
   Pair[k[1], k[3]] -> -(T   - mOfLeg[1]^2 - mOfLeg[3]^2)/2,
   Pair[k[2], k[3]] -> -(U   - mOfLeg[2]^2 - mOfLeg[3]^2)/2,
   Pair[k[3], k[5]] ->  (S35 - mOfLeg[3]^2 - mOfLeg[5]^2)/2
};
result = result /. toMandelstam /. MassFu[i_ /; i<3]|MassFd[_] ->0;
result = result /. Den[x_, y_] :> 1/(x-y);

Power[g3, n_] ^:= (alphas*4*Pi)^(n/2);
Conjugate[alphas] ^= alphas;

(* For a 2->3 process there are 5 independent Mandelstam variables.
   For example, one can choose set {S, T, T15, T25, S35} (in FC notation).
   The final expression seems actually shorter after substituting U *)
result = result /. U -> mOfLeg[4]^2 + 2*mOfLeg[3]^2 - S - T - S35 - T15 - T25;
result = result/(alphas*Pi)^3;

(* In FormCalc notation Eps represents -I*eps(....) meaning that Eps is complex
   Since squared 2->3 ME is real, terms with a single Eps have to vanish *)
Print["Trying to simplify Eps tensors..."];
result =
   Simplify[
      result,
      ComplexityFunction -> (100*Count[#, _Eps, {0, Infinity}] + LeafCount[#])&,
      TimeConstraint->Infinity
   ];
If[!FreeQ[result, _Eps],
   Print[result];
   If[FreeQ[result, _Complex], Print["The amplitude is strictly real"]];
   Print["Error: Eps tensor doesn't vanish!"];
   Quit[1],
   Print["Eps tensors vanished"]
];

(* Apparently helicity is averaged, colour and polarization is summed.
   So for gu->sqsq*u we have a factor of
   colour: 1/(3*8)
   polarization: 1/2
   helicity: 2 *)
colAverage =
   Switch[processNr,
      1, 1/9,
      2, 1/8^2,
      3, 1/3^2,
      4, 1/9,
      5, 1/9,
      6, 1/9,
      7, 1/8^2,
      9|10|11|12|13|14, 1/8*1/3,
      _, Print["Unknown color/polarization/helicity average"]; Quit[1]
   ];
polAverage =
   Switch[processNr,
      1, 1,
      2, 1/2^2,
      3, 1,
      4, 1,
      5, 1,
      6, 1,
      7, 1/2^2,
      9|10|11|12|13|14, 1/2,
      _, Print["Unknown color average"]; Quit[1]
   ];
helSum =
   Switch[processNr,
      1|2|3|4|5|6|7, 1,
      9|10|11|12|13|14, 2,
      _, Print["Unknown helicity sum"]; Quit[1]
   ];

result = result*colAverage*helSum*polAverage;

(*
Needs["MultivariateApart`", FileNameJoin[{ParentDirectory @ DirectoryName @ $InputFileName, "MultivariateApart.wl"}]];
result = MultivariateApart[result /. MassSu[_]|MassSd[_]->MassSq];
*)
Needs["Apart`", FileNameJoin[{ParentDirectory @ DirectoryName @ $InputFileName, "Apart.wl"}]];
result = ApartAll[result /. MassSu[_]|MassSd[_]->MassSq /. S35 -> S35p + MassSq^2, {T15, T25, S35p}] /. ApartIR[pcs_,cs_,np_,vars_] :> Times@@((((#1.vars&)/@pcs+cs))^np);
Put[result, "gg_suLsuLdaggerg2.m"];
Quit[1];

RealSquared0 = result //. momReplaceRules;

replaceMandelstam = {
   Pair[k[1], k[3]] -> -(t-TheMass[process[[2,1]]]^2)/2,
   Pair[k[2], k[4]] -> -(t-TheMass[process[[2,2]]]^2)/2,
   Pair[k[2], k[3]] -> -(u-TheMass[process[[2,1]]]^2)/2,
   Pair[k[1], k[4]] -> -(u-TheMass[process[[2,2]]]^2)/2,
   Pair[k[_], k[5]] -> 0,
   Pair[k[1], k[2]] -> s/2,
   Pair[k[3], k[4]] -> (s - TheMass[process[[2,1]]]^2 - TheMass[process[[2,2]]]^2)/2,
   Den[mom2_, mass2_] :> 1/(mom2 - mass2)
};

Needs["MultivariateApart`", FileNameJoin[{ParentDirectory @ DirectoryName @ $InputFileName, "MultivariateApart.wl"}]];

Print[RealSquared0];
Quit[1];

Print[1/s35^2]
s35s35 = FullSimplify[Coefficient[Expand[Apart[RealSquared0]], Den[TheMass[process[[2,1]]]^2 + 2*Pair[k[3], k[5]], TheMass[process[[2,1]]]^2]^2] //. replaceMandelstam];
Print[s35s35];
Print[1/s45^2]
s45s45 = FullSimplify[Coefficient[Expand[Apart[RealSquared0]], Den[TheMass[process[[2,2]]]^2 + 2*Pair[k[4], k[5]], TheMass[process[[2,2]]]^2]^2] //. replaceMandelstam];
Print[s45s45];
Print[1/s35 1/s45]
s35s45 = FullSimplify[
   Coefficient[
   Expand[Apart[RealSquared0]],
   Den[TheMass[process[[2,1]]]^2 + 2*Pair[k[3], k[5]], TheMass[process[[2,1]]]^2]*Den[TheMass[process[[2,2]]]^2 + 2*Pair[k[4], k[5]], TheMass[process[[2,2]]]^2]
   ] //. replaceMandelstam];
Print[s35s45];
Print[1/t15 1/t25]
t15t25 = FullSimplify[Coefficient[Expand[RealSquared0], Den[-2*Pair[k[1], k[5]], 0]*Den[-2*Pair[k[2], k[5]], 0]] //. replaceMandelstam];
Print[t15t25];
Print[1/t15 1/s35]
t15s35 = FullSimplify[Coefficient[Expand[Apart[RealSquared0]], Den[-2*Pair[k[1], k[5]], 0]*Den[TheMass[process[[2,1]]]^2 + 2*Pair[k[3], k[5]], TheMass[process[[2,1]]]^2]] //. replaceMandelstam];
Print[t15s35];
Print[1/t25 1/s35]
t25s35 = FullSimplify[Coefficient[Expand[Apart[RealSquared0]], Den[-2*Pair[k[2], k[5]], 0]*Den[TheMass[process[[2,1]]]^2 + 2*Pair[k[3], k[5]], TheMass[process[[2,1]]]^2]] //. replaceMandelstam];
Print[1/t15 1/s45]
t15s45 = FullSimplify[Coefficient[Expand[Apart[RealSquared0]], Den[-2*Pair[k[1], k[5]], 0]*Den[TheMass[process[[2,2]]]^2 + 2*Pair[k[4], k[5]], TheMass[process[[2,2]]]^2]] //. replaceMandelstam];
Print[t15s45];
Print[1/t25 1/s45]
t25s45 = FullSimplify[Coefficient[Expand[Apart[RealSquared0]], Den[-2*Pair[k[2], k[5]], 0]*Den[TheMass[process[[2,2]]]^2 + 2*Pair[k[4], k[5]], TheMass[process[[2,2]]]^2]] //. replaceMandelstam];
Print[t25s45];
Print[1/t15 1/t15]
t15t15 = FullSimplify[Coefficient[Expand[Apart[RealSquared0]], Den[-2*Pair[k[1], k[5]], 0]^2] //. replaceMandelstam];
Print[t15t15];
Print[1/t25 1/t25]
t25t25 = FullSimplify[Coefficient[Expand[Apart[RealSquared0]], Den[-2*Pair[k[2], k[5]], 0]^2] //. replaceMandelstam];
Print[t25t25];

(*Print[FullSimplify[Normal@Series[t15t25 /. ML2 ->0 /. Dminus4 -> -2*eps, {eps, 0, 1}]/((t^2+u^2)/s^2-eps)  /. t->-s-u]]*)

$Assumptions = eps<0 && b>=0 && b<=1 && dS>0;

RealSquared2 = RealSquared0  //. {
     Pair[k[3], k[5]] -> s35/2, Pair[k[4], k[5]] -> s45/2,
     Pair[k[1], k[5]] -> -t15/2, Pair[k[2], k[5]] -> -t25/2,
     Pair[k[1], k[4]] -> Pair[k[2], k[3]],
     Pair[k[2], k[4]] -> Pair[k[1], k[3]], Pair[k[1], k[2]] -> s12/2,
     Pair[k[3], k[4]] -> 1/2*(s12 - (Plus@@(Power[TheMass[#], 2]& /@ process[[2]]))),
     Pair[k[1], k[3]] -> -t/2, Pair[k[2], k[3]] -> -u/2
     } /. Dminus4 -> -2 eps;

RealSquared0 = RealSquared0 /. Den[p2_, m2_] :> 1/(p2 - m2);
RealSquared0 = RealSquared0 /. {
   S35 -> S35p + TheMass[process[[2,1]]]^2,
   S45 -> S45p + TheMass[process[[2,2]]]^2,
   T14 -> T14p + TheMass[process[[2,2]]]^2,
   T24 -> T24p + TheMass[process[[2,2]]]^2,
   S34 -> S34p + TheMass[process[[2,1]]]^2 + TheMass[process[[2,2]]]^2,
   T -> Tp + TheMass[process[[2,1]]]^2,
   U -> Up + TheMass[process[[2,1]]]^2
};
(*RealSquared0 = Expand@Normal@Series[RealSquared0, {S35p, 0, 0}, {S45p, 0, 0}, {T15, 0, 0}, {T25, 0, 0}];*)

mode = "normal";
(*
   Eq. 3.12
      s35 = Sqrt[s12]E5(1-beta Sin[th] \Sin[th1]Cos[th2]-beta Cos[th]Cos[th1])
      s45 = Sqrt[s12]E5(1+beta Sin[th] \Sin[th1]Cos[th2]+beta Cos[th]Cos[th1]);
   Eq. 3.75
      t15 = -Sqrt[s12] E5(1 - Cos[th1])
      t25 = -Sqrt[s12] E5(1 + Cos[th1]);
      A = 1; C= \[PlusMinus] beta Sin[th]; B = \[PlusMinus] \beta Cos[th]
*)

(* DoubleEnergyIntegral has 1/eps pole,
   SingleEnergyIntegral does not *)
DoubleEnergyIntegral =
   1/Pi*(4/s12)^-eps*Integrate[
      E5^(1-2*eps)/E5^2,
      {E5, 0, dS Sqrt[s12]/2}
   ];

colPolAverage =
   Switch[processNr,
      1, 1/9,
      2, 1/8^2*1/2^2,
      3, 1/3^2,
      4, 1/9,
      5, 1/9,
      6, 1/9,
      (* 7, 1/8^2*1/(2*(1-eps))^2 *)
      7, 1/8^2*1/2^2
   ];

I11n = Pi/(
   a (A + B)) (-(1/eps) +
     Log[(A + B)^2/(
      A^2 - B^2 -
       CCC^2)] - eps (Log[(A - Sqrt[B^2 + CCC^2])/(A + B)]^2 -
        1/2 Log[(A + Sqrt[B^2 + CCC^2])/(A - Sqrt[B^2 + CCC^2])]^2 +
        2 PolyLog[
          2, -((B + Sqrt[B^2 + CCC^2])/(A - Sqrt[B^2 + CCC^2]))] -
        2 PolyLog[2, (B - Sqrt[B^2 + CCC^2])/(A + B)]));

prefactor =(
       (*
       choose how to treat O(eps) prefactors which cancel \
         between real and virtual corrections *)
       Switch[mode,
         "debug",
            Print["keeping all O(eps) prefactors"];
            1/(2 s12) (2^(2 eps)/(
            16 Pi) ((4 Pi mu2)/s12)^eps beta^(
            1 - 2 eps) Sin[th]^(1 - 2 eps)/
            Gamma[1 - eps]) ((4 Pi mu2)/
            s12)^eps Gamma[1 - eps]/
            Gamma[1 - 2 eps],
          "normal",
             Print["excluding all O(eps) prefactors"];
             1/(2*s12)*(beta/(16*Pi)*Sin[th])*(mu2/s12)^eps
       ]*
        1/(2*(2*Pi)^2) *
        (* average over initial states *)
        colPolAverage
        (* final state symmetry factor *)
        * Switch[processNr, 1, 1/2, 2, 1/2, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1])

(* angular integrals *)
replaceAngularIntegrals = {
   (* this substitutions don't have in general beta->1 limit *)
   (* Eq. 3.14 *)
   s35^-2 ->
      1/(2*TheMass[process[[2,1]]]^2) (-1/(2*eps) + Log[dS] - 1/(2*beta)*Log[(1+beta)/(1-beta)]),
   (* Eq. 3.14 *)
   s45^-2 ->
      1/(2*TheMass[process[[2,2]]]^2) (-1/(2*eps) + Log[dS] - 1/(2*beta)*Log[(1+beta)/(1-beta)]),
   (* Eq. 3.15 *)
   s35^-1 s45^-1 ->
      1/(s12*beta)*(-1/(2 eps)*Log[(1+beta)/(1-beta)] - PolyLog[2, (2*beta)/(1+beta)] - 1/4*Log[(1+beta)/(1-beta)]^2 + Log[dS]*Log[(1+beta)/(1-beta)]),
   (* Eq. 3.76 *)
   t15^-1 t25^-1 ->
      1/(2*s12)*(1/eps^2 - 2/eps*Log[dS] + 2*Log[dS]^2),

   (* my own derivations *)
   s35^-1 t15^-1 ->
      DoubleEnergyIntegral/-s12*I11n /. a -> 1 /. A -> 1 /. B -> -Cos[th] beta /. CCC -> -Sin[th] beta,
   s45^-1 t15^-1 -> DoubleEnergyIntegral/-s12*I11n /. a -> 1 /. A -> 1 /. B -> Cos[th] beta /. CCC -> Sin[th] beta,
   (* is there really a t/
      u symmetry and does it work like I do it below? *)

   s35^-1 t25^-1 ->
      DoubleEnergyIntegral/-s12*I11n /. a -> 1 /. A -> 1 /. B ->  Cos[th] beta /. CCC ->  Sin[th] beta,
   s45^-1 t25^-1 ->
      DoubleEnergyIntegral/-s12*I11n /. a -> 1 /. A -> 1 /. B -> -Cos[th] beta /. CCC -> -Sin[th] beta,

   (* this integrals are only needed for gg initial state *)
   t15^-2 ->
        -DoubleEnergyIntegral/s12 * Pi * (1 - eps)*0,
   t25^-2 ->
        -DoubleEnergyIntegral/s12 * Pi * (1 - eps)*0,

   (* 1/(single Mandelstam variable) integrals go to 0 as dS->0 *)
   s35^-1 | s45^-1 | t15^-1 | t25^-1 -> 0
};

res = s35s35 * s35^-2 + s45s45 * s45^-2 + s35s45 * s35^-1 s45^-1 + t15t25 * t15^-1 t25^-1 + t15s35 * s35^-1 t15^-1 + t15s45 * s45^-1 t15^-1 + t25s35 * s35^-1 t25^-1 + t25s45 * s45^-1 t25^-1 + t15t15 * t15^-2*0 + t25t25 * t25^-2 *0//. replaceAngularIntegrals;

(* see https://reference.wolfram.com/language/Developer/ref/PolyLogSimplify.html *)
Needs["Developer`"];
(* PolyLogSimplify should be auto usec by fullsimplify 
res = PolyLogSimplify[res];
*)

RealSquared2 = If[FreeQ[#, Alternatives @@ {s35, s45, t15, t25}], 0, #]& /@ RealSquared2;
(*Print["Ania"];
(Print[#];Print[# /. replaceAngularIntegrals]; Print@FreeQ[#, Alternatives @@ {s35, s45, t15, t25}])& /@ RealSquared2;*)

RealSquared3 = RealSquared2 /. replaceAngularIntegrals /. s34 -> s12 /.
   t -> -(s12/2) (1 + TheMass[process[[2,1]]]^2/s12 - TheMass[process[[2,2]]]^2/s12 - beta*Cos[th]) + TheMass[process[[2,1]]]^2 /.
   u -> -(s12/2) (1 + TheMass[process[[2,2]]]^2/s12 - TheMass[process[[2,1]]]^2/s12 + beta*Cos[th]) + TheMass[process[[2,2]]]^2 ;

(* O(eps^0) term *)
res0ord = SeriesCoefficient[Series[prefactor*res /. Dminus4 -> -2*eps, {eps, 0, 0}], 0];
res0ord = res0ord /. t -> -(s12/2) (1 + TheMass[process[[2,1]]]^2/s12 - TheMass[process[[2,2]]]^2/s12 - beta*Cos[th])  + TheMass[process[[2,1]]]^2 /.
   u -> -(s12/2) (1 + TheMass[process[[2,2]]]^2/s12 - TheMass[process[[2,1]]]^2/s12 + beta*Cos[th]) + TheMass[process[[2,1]]]^2 /. TheMass[process[[2,1]]]|TheMass[process[[2,2]]] -> 1/2*Sqrt[(1-b^2)*s] ;

replace = {Pi->pi, s12->s, Sqrt[x_]->sqrt[x], ArcTanh[x_]:>atanh[x], Log[x_]:>log[x]};


res0ord = FullSimplify[res0ord];

Unprotect[Power];
Power[MassGlu, 2] ^= MassGlu2;
Power[MassSq, 2] ^= MassSq2;
Power[MassTop, 2] ^= MassTop2;
Power[mu, 2] ^= mu2;
Power[mu, -2] ^= 1/mu2;
Power[alphas, 3] ^= alphas3;
Power[MasssigmaO, 2] ^= MasssigmaO2;
Format[Power[x_, i_Integer], CForm] := "pow<i>(x)";
Protect[Power];

(* write final expression to a file *)
WriteString[
   FileNameJoin[{"src", "models", model, ToLowerCase[model] <> "_" <> fileName <> "_soft.cpp"}],
"double " <> model <> "::matrixSoft_" <> fileName <> "(const double alphas, const double s, const double th, const double dS, const double mu) const {
   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/s);
   const double lndS = std::log(dS);
   const double cth = std::cos(th);
   const double sth = std::sin(th);\n" <>
   If[!FreeQ[res0ord, MassGlu2], "   const double MassGlu2 = pow<2>(MassGlu);\n", ""] <>
   If[!FreeQ[res0ord, MassSq2],  "   const double MassSq2 = pow<2>(MassSq);\n", ""] <>
   If[!FreeQ[res0ord, MassTop2],  "   const double MassTop2 = pow<2>(MassTop);\n", ""] <>
   If[!FreeQ[res0ord, MasssigmaO2],  "   const double MasssigmaO2 = pow<2>(MasssigmaO);\n", ""] <>
   If[!FreeQ[res0ord, mu2],  "   const double mu2 = pow<2>(mu);\n", ""] <>
   If[!FreeQ[res0ord, alphas3],  "   const double alphas3 = pow<3>(alphas);\n", ""] <>
"   const double res = " <>
   StringReplace[
      ToString@CForm[
         N[
            res0ord /.
               Log[dS] -> lndS /.
               PolyLog[2, x_] :> Li2[x] /.
               Cos[th] -> cth /. Sin[th] -> sth /.
               beta->b /. th -> th /. Cos->cos //. replace,
            18
         ]
      ],
      "Li2" -> "polylogarithm::Li2"
   ] <> ";\n" <>
"   return res;\n}"
];

