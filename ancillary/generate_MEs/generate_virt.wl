<<FeynArts.m;
<<FormCalc.m;

(* default is 1 for one-loop which is fine for UV divergences but not for IR
   ones *)
$Dminus4MaxPower = 2;

(* topologies *)
treeTop = CreateTopologies[0, 2->2];
virtTop = CreateTopologies[1, 2->2, ExcludeTopologies->{Tadpoles,   WFCorrections}];
ctTop = CreateCTTopologies[1, 2->2, ExcludeTopologies->{TadpoleCTs, WFCorrectionCTs}];

(* process *)
model = $ScriptCommandLine[[2]];
processNr = ToExpression[$ScriptCommandLine[[3]]];

SetOptions[InsertFields, Model -> model <> "EWSB_myRenConst", ExcludeParticles -> {V[n_ /; n != 5]}];

Switch[processNr,
   3, process = {F[3, {1}], F[3, {1}]} ->  {S[13, {1}], S[13, {4}]}; fileName = "uu_suLsuR",
   5, process = {F[3, {1}], -F[3, {1}]} -> {S[13,{1}], -S[13,{1}]};  fileName = "uubar_suLsuLdagger",
   6, process = {F[4, {1}], -F[4, {1}]} -> {S[13,{1}], -S[13,{1}]};  fileName = "ddbar_suLsuLdagger",
   7, process = {V[5], V[5]} -> {S[13,{1}], -S[13,{1}]}; fileName = "gg_suLsuLdagger",
   8, process = {-F[3, {1}], -F[3, {1}]} ->  {-S[13, {1}], -S[13, {4}]}; fileName = "ubarubar_suLdaggersuRdagger",
   9, process = {F[3, {1}], F[4, {1}]} ->  {S[13, {1}], S[14, {4}]}; fileName = "ud_suLsdR",
   _, Print["Unknown process ", processNr]
];

colPolAverage =
   Switch[processNr,
      1, 1/9,
      2, 1/8^2 1/2^2,
      3|5|6|8|9, 1/3^2,
      (* 7, 1/8^2*1/(Dminus4+2)^2, *)
      7, 1/8^2*1/2^2,
      _, Print["Unknown color/polarization averaging factor"]; Quit[1]
   ];

finalStateSymmetry =
   Switch[processNr,
      1|2, 1/2,
      3|4|5|6|7|8|9, 1,
      _, Print["Unknown final state symmetry factor"]; Quit[1]
   ];

Neglect[MassFu[1]] = 0;
Neglect[MassFu[2]] = 0;
Neglect[MassFd[_]] = 0;
Neglect[g1] = 0;
Neglect[g2] = 0;
Neglect[Yu[__]] = 0;
Neglect[Yd[__]] = 0;

(* insertions *)
bornIns = InsertFields[treeTop, process];
ctIns = InsertFields[ctTop, process];
virtIns = InsertFields[virtTop, process];

(* keep O(eps) terms from everything *)
SetOptions[CalcFeynAmp, FermionChains->VA, Dimension->0, FermionOrder->None];
SetOptions[HelicityME, Dimension->0];
SetOptions[PolarizationSum, Dimension->4, GaugeTerms->True];

(* amplitudes *)
born = CalcFeynAmp[CreateFeynAmp[bornIns]];
virt = CalcFeynAmp[CreateFeynAmp[virtIns]];
counter = CalcFeynAmp[CreateFeynAmp[ctIns]];

mandelstamReplaceRule = U -> Plus@@(Power[TheMass[#], 2]& /@ process[[2]]) - S - T
momReplaceRules = {
   Pair[eta[i_Integer], k[4]] :> Pair[eta[i], k[1]] + Pair[eta[i], k[2]] - Pair[eta[i], k[3]],
   Eps[k[i1_], k[i2_], k[i3_], k[4]] :> Eps[k[i1], k[i2], k[i3], k[1]] + Eps[k[i1], k[i2], k[i3], k[2]] - Eps[k[i1], k[i2], k[i3], k[3]],
   Eps[k[i1_], k[i2_], k[i3_], k[i4_]] :> Signature[{i1, i2, i3, i4}] * Eps@@(k /@ Sort[{i1, i2, i3, i4}])
};

_Hel=0;

Print["Born..."];
colour   = ColourME[born];
helicity = HelicityME[born];
born2    = SquaredME[born];
(* needed for hard-colinear counterterm *)
born2 =
   FullSimplify[
      FullSimplify[PolarizationSum[born2[[1]] /. born2[[2]] /. helicity] /. colour //. SubExpr[] //. Abbr[] //. momReplaceRules /. Den[x_,y_]:>1/(x-y), Assumptions-> U==Plus@@(Power[TheMass[#], 2]& /@ process[[2]]) - S - T]
   ];

Print["Un-renormalized virtuals..."];
colour                 = ColourME[virt, born];
helicity               = HelicityME[virt, born];
unrenormalizedVirtBorn = SquaredME[born, virt];
unrenormalizedVirtBorn = PolarizationSum[unrenormalizedVirtBorn[[1]] /. unrenormalizedVirtBorn[[2]] /. helicity] /. colour //. SubExpr[] //. Abbr[];

Print["CTs..."];
colour   = ColourME[born, counter];
helicity = HelicityME[born, counter];
ctBorn   = SquaredME[born, counter];
ctBorn   = PolarizationSum[ctBorn[[1]] /. ctBorn[[2]] /. helicity] /. colour //. SubExpr[] //. Abbr[];

renConstSubs = (# -> RenConst[#])& /@ {
   dZSu1[1], dZGlL1, dZgs1, dMGl1, dZGG1, dMSu1[1], dZFu1[1], dZFd1[1], dZSu1[4], dZGlR1
};
ctBorn = ctBorn //. renConstSubs;

total = ExpandSums[2*(unrenormalizedVirtBorn + ctBorn)];

total = total /. MassFu[1|2]|MassFd[_]|Yu[_,_]|Yd[_,_]|g1|g2 ->0 /. Den[x_,y_]:>1/(x-y) /. MDO->MassGlu //. momReplaceRules;

WriteString["stdout", "Checking UV finitenes... "];
uvDiv =
   FullSimplify[
      UVDivergentPart[total /. Dminus4->0 /. mandelstamReplaceRule]
   ];
If[uvDiv =!= 0,
   Print["fail\n", uvDiv];
   Quit[1];,
   WriteString["stdout", "ok\n"];
];

If[!FreeQ[total, _eta],
   WriteString["stdout", "Checking gauge invariance... "];
   (* 7 is *probably* to complicated to prove cancellation or eta terms *)
   If[processNr == 7,
      total = total /. eta[1] -> k[2] /. eta[2] -> k[1] /. Pair[k[i_], k[j_]] /; i>j :> Pair[k[j], k[i]] /. Pair[k[2], k[2]]|Pair[k[1],k[1]] -> 0 /. Pair[k[1], k[3]] -> -1/2*(T - TheMass[process[[2,1]]]^2) /. Pair[k[2], k[3]] -> -1/2*(U - TheMass[process[[2,1]]]^2) /. Pair[k[1], k[2]] -> S/2 /. Pair[k[1], k[4]] -> -1/2*(U-TheMass[process[[2,2]]]^2) /. Pair[k[2], k[4]] -> -1/2*(T-TheMass[process[[2,1]]]^2) //. momReplaceRules;
   ];
   (*
   total =
      Simplify[
         total,
         ComplexityFunction -> Count[#, _eta, {0, Infinity}] + LeafCount[#]&,
         Assumptions -> U == Plus@@(Power[TheMass[#], 2]& /@ process[[2]]) - S - T,
         TimeConstraint->Infinity
      ];
      *)

   If[!FreeQ[total, _eta],
      Print["fail"],
      Print["ok"]
   ];
];

If[!FreeQ[result, _Eps],
   Print["Trying to simplify Eps tensors..."];
   result =
      Simplify[
         result,
         ComplexityFunction -> (100*Count[#, _Eps, {0, Infinity}] + LeafCount[#])&,
         TimeConstraint->Infinity
   ];
];

(* symmetry factors *)
total = total*colPolAverage*finalStateSymmetry;

(* WRITING RESULTS TO FILES *)

Unprotect[Power];
Format[Power[x_, i_Integer], CForm] := Format["pow<" <> ToString[i] <> ">(" <> ToString@CForm[x] <> ")", OutputForm];
Format[Power[x_, i_Integer], CForm] := Format["pow<" <> ToString[i] <> ">(" <> ToString@CForm[x] <> ")", OutputForm];
squareAbbreviation = {MassGlu, MassSq, MassTop, MasssigmaO, MassphiO, mu};
(Format[Power[#, 2], CForm] := Format[ToString[#] <> "2", OutputForm])& /@ squareAbbreviation;
Format[Power[mu, -2], CForm] := Format["1/mu2", OutputForm];
Format[Power[x_ /; MemberQ[squareAbbreviation, x], i_Integer /; EvenQ[i]], CForm] := Format["pow<" <> ToString[i/2] <> ">(" <> ToString@CForm[x^2] <> ")", OutputForm];
Protect[Power];

Power[g3, n_] ^:= (alphas*4*Pi)^(n/2);
Conjugate[alphas] ^= alphas;

abbrFinalStateMasses = {TheMass[process[[2,1]]]->MassSq, TheMass[process[[2,2]]]->MassSq};
mandelstamReplaceRule = mandelstamReplaceRule /. abbrFinalStateMasses;
born2 = born2 /. abbrFinalStateMasses;
(*
flux = 1/(2*S);
phaseSpaceDminus4pow0 = 1/(8*S*Pi);
*)
(*phaseSpaceDminus4pow1 = Pi^(-Dminus4/2)/2^Dminus4 *1/Gamma[1 + Dminus4/2]*((T*U-TheMass[process[[2,1]]]^2*TheMass[process[[2,2]]]^2)/S)^(Dminus4/2)*)
phaseSpaceDminus4pow1 = 1;

born2=born2*phaseSpaceDminus4pow1*colPolAverage*finalStateSymmetry;
dDimBorn = FullSimplify@Series[born2, {Dminus4, 0, 2}];

total = phaseSpaceDminus4pow1*total /. abbrFinalStateMasses /. MassSu[_]|MassSd[_]->MassSq /. MassFu[3] -> MassTop /. MDO -> MassGlu;

(* Re will not simplify expressions like
   Re[1/(-MassGlu^2 + T)] unless it knows that T!=MassGlu^2 *)
totalRef =
   Refine[
      Re[total],
      Assumptions->And@@(Element[#, Reals]& /@ {S,T,U,Dminus4,MassGlu,MassSq,MassTop,FiniteGs,alphas,Divergence})&& And@@(#>0& /@ {MassSq, MasssigmaO, MassphiO, MassTop, MassGlu, mu})&&T!=MassGlu^2&&U!=MassGlu^2&&S>0&&Dminus4!=-2&&S!=MassphiO^2&&T!=MassSq^2&&U!=MassSq^2&&T*U>MassSq^4
   ];

Print["Expanding B*V to order Dminus4^", $Dminus4MaxPower];
(* Include Dminus4 offset of 2-body phase space.
   4d phase space will be included in the integration routine *)
totalSeries = Series[totalRef, {Dminus4, 0, $Dminus4MaxPower}];
If[SeriesCoefficient[total, 2] =!= 0,
   Print["INFO: There are O(Dminus4^2) terms"]
];

Print["Simplifying final expression..."];

(* WRITE BORN ME TO ORD(eps^2) *)
filePathBorn = FileNameJoin[{"src", "models", model, ToLowerCase[model] <> "_" <> fileName <> "_born_matrix.cpp"}];
Print["Writing ME to " <> filePathBorn];
WriteString[
   OpenWrite[filePathBorn, CharacterEncoding -> "UTF8"],
"double " <> model <> "::matrixTree_" <> fileName <> "(const double alphas, const double S, const double T, const int Dminus4SeriesCoeff) const {\n" <>
   If[!FreeQ[dDimBorn, MassGlu^n_ /; EvenQ[n]], "   const double MassGlu2 = pow<2>(MassGlu);\n", ""] <>
   If[!FreeQ[{dDimBorn, Plus@@(Power[TheMass[#], 2]& /@ process[[2]]) /. abbrFinalStateMasses}, MassSq^n_ /; EvenQ[n]],  "   const double MassSq2 = pow<2>(MassSq);\n", ""] <>
   If[!FreeQ[dDimBorn, MassTop^n_ /; EvenQ[n]],  "   const double MassTop2 = pow<2>(MassTop);\n", ""] <>
   If[!FreeQ[dDimBorn, MasssigmaO^n_ /; EvenQ[n]],  "   const double MasssigmaO2 = pow<2>(MasssigmaO);\n", ""] <>
   If[!FreeQ[dDimBorn, MassphiO^n_ /; EvenQ[n]],  "   const double MassphiO2 = pow<2>(MassphiO);\n", ""] <>
   If[!FreeQ[dDimBorn, U],  "   const double U = 2.*MassSq2 - T - S;\n", ""] <>
"   double res = 0.;
   if (Dminus4SeriesCoeff == 0) {
      res = " <> StringReplace[ToString[CForm[If[FreeQ[born2, Dminus4], born2, SeriesCoefficient[dDimBorn, 0]]]], {"Power"->"pow", "Log"->"log", "Pi"->"pi", "EulerGamma"->"euler"}] <> ";\n" <>
"   }
   else if (Dminus4SeriesCoeff == 1) {
      res = " <> StringReplace[ToString[CForm[If[FreeQ[born2, Dminus4], 0, SeriesCoefficient[dDimBorn, 1]]]], {"Power"->"pow", "Log"->"log", "Pi"->"pi", "EulerGamma"->"euler"}] <> ";\n" <>
"   }
   else if (Dminus4SeriesCoeff == 2) {
      res = " <> StringReplace[ToString[CForm[If[FreeQ[born2, Dminus4], 0, SeriesCoefficient[dDimBorn, 2]]]], {"Power"->"pow", "Log"->"log", "Pi"->"pi", "EulerGamma"->"euler"}] <> ";\n" <>
"   }
   if(std::isinf(res)) {
      const double beta = std::abs(1- 4*MassSq2/S);
      if (beta < 1e-10) {
         spdlog::get(\"console\")->warn(\"Infinity in the 2->2 |ME|^2 at threshold (\[Beta] = {}). Returning 0.\", beta);
         return 0.;
      }
      else {
         spdlog::get(\"console\")->error(\"Infinity in the 2->2 |ME|^2 for \[Beta] = {}\", beta);
         return res;
      }
   }
   else {
      return res;
   }\n}\n"
];

(* WRITE VIRT ME TO ORD(eps^2) *)
Print["Dminus4^0"];
total0 =
   FullSimplify[
      SeriesCoefficient[totalSeries, 0],
      Assumptions -> U == 2*MassSq^2 - S - T,
      TimeConstraint->Infinity
   ];
Print["Dminus4^1"];
totalDm4 =
   FullSimplify[
      SeriesCoefficient[totalSeries, 1],
      Assumptions -> U == 2*MassSq^2 - S - T,
      TimeConstraint->Infinity
   ];
Print["Dminus4^2"];
totalDm4Sq =
   FullSimplify[
      SeriesCoefficient[totalSeries, 2],
      Assumptions -> U == 2*MassSq^2 - S - T,
      TimeConstraint->Infinity
   ];

filePath = FileNameJoin[{"src", "models", model, ToLowerCase[model] <> "_" <> fileName <> "_virt_matrix.cpp"}];
Print["Writing ME to " <> filePath];
WriteString[
   filePath,
"double " <> model <> "::matrixVirt_" <> fileName <> "(const double alphas, const double S, const double T,
   const double FiniteGs, const int Dminus4SeriesCoeff, const int lambda, const double mu) const {\n" <>
   If[!FreeQ[total, MassGlu^2], "   const double MassGlu2 = pow<2>(MassGlu);\n", ""] <>
   If[!FreeQ[total, MassSq^2],  "   const double MassSq2 = pow<2>(MassSq);\n", ""] <>
   If[!FreeQ[total, MassTop^2],  "   const double MassTop2 = pow<2>(MassTop);\n", ""] <>
   If[!FreeQ[total, MasssigmaO^2],  "   const double MasssigmaO2 = pow<2>(MasssigmaO);\n", ""] <>
   If[!FreeQ[total, MassphiO^2],  "   const double MassphiO2 = pow<2>(MassphiO);\n", ""] <>
"   const double mu2 = pow<2>(mu);
   const double U = 2*MassSq2 - S - T;
   setmudim(mu2);
   setlambda(lambda);\n" <>
   If[!FreeQ[total0, Divergence] || !FreeQ[totalDm4, Divergence], "   const double Divergence = 1 ? lambda == -1 : 0;\n", ""] <>
"   if (Dminus4SeriesCoeff == 0) {
      return " <> StringReplace[ToString[CForm[total0 /. Pi->pi]], {"Power"->"pow", "Log"->"log", "EulerGamma"->"euler"}] <> ";\n" <>
"   }
   else if (Dminus4SeriesCoeff == 1) {
      return " <> StringReplace[ToString[CForm[totalDm4 /. Pi->pi]], {"Power"->"pow", "Log"->"log", "EulerGamma"->"euler"}] <> ";\n" <>
"   }
   else if (Dminus4SeriesCoeff == 2) {
      return " <> StringReplace[ToString[CForm[totalDm4Sq /. Pi->pi]], {"Power"->"pow", "Log"->"log", "EulerGamma"->"euler"}] <> ";\n" <>
"   }\n}\n"
];

(* VALIDATION *)

Install[FileNameJoin[{$HomeDirectory, "HEP-software/gcc/LoopTools-2.16/bin/LoopTools"}]];

numerics = {
  MassSu[_] -> 1500, MassSd[_] -> 1500, MassSq -> 1500,
  MassGlu -> 1000, MDO -> 1000,
  S -> 6000^2, T -> -22208172,
  mu->1500,
  MassTop->172,
  MasssigmaO->5000,
  MassphiO->Sqrt[5000^2 + 4*1000^2]
};

Block[{U = 2*MassSq^2 - S - T},
normalization = born2*alphas/(2*Pi) /. numerics /. FiniteGs|Dminus4->0;

SetMudim[1500^2];
SetUVDiv[1];
SetLambda[-2];
dpole = Simplify[total0/normalization /. numerics /. FiniteGs->0];
Print["tree: ", N[born2 /. FiniteGs|Dminus4->0 /. numerics] /. alphas->0.1184];
Print["double pole: ", dpole];
dtospole = totalDm4*(-2)/normalization /. numerics /. FiniteGs|Divergence->0;

SetLambda[-1];
spole = Expand[(total0/normalization /. numerics /. FiniteGs->0 /. Divergence->1) + A*dtospole];
Print["single pole: ", spole];
stofinite = Expand[totalDm4/normalization /. numerics /. FiniteGs -> 0 /. Divergence->1]*(-2);

SetLambda[0];
finite = Expand[(total0/normalization /. numerics /. FiniteGs -> 1 /. Divergence->0) + dpole*Pi^2/6*0 + stofinite];
Print["finite: ", finite];
];
