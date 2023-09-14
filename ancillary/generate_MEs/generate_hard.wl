<<FeynArts`
<<FormCalc`

If[Length[$ScriptCommandLine] =!= 4,
   Print["Need exactly 3 arguments. Got ", $ScriptCommandLine];
   Quit[1]
];

model = $ScriptCommandLine[[2]];
processNr = ToExpression[$ScriptCommandLine[[3]]];
diagRemoval =
   Switch[$ScriptCommandLine[[4]],
      "NoDR", False,
      "DR", True,
      _, Print["Diagram removal yes/no"]; Quit[1]
   ];

Switch[model,
   (* exclude W, Z, gamma *)
   "MSSM"|"MRSSM", SetOptions[InsertFields, Model -> model <> "EWSB_myRenConst", InsertionLevel -> {Particles}, ExcludeParticles -> {V[n_ /; n != 5]}],
   (* exclude W, Z *)
   "SMQEDQCD", SetOptions[InsertFields, Model -> "SMQCD", InsertionLevel -> {Particles}, ExcludeParticles -> {V[2], V[3]}],
   "SMQCD", SetOptions[InsertFields, Model -> "SMQCD", InsertionLevel -> {Particles}, ExcludeParticles -> {V[n_ /; n != 5]}],
   (*
SetOptions[InsertFields, Model -> "MyMRSSM",
         ExcludeParticles -> {S[n_ /; n != 26],(*V[n_/;n\[NotEqual]5],*)
          F[n_ /; n == 15]}];
      SetOptions[InsertFields, Model -> "MyMRSSM"];
      *)
   _, Print["Unknown model ", model]; Quit[1]
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
If[diagRemoval,
   ins = DiagramSelect[ins, (Not[SChannelQ[F[15]][##]] && Not[SChannelQ[-F[15]][##]])&];
];

Neglect[MassFu[1]]=0;
Neglect[MassFu[2]]=0;
Neglect[MassFd[_]]=0;

mOfLeg[i_Integer] := TheMass[Flatten[List@@process][[i]]];
Square[mOfLeg[3]];
Square[mOfLeg[4]];

SetOptions[CalcFeynAmp, FermionChains->Chiral];
If[diagRemoval,
   SetOptions[PolarizationSum, GaugeTerms->True],
   SetOptions[PolarizationSum, GaugeTerms->Off]
];

(* For gg->suLsuL*g with Transverse->False FormCalc crashes because the
   expression is larger than a hardcoded maximal size. One needs to edit file
   ReadForm.tm in FormCalc directory and change TERMBUF to lets say 6000000 *)
amp =
   If[diagRemoval,
      CalcFeynAmp[CreateFeynAmp[ins], Invariants->False, Normalized->True, Transverse->True],
      CalcFeynAmp[CreateFeynAmp[ins], Invariants->False, Normalized->False, Transverse->False]
   ];
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
   Eps[k[i1_], k[i2_], k[i3_], k[4]] :> Eps[k[i1], k[i2], k[i3], k[1]] + Eps[k[i1], k[i2], k[i3], k[2]] - Eps[k[i1], k[i2], k[i3], k[3]] - Eps[k[i1], k[i2], k[i3], k[5]],
   Eps[k[i1_], k[i2_], k[4], k[5]] :> Eps[k[i1], k[i2], k[3], k[5]] + Eps[k[i1], k[i2], k[4], k[5]] - Eps[k[i1], k[i2], k[1], k[5]] - Eps[k[i1], k[i2], k[2], k[5]],
   Eps[k[i1_], k[i2_], k[i3_], k[i4_]] :> Signature[{i1, i2, i3, i4}] * Eps@@(k /@ Sort[{i1, i2, i3, i4}]),
   Pair[k[i_Integer], k[i_Integer]] :> mOfLeg[i]^2
};
If[diagRemoval,
   result = result /. eta[1] -> k[2]
];
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
If[!FreeQ[result, _Eps],
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

Needs["MultivariateApart`", FileNameJoin[{ParentDirectory @ DirectoryName @ $InputFileName, "MultivariateApart.wl"}]];
result = MultivariateApart[result /. MassSu[_]|MassSd[_]->MassSq];

Unprotect[Power];
Format[Power[x_, i_Integer], CForm] := Format["pow<" <> ToString[i] <> ">(" <> ToString@CForm[x] <> ")", OutputForm];
Format[Power[x_, i_Integer], CForm] := Format["pow<" <> ToString[i] <> ">(" <> ToString@CForm[x] <> ")", OutputForm];
squareAbbreviation = {MassGlu, MassSq, MassTop, MasssigmaO, MassphiO, mu};
(Format[Power[#, 2], CForm] := Format[ToString[#] <> "2", OutputForm])& /@ squareAbbreviation;
Format[Power[x_ /; MemberQ[squareAbbreviation, x], i_Integer /; EvenQ[i]], CForm] := Format["pow<" <> ToString[i/2] <> ">(" <> ToString@CForm[x^2] <> ")", OutputForm];
Protect[Power];

filePath = FileNameJoin[{"src", "models", model, ToLowerCase[model] <> "_" <> procName <> "_hard" <> If[diagRemoval, "-DR", ""] <> ".cpp"}];
Print["Writing to " <> filePath];

(* write final expression to a file *)
WriteString[
   filePath,
"double " <> model <> "::matrixHard_" <> procName <> If[diagRemoval, "_DR", ""] <> "(const double alphas, std::array<std::array<double, 4>, 5> const &p) const {
   const double k12 = p[0][0]*p[1][0]-p[0][3]*p[1][3];
   const double k13 = p[0][0]*p[2][0]-p[0][3]*p[2][3];
   const double k35 = p[2][0]*p[4][0]-p[2][1]*p[4][1]-p[2][2]*p[4][2]-p[2][3]*p[4][3];
   const double T15 = -2*(p[0][0]*p[4][0]-p[0][3]*p[4][3]);
   const double T25 = -2*(p[1][0]*p[4][0]-p[1][3]*p[4][3]);
   const double S   = 2*k12;
   const double T   = pow<2>(MassSq) - 2*k13;
   const double S35 = pow<2>(MassSq) + 2*k35;\n" <>
   If[!FreeQ[result, MassGlu^2], "   const double MassGlu2 = pow<2>(MassGlu);\n", ""] <>
   If[!FreeQ[result, MassSq^2],  "   const double MassSq2 = pow<2>(MassSq);\n", ""] <>
"   const double res = " <>
      ToString@CForm[result] <> ";\n" <>
"   return pow<3>(alphas*pi)*res;\n}\n"
];
