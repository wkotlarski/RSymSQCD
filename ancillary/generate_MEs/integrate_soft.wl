expr = Get[$ScriptCommandLine[[2]]];

findPower[expr_]:=Module[{nom, denom, x},
   {nom, denom} = NumeratorDenominator[expr] /. -S35p - T15 - T25 -> -x /. S35p + T15 + T25 -> x /. S35p -> x /. T15 -> x /. T25 ->x;
   If[FreeQ[#, x], 0, (List@@Series[#, {x,0,0}])[[-3]]]& /@ {nom, denom}
];

expr = If[findPower[#][[1]] == 0 && findPower[#][[2]]==2, #, 0]& /@ expr;

$Assumptions = eps<0 && beta>=0 && beta<=1 && dS>0 && th>=0 && th<=Pi && S>0 && mu2>0;
processNr = 7;
(* DoubleEnergyIntegral has 1/eps pole,
   SingleEnergyIntegral does not *)
DoubleEnergyIntegral =
   1/Pi*(4/S)^-eps*Integrate[
      E5^(1-2*eps)/E5^2,
      {E5, 0, dS Sqrt[S]/2}
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

flux = 1/(2*S);
bornPS = 2^(2*eps)/(16*Pi) * ((4*Pi)/S)^eps * beta^(1-2*eps) * Sin[th]^(1-2*eps)/Gamma[1 - eps];
threeBody = ((4*Pi*mu2)/S)^eps*Gamma[1 - eps]/Gamma[1 - 2 eps];

prefactor = (
    (flux * bornPS * threeBody /. eps->0)*0+
     (* flux*(beta/(16*Pi)*Sin[th])* *) (mu2/S)^eps *
        1/(2*(2*Pi)^2) *
        (* final state symmetry factor *)
        Switch[processNr, 1, 1/2, 2, 1/2, 3, 1, 4, 1, 5, 1, 6, 1, 7, 1]
)

(* angular integrals *)
replaceAngularIntegrals = {
   (* this substitutions don't have in general beta->1 limit *)
   (* Eq. 3.14 *)
   (* Eq. 3.14 *)
   (-S35p - T15 - T25)^-2|(S35p + T15 + T25)^-2 ->
      1/(2*MassSq^2) (-1/(2*eps) + Log[dS] - 1/(2*beta)*Log[(1+beta)/(1-beta)]),
   S35p^-2 ->
      1/(2*MassSq^2) (-1/(2*eps) + Log[dS] - 1/(2*beta)*Log[(1+beta)/(1-beta)]),
   (* Eq. 3.15 *)
   (-S35p - T15 - T25)^-1 S35p^-1 ->
      1/(S*beta)*(-1/(2 eps)*Log[(1+beta)/(1-beta)] - PolyLog[2, (2*beta)/(1+beta)] - 1/4*Log[(1+beta)/(1-beta)]^2 + Log[dS]*Log[(1+beta)/(1-beta)]),
   (S35p + T15 + T25)^-1 S35p^-1 ->
      -1/(S*beta)*(-1/(2 eps)*Log[(1+beta)/(1-beta)] - PolyLog[2, (2*beta)/(1+beta)] - 1/4*Log[(1+beta)/(1-beta)]^2 + Log[dS]*Log[(1+beta)/(1-beta)]),
   (* Eq. 3.76 *)
   T15^-1 T25^-1 ->
      1/(2*S)*(1/eps^2 - 2/eps*Log[dS] + 2*Log[dS]^2),

   (* my own derivations *)
   S35p^-1 T15^-1 ->
      DoubleEnergyIntegral/-S*I11n /. a -> 1 /. A -> 1 /. B -> -Cos[th] beta /. CCC -> -Sin[th] beta,
   (-S35p - T15 - T25)^-1 T15^-1 -> DoubleEnergyIntegral/-S*I11n /. a -> 1 /. A -> 1 /. B -> Cos[th] beta /. CCC -> Sin[th] beta,
   (S35p + T15 + T25)^-1 T15^-1 -> -DoubleEnergyIntegral/-S*I11n /. a -> 1 /. A -> 1 /. B -> Cos[th] beta /. CCC -> Sin[th] beta,
   (* is there really a t/
      u symmetry and does it work like I do it below? *)

   S35p^-1 T25^-1 ->
      DoubleEnergyIntegral/-S*I11n /. a -> 1 /. A -> 1 /. B ->  Cos[th] beta /. CCC ->  Sin[th] beta,
   (-S35p - T15 - T25)^-1 T25^-1 ->
      DoubleEnergyIntegral/-S*I11n /. a -> 1 /. A -> 1 /. B -> -Cos[th] beta /. CCC -> -Sin[th] beta,
   T25^-1 * (S35p + T15 + T25)^-1  ->
      -DoubleEnergyIntegral/-S*I11n /. a -> 1 /. A -> 1 /. B -> -Cos[th] beta /. CCC -> -Sin[th] beta,

   (* this integrals are only needed for gg initial state *)
   T15^-2 ->
        -DoubleEnergyIntegral/S * Pi * (1 - eps)*0,
   T25^-2 ->
        -DoubleEnergyIntegral/S * Pi * (1 - eps)*0
};

expr = (expr /. replaceAngularIntegrals /. T15|T25|S35p->0);

expr = expr /.
   U->2*MassSq^2-S-T /.
   T -> -(S/2) (1 + MassSq^2/S - MassSq^2/S - beta*Cos[th]) + MassSq^2 /.
   MassSq->Sqrt[S - b^2*S]/2;

(* O(eps^0) term *)
resAllOrd = Series[prefactor*expr /. Dminus4 -> -2*eps, {eps, 0, 0}];

replace = {Pi->pi, Sqrt[x_]->sqrt[x], ArcTanh[x_]:>atanh[x], Log[x_]:>log[x]};

Needs["Developer`"];
resAllOrd = PolyLogSimplify[resAllOrd];
resAllOrd = FullSimplify[resAllOrd];

Unprotect[Power];
Format[Power[x_, i_Integer], CForm] := Format["pow<" <> ToString[i] <> ">(" <> ToString@CForm[x] <> ")", OutputForm];
Format[Power[x_, i_Integer], CForm] := Format["pow<" <> ToString[i] <> ">(" <> ToString@CForm[x] <> ")", OutputForm];
squareAbbreviation = {MassGlu, MassSq, MassTop, MasssigmaO, MassphiO, mu};
(Format[Power[#, 2], CForm] := Format[ToString[#] <> "2", OutputForm])& /@ squareAbbreviation;
Format[Power[x_ /; MemberQ[squareAbbreviation, x], i_Integer /; EvenQ[i]], CForm] := Format["pow<" <> ToString[i/2] <> ">(" <> ToString@CForm[x^2] <> ")", OutputForm];
Protect[Power];

model="MRSSM";
tempName=Last@StringSplit[$ScriptCommandLine[[2]], "/"];
tempName=StringSplit[tempName, "_"];
fileName=tempName[[1]] <> "_" <> tempName[[2]];

(* write final expression to a file *)
For[i=-2, i<=0, i++,
   suffix=Switch[i, -2, "dp", -1, "sp", 0, "finite"];
   resOrd = SeriesCoefficient[resAllOrd, i];
WriteString[
   FileNameJoin[{"src", "models", model, ToLowerCase[model] <> "_" <> fileName <> "_soft_" <> suffix <> ".cpp"}],
"double " <> model <> "::matrixSoft_" <> fileName <> "_" <> suffix <> "(const double alphas, const double S, const double th" <> If[i==-2, "", ", const double dS, const double mu"] <> ") const {\n" <>
   If[!FreeQ[resOrd, beta], "   const double b = std::sqrt(1. - 4.*pow<2>(MassSq)/S);\n", ""] <>
   If[!FreeQ[resOrd, Log[dS]], "   const double lndS = std::log(dS);\n", ""] <>
   If[!FreeQ[resOrd, Cos[th]], "   const double cth = std::cos(th);\n" , ""] <>
   If[!FreeQ[resOrd, Sin[th]], "   const double sth = std::sin(th);\n" , ""] <>
   If[!FreeQ[resOrd, MassGlu^2], "   const double MassGlu2 = pow<2>(MassGlu);\n", ""] <>
   If[!FreeQ[resOrd, MassTop^2],  "   const double MassTop2 = pow<2>(MassTop);\n", ""] <>
   If[!FreeQ[resOrd, MasssigmaO^2],  "   const double MasssigmaO2 = pow<2>(MasssigmaO);\n", ""] <>
   If[!FreeQ[resOrd, mu2],  "   const double mu2 = pow<2>(mu);\n", ""] <>
   If[!FreeQ[resOrd, alphas3],  "   const double alphas3 = pow<3>(alphas);\n", ""] <>
   If[!FreeQ[resOrd, MassSq^2],  "   const double MassSq2 = pow<2>(MassSq);\n", ""] <>
"   const double res = " <>
   StringReplace[
      ToString@CForm[
            SeriesCoefficient[resAllOrd, i] /.
               Log[dS] -> lndS /.
               PolyLog[2, x_] :> Li2[x] /.
               Cos[th] -> cth /. Sin[th] -> sth /. Sin -> sin /.
               beta->b /. Cos->cos //. replace
      ],
      {"Li2" -> "polylogarithm::Li2", "EulerGamma" -> "euler"}
   ] <> ";\n" <>
"   return pow<3>(alphas*pi)*res;\n}"
]
];
