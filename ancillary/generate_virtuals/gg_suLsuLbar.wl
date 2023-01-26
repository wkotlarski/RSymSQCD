<<FeynArts.m;
<<FormCalc.m;

SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]
$PaintSE = MkDir["diagrams"];
DoPaint[diags_, file_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, file <> ".ps"], #]&)];

(* topologies *)
treeTop = CreateTopologies[0, 2->2];
virtTop = CreateTopologies[1, 2->2, ExcludeTopologies->{Tadpoles,WFCorrections}];
ctTop = CreateCTTopologies[1, 2->2, ExcludeTopologies->{TadpoleCTs,WFCorrectionCTs}];

(* process *)
SetOptions[InsertFields, Model->"MRSSMEWSB_myRenConst", InsertionLevel->Particles, ExcludeParticles->{V[2], V[3]}];
process = {V[5], V[5]} -> {S[13,{1}], -S[13,{1}]};
Neglect[MassFu[1]] = Neglect[Mass[Fu[2]]] = Neglect[MassFd[_]] = Neglect[g1] = Neglect[g2] = Neglect[Yu[__]] = Neglect[Yd[__]] = 0;

(* insertions *)
bornIns = InsertFields[treeTop, process];
DoPaint[bornIns, "born"];
ctIns = InsertFields[ctTop, process];
DoPaint[ctIns, "ct"];
virtIns = InsertFields[virtTop, process];
DoPaint[virtIns, "virt"];

(* amplitudes *)
SetOptions[CalcFeynAmp, FermionChains->VA, Dimension->0, FermionOrder->None];
SetOptions[PolarizationSum, Dimension->0, GaugeTerms->False];
born = CalcFeynAmp[CreateFeynAmp[bornIns]];
virt = CalcFeynAmp[CreateFeynAmp[virtIns]];
counter = CalcFeynAmp[CreateFeynAmp[ctIns]];

_Hel=0;
colour=ColourME[born];
born2=SquaredME[born];
born2=Simplify[born2[[1]]/.born2[[2]] /. colour /. Den[x_,y_]:>1/(x-y)];
born2=FullSimplify[PolarizationSum[born2] //. SubExpr[] //. Abbr[] /. U -> 2*MassSu[1]^2 - S -T /. Pair[eta[n_], k[4]] :> Pair[eta[n], k[1]] + Pair[eta[n], k[2]] - Pair[eta[n], k[3]]];

colour=ColourME[born,virt];
virt2=SquaredME[born,virt];
virt2=PolarizationSum[virt2[[1]]/.virt2[[2]]] /. colour //. SubExpr[] //. Abbr[] /. MassFu[1|2]|MassFd[_]|Yu[_,_]|Yd[_,_]|g1|g2 -> 0 /. U -> 2*MassSu[1]^2 - S -T/. Den[x_,y_]:>1/(x-y)/. Pair[eta[n_], k[4]] :> Pair[eta[n], k[1]] + Pair[eta[n], k[2]] - Pair[eta[n], k[3]];
virt2=virt2;

colour=ColourME[born, counter];
CT2=SquaredME[born, counter];
CT2=PolarizationSum[CT2[[1]] /.CT2[[2]]]  /. colour //. SubExpr[] //. Abbr[] /. MassFu[1|2]|MassFd[_]|Yu[_,_]|Yd[_,_]|g1|g2 ->0 /. MassFd[_]->0 /. U -> 2*MassSu[1]^2 - S -T/. Den[x_,y_]:>1/(x-y)/. Pair[eta[n_], k[4]] :> Pair[eta[n], k[1]] + Pair[eta[n], k[2]] - Pair[eta[n], k[3]];
CT2=CT2;

renConstSubs = (# -> RenConst[#])& /@ {
   dZSu1[1], dZGlL1, dZgs1, dMGl1, dZGG1, dMSu1[1] 
};
renConstSubs = renConstSubs/. MassFu[1|2]|MassFd[_]|Yu[_,_]|Yd[_,_]|g1|g2 ->0

CT2 = CT2/.renConstSubs;

total = 2*(virt2 + CT2);

uvDiv = FullSimplify@UVDivergentPart[total/. Dminus4->0];
If[uvDiv =!= 0,
  Print["UV divergence: ", uvDiv];
  Quit[1];
];

Install["/home/kotlarskiw/HEP-software/gcc/LoopTools-2.16/bin/LoopTools"];

numerics = {
  g3 -> 0.1184, 
  MassSu[_] -> 1500, MassSd[_] -> 1500, MassSq -> 1500,
  MassGlu -> 1000, MDO -> 1000,
  S -> 6000^2, T -> -22208172, 
  mu->1500, 
  MassFu[3]->172, 
  MasssigmaO->5000, 
  MassphiO->Sqrt[5000^2 + 4*1000^2]
};

normalization = born2*g3^2/(8*Pi^2) /. numerics /. FiniteGs|Dminus4->0;

SetMudim[1500^2];
SetLambda[-2];
dpole = Simplify[total/normalization /. numerics /. FiniteGs|Dminus4->0];
Print["double pole: ", dpole];
dtospole = Coefficient[total, Dminus4]*(-2)/normalization /. numerics /. FiniteGs->0;

SetLambda[-1];
spole = Re[Simplify[total/normalization /. numerics /. FiniteGs|Dminus4->0] + dtospole];
Print["single pole: ", spole];
stofinite = Simplify[Coefficient[total, Dminus4]/normalization /. numerics /. FiniteGs|Dminus4 -> 0]*(-2);

SetLambda[0];
finite = Re[Simplify[total/normalization /. numerics /. FiniteGs -> 1 /. Dminus4->0] + dpole*Pi^2/6*0 + stofinite];
Print["finite: ", finite];
