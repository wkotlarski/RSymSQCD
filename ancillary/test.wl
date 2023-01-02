<<FeynArts`
<<FormCalc`

(* topologies *)
top = CreateTopologies[0, 2 -> 3];

(* insertions *)
ins = InsertFields[
   top,
   {V[5], F[3, {1}]} -> {S[13, {1}], -S[13, {1}], F[3, {1}]},
   Model->"MRSSMEWSB_myRenConst"
];
(* diagram removal *)
ins = DiagramSelect[ins, (Not[SChannelQ[F[15]][##]]&& Not[SChannelQ[-F[15]][##]])&];

(* paint diagrams *)
$PaintSE = MkDir["our.diagrams"];
DoPaint[diags_, file_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, file <> ".ps"], #]&), FieldNumbers->True]
SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]
DoPaint[top, "topos"];
DoPaint[ins, "ins"];

Neglect[MassFu[1]]=0;

amp = CreateFeynAmp[ins];
result = CalcFeynAmp[amp, FermionChains->Chiral, MomElim->5, DotExpand->True];

result = PolarizationSum[result] //. Abbr[] //. Subexpr[] /. MassFu[1] ->0 /. Den[x_, y_] :> 1/(x-y);

result = (
  result /. eta[1]->k[2]
  /. (Pair[k[i_Integer], k[j_Integer]] /; j<i) :> Pair[k[j], k[i]]
  /. Pair[k[2], k[2]]|Pair[k[1], k[1]]|Pair[k[5], k[5]] -> 0
  /. Pair[k[3],k[3]]|Pair[k[4],k[4]] :> MassSq^2
  /. MassSu[1] -> MassSq
  /. Pair[k[1], k[2]] -> S/2
  /. Pair[k[1], k[3]] -> (MassSq^2 - T)/2
  /. Pair[k[1], k[4]] -> (MassSq^2 - T14)/2
  /. Pair[k[2], k[3]] -> (MassSq^2 - U)/2
  /. Pair[k[2], k[4]] -> (MassSq^2 - T24)/2
  /. Eps[k[i1_], k[i2_], k[i3_], k[i4_]] :> Signature[{i1, i2, i3, i4}] * Eps@@(k /@ Sort[{i1, i2, i3, i4}]));

If[FullSimplify@Coefficient[result, Eps[k[1], k[2], k[3], k[4]]] =!= 0,
   Print["Error: Eps tensor doesn't wanish!"];
   Quit[1]
];

result = FullSimplify[result /. Eps[k[1], k[2], k[3], k[4]] -> 0];

Put[CForm[result], "gu_suLsuLdaggeru_expr.cpp"];

