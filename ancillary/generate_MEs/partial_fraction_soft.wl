Get[FileNameJoin[{DirectoryName @ $InputFileName, "Apart.m"}]];

origExpr = Expand[Get[$ScriptCommandLine[[2]]] /. S35 -> S35p + MassSq^2];

findPower[expr_]:=Module[{res, x},
   res = expr /. -S35p - T15 - T25 -> -x /. S35p + T15 + T25 -> x /. S35p -> x /. T15 -> x /. T25 ->x;
   If[FreeQ[res, x], 0, (List@@Series[res, {x,0,0}])[[-3]]]
];

expr = Plus@@(If[findPower[#]<-2, ApartAll[#, {S35p,T25,T15}] /. ApartIR[pcs_,cs_,np_,vars_] :> Times@@((((#1.vars&)/@pcs+cs))^np), #]& /@ List@@origExpr);
Print["Leftover bad terms...."];
expr = Plus@@(If[findPower[#]<-2, Print[#];Apart[#], #]& /@ expr);

expr = Simplify@Select[expr, findPower[#]<-2&] + Select[expr, findPower[#]>=-2&];

minPow = Min[findPower /@ List@@Expand[expr]];
Print["INFO: minimal power ", minPow];
If[minPow < -2, Quit[1]];

If[!Simplify[origExpr == expr],
   Print["Expressions differ"];
   Quit[1];
];

fileName = Last@StringSplit[$ScriptCommandLine[[2]], "/"];
fileName = StringSplit[fileName, "."];
Put[Expand[expr], FileNameJoin[{DirectoryName @ $InputFileName, "soft_mes", fileName[[1]] <> "-fractioned.m"}]];
