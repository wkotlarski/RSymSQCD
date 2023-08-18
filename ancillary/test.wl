Needs["MultivariateApart`"];

{expr, dens, qis} = AbbreviateDenominators[(2y-x)/(y(x+y)(y-x))];
(* {-(q1 q2 q3 (-x + 2 y)), {x - y, y, x + y}, {q1, q2, q3}} *)


ord = ApartOrder[dens, qis];

gb = ApartBasis[dens, qis, ord];

Print[ApartReduce[expr, gb, ord]];

