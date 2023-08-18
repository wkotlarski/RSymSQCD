SetOptions[Paint, PaintLevel -> {Classes}, ColumnsXRows -> {4, 5}]
$PaintSE = MkDir["diagrams"];
DoPaint[diags_, file_, opt___] := Paint[diags, opt,
  DisplayFunction -> (Export[ToFileName[$PaintSE, file <> ".ps"], #]&)];
