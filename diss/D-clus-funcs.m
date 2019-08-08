ik[key_String]:=key/*Import
ikk[key_String]:=<|
	"evs"->"ndeevs"/*Import/*(AssociationThread[Table[tssf["i`1`",i],{i,Length[#]}],h2mev[#]]&),
	"etr"->"ndeevs"/*Import/*h2mev/*(AssociationThread[Table[tssf["j`1`",j],{j,Length[#]}],#-#[[1]]]&),
	"mx"->"ndeparams"/*Import/*"mx",
	"my"->"ndeparams"/*Import/*"my",
	"m0"-><|"mx"->"ndeparams"/*Import/*"mx","my"->"ndeparams"/*Import/*"my"|>/*(#mx*#my/((Sqrt[#mx]+Sqrt[#my])^2)&),
	"mxy"-><|"mx"->"ndeparams"/*Import/*"mx","my"->"ndeparams"/*Import/*"my"|>/*(#mx/#my&),
	"fox"->"fox"/*Import/*(AssociationThread[Table[tssf["i`1`",i],{i,12}],Map[AssociationThread[Table[tssf["j`1`",j],{j,13-Length[#],12}],#]&,#]]&),
	"foy"->"foy"/*Import/*(AssociationThread[Table[tssf["i`1`",i],{i,12}],Map[AssociationThread[Table[tssf["j`1`",j],{j,13-Length[#],12}],#]&,#]]&)
|>[key]
h2mev[quant_]:=UnitConvert[Quantity[quant,"Hartrees"],"Millielectronvolts"]

(* Of all the built-in mathematica functions, not having a light/dark function is absurd to me *)
lida[col_,ld_]:=If[ld>=0.5,Darker[col,2*(ld-0.5)],Lighter[col,2*(0.5-ld)]]

(* Easier memory monitoring *)
bmb[bytes_]:=N@UnitConvert[Quantity[bytes,"Bytes"],"Megabytes"]
bgb[bytes_]:=N@UnitConvert[Quantity[bytes,"Bytes"],"Gigabytes"]
mumb:=bmb[MemoryInUse[]]
mugb:=bmb[MemoryInUse[]]
