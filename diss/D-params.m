(* Physical Constants *)
ce0 = Quantity[ "ElectricConstant" ];
ckk = 1/( 4 * Pi * ce0 );
cee = Quantity[ "ElementaryCharge" ];
cm0 = Quantity[ "ElectronMass" ];
chh = Quantity[ "PlanckConstant" ];
chb = Quantity[ "ReducedPlanckConstant" ];
ccc = Quantity[ "SpeedOfLight" ];
ckb = Quantity[ "BoltzmannConstant" ];

(* Units *)
qbhr = Quantity[ "BohrRadius" ];
qang = Quantity[ "Angstroms" ];
qnm = Quantity[ "Nanometers" ];
qmet = Quantity[ "Meters" ];
qsec = Quantity[ "Seconds" ];
qmev = Quantity[ "Millielectronvolts" ];
qev = Quantity[ "Electronvolts" ];
qvol = Quantity[ "Volts" ];
qhart = Quantity[ "Hartrees" ];
qconc = Quantity[ "Meters" ^ ( -2 ) ];
qfreq = Quantity[ "Seconds" ^ ( -1 ) ];

(* Common parameters *)
lBN = UnitConvert[ Quantity[ 0.333, "Nanometers" ], "BohrRadius" ];
n0 = ((5 * 10^(15)) * qconc);
g0bn = ((10^13) * qfreq);
g0fs = ((10^14) * qfreq);

macheps = $MachineEpsilon;
pio2 = Pi / 2;
finpi = pio2 - macheps;

(* Default function options *)
(* ndeopts... etc. *)

(* Raw Material Parameter Dataset/*{{{*/*)
rawmpds = Dataset[<|
	"t"-><|
		"mos2"	->	<|
				"hi"	->	<| "mu"->	0.28 * cm0,	"chi"->	6.600 * qang,	"l"->	6.180 * qang,	"eps"->	10|>,
				"lo"	->	<| "mu"->	0.16 * cm0,	"chi"->	7.112 * qang,	"l"->	6.180 * qang,	"eps"->	10|>
		|>,
		"mose2"	->	<|
				"hi"	->	<| "mu"->	0.31 * cm0,	"chi"->	8.230 * qang,	"l"->	6.527 * qang,	"eps"->	10|>,
				"lo"	->	<| "mu"->	0.27 * cm0,	"chi"->	8.461 * qang,	"l"->	6.527 * qang,	"eps"->	10|>
		|>,
		"ws2"	->	<|
				"hi"	->	<| "mu"->	0.23 * cm0,	"chi"->	6.030 * qang,	"l"->	6.219 * qang,	"eps"->	10|>,
				"lo"	->	<| "mu"->	0.15 * cm0,	"chi"->	6.393 * qang,	"l"->	6.219 * qang,	"eps"->	10|>
		|>,
		"wse2"	->	<|
				"hi"	->	<| "mu"->	0.27 * cm0,	"chi"->	7.180 * qang,	"l"->	6.575 * qang,	"eps"->	10|>,
				"lo"	->	<| "mu"->	0.15 * cm0,	"chi"->	7.571 * qang,	"l"->	6.575 * qang,	"eps"->	10|>
		|>
	|>,
	"x"-><|
		"fssi"	->	<| "eg" -> 1.9 * qmev	, "d0" -> 0.460 * qang, "vf" -> Quantity[6.5*10^5, "Meters"/"Seconds"]	, "l" -> 4.0 * qang , "eps" -> 11.9 |>,
		"fsge"	->	<| "eg" -> 33 * qmev	, "d0" -> 0.676 * qang, "vf" -> Quantity[6.2*10^5, "Meters"/"Seconds"]	, "l" -> 4.5 * qang , "eps" -> 16.0 |>,
		"fssn"	->	<| "eg" -> 101 * qmev	, "d0" -> 0.850 * qang, "vf" -> Quantity[5.5*10^5, "Meters"/"Seconds"]	, "l" -> 5.0 * qang , "eps" -> 24.0 |>,
		"shsi"	->	<| "eg" -> 13.5 * qmev, "d0" -> 0.460 * qang, "vf" -> Quantity[4.33*10^5, "Meters"/"Seconds"]	, "l" -> 3.33 * qang, "eps" -> 11.9 |>,
		"bhsi"	->	<| "eg" -> 19 * qmev	, "d0" -> 0.460 * qang, "vf" -> Quantity[5.06*10^5, "Meters"/"Seconds"]	, "l" -> 3.33 * qang, "eps" -> 11.9 |>
	|>,
	"p"-><|
		"mla"		->	<| "mex" -> 0.16 * cm0	, "mey" -> 1.24 * cm0		, "mhx" -> 0.15 * cm0		, "mhy" -> 4.92 * cm0		, "cite" -> "Peng2014" 	|>,
		"mlb"		->	<| "mex" -> 0.10 * cm0	, "mey" -> 1.3 * cm0		, "mhx" -> 0.20 * cm0		, "mhy" -> 2.8 * cm0		, "cite" -> "Tran2014a"	|>,
		"mlc"		->	<| "mex" -> 0.1990 * cm0, "mey" -> 0.7527 * cm0	, "mhx" -> 0.1678 * cm0	, "mhy" -> 5.3525 * cm0	, "cite" -> "Paez2016" 	|>,
		"mld"		->	<| "mex" -> 0.17 * cm0	, "mey" -> 1.12 * cm0		, "mhx" -> 0.15 * cm0		, "mhy" -> 6.35 * cm0		,	"cite" -> "Qiao2014" 	|>,
		"bld"		->	<| "mex" -> 0.18 * cm0	, "mey" -> 1.13 * cm0		, "mhx" -> 0.15 * cm0		, "mhy" -> 1.81 * cm0		, "cite" -> "Qiao2014" 	|>,
		"chi"		-> ( 4.1 * qang ),
		"l"			-> ( 5.41 * qang ),
		"eps"		-> 11
	|>
|>];
(* }}}*)

genNDE = NDEigensystem[
	{
		transdiffeq,
		dbc
	},
	func,
	reg,
	nmax,
	opts
]

(* Association with the diff eq's for each material *)
gendiffeq = <|
	"t" -> - ( f''[r] * ( 1 / ( 2 * mu ))) - ( f'[r] * ( 1 / ( 2 * r * mu ))) + ( f[ r ] * ( genrk["t"] + ( angm^2 / ( 2 * r^2 * mu )) + shift)),
	"x" -> - ( f''[r] * ( 1 / ( 2 * mu ))) - ( f'[r] * ( 1 / ( 2 * r * mu ))) + ( f[ r ] * ( genrk["t"] + ( angm^2 / ( 2 * r^2 * mu )) + shift)),
	"p" -> - ( 1 / 2 ) * ( ( 1 / mx ) * D[ f[ x, y ], { x, 2 }, { y, 0 } ] + ( 1 / my ) * D[ f[ x, y ], { x, 0 }, { y, 2 } ] ) + ( potfun[ nbn, kappa, matkeys ] + shift ) * f[ x, y ]
|>

(* Association of potentials *)
genrk = <|
	"t" -> - ( Pi / ( 2 * kappa * rho0 )) * ( StruveH[ 0, Sqrt[ r^2 + DD^2 ] / rho0 ] - BesselY[ 0, Sqrt[ r^2 + DD^2 ] / rho0 ] ),
	"x" -> - ( Pi / ( 2 * kappa * rho0 )) * ( StruveH[ 0, Sqrt[ r^2 + DD^2 ] / rho0 ] - BesselY[ 0, Sqrt[ r^2 + DD^2 ] / rho0 ] ),
	"p" -> - ( Pi / ( 2 * kappa * rho0 )) * ( StruveH[ 0, Sqrt[ x^2 + y^2 + DD^2 ] / rho0 ] - BesselY[ 0, Sqrt[ x^2 + y^2 + DD^2 ] / rho0 ] )
|>

(* Association of diff eq transformation rules *)
gentransrules = <|
	"t" -> (Simplify[ # /. f -> (g[ArcTan[#]] &) /. r-> (Tan[s]), 0 < s < Pi/2]&),
	"x" -> (Simplify[ # /. f -> (g[ArcTan[#]] &) /. r-> (Tan[s]), 0 < s < Pi/2]&),
	"p" -> (Simplify[ # /.{f -> (g[ArcTan[#1/c], ArcTan[#2/c]] &) } /. { x -> ( c * Tan[\[Xi]]), y -> ( c * Tan[\[Psi]]) }, { -(Pi/2) < \[Xi] < (Pi/2), -(Pi/2) < \[Psi] < (Pi/2) } ]&)
|>

genndeotherstuff = <| (* [[1]]: DirichletCondition || [[2]]: 2nd arg to NDE || [[3]]: Region (3rd arg to NDE) *)
	"t" -> { DirichletCondition[g[s]==0, s==Pi/2], g[s], {s, 0, Pi/2} },
	"x" -> { DirichletCondition[g[s]==0, s==Pi/2], g[s], {s, 0, Pi/2} },
	"p" -> { DirichletCondition[u[\[Xi],\[Psi]] == 0, Abs[\[Xi]] == finpi || Abs[\[Psi]] == finpi], u[\[Xi],\[Psi]], {\[Xi], \[Psi]} \[Element] Rectangle[{-finpi,-finpi},{finpi,finpi}] }
|>

anisomu[me_,mh_]:=(me*mh/(me+mh))

(* Functions to prep the raw MPDS for NDEigensystem *)
xprep[ ep_Quantity ] := (<|
	"ega" -> UnitConvert[ Abs[ #["eg"] + ( #["d0"] * cee * ep ) ], "Millielectronvolts" ],
	"mua" -> UnitConvert[ Abs[ #["eg"] + ( #["d0"] * cee * ep ) ] / ( 2 * ( #[ "vf" ]^2 ) ), "ElectronMass" ],
	"egb" -> UnitConvert[ Abs[ #["eg"] - ( #["d0"] * cee * ep ) ], "Millielectronvolts" ],
	"mub" -> UnitConvert[ Abs[ #["eg"] - ( #["d0"] * cee * ep ) ] / ( 2 * ( #[ "vf" ]^2 ) ), "ElectronMass" ],
	"rho0" ->Function[ { kappa }, #["l"] * #["eps"] / kappa ]
|>&)
xprep[ ep_Real | ep_Integer ] := (<|
	"ega" -> UnitConvert[ Abs[ #["eg"] + ( #["d0"] * cee * Quantity[ ep, "Volts"/"Angstroms" ] ) ], "Millielectronvolts" ],
	"mua" -> UnitConvert[ Abs[ #["eg"] + ( #["d0"] * cee * Quantity[ ep, "Volts"/"Angstroms" ] ) ] / ( 2 * ( #[ "vf" ]^2 ) ), "ElectronMass" ],
	"egb" -> UnitConvert[ Abs[ #["eg"] - ( #["d0"] * cee * Quantity[ ep, "Volts"/"Angstroms" ] ) ], "Millielectronvolts" ],
	"mub" -> UnitConvert[ Abs[ #["eg"] - ( #["d0"] * cee * Quantity[ ep, "Volts"/"Angstroms" ] ) ] / ( 2 * ( #[ "vf" ]^2 ) ), "ElectronMass" ],
	"rho0" ->Function[ { kappa }, #["l"] * #["eps"] / kappa ]
|>&)

pprep[ muind_String ] := <|
	"mux" -> muind /* (anisomu[ #mex, #mhx ]&),
	"muy" -> muind /* (anisomu[ #mey, #mhy ]&),
	"chi" -> "chi",
	"l"		-> "l",
	"eps" -> "eps"
|>

(* Other calculation parameters *)
envs = <|
	"fs" -> 1,
	"su" -> ( 1 + 3.9 ) / 2,
	"hu" -> ( 1 + 4.89 ) / 2,
	"sh" -> ( 3.9 + 4.89 ) / 2,
	"hh" -> 4.89
|>;

(* xgetmu[ ep_ ] := *)

nbntod[ n_Integer ] := ( n * lBN )
nbntod[ n_Integer, matkeys__ ] := If[ n == -1, Quantity[ 0, "BohrRadius" ], ( n * lBN ) + matpars[ matkeys, "l" ] ]

(* Define potential functions that can pull from the matpar dataset *)
vkeld[ nbn_Integer, kappa_, matkeys__] := With[
	{
		rho0 = If[
			( List[ matkeys ] )[[ 1 ]] == "x",
			QuantityMagnitude[ matpars[ matkeys, "l" ] * matpars[ matkeys, "eps" ] / ( 2 * kappa ), "BohrRadius" ],
			QuantityMagnitude[ ( 2 * Pi * matpars[ matkeys, "chi" ] ) / kappa, "BohrRadius" ]
		],
		DD = QuantityMagnitude[ nbntod[ nbn, matkeys ], "BohrRadius" ]
	},
	If[
		( List[ matkeys ] )[[ 1 ]] == "p",
		- ( Pi / ( 2 * kappa * rho0 )) * ( StruveH[ 0, Sqrt[ x^2 + y^2 + DD^2 ] / rho0 ] - BesselY[ 0, Sqrt[ x^2 + y^2 + DD^2 ] / rho0 ] ),
		- ( Pi / ( 2 * kappa * rho0 )) * ( StruveH[ 0, Sqrt[ r^2 + DD^2 ] / rho0 ] - BesselY[ 0, Sqrt[ r^2 + DD^2 ] / rho0 ] )
	]
]

(* Set up the differential equation (schrodinger) *)
makeDE[ angm_, nbn_Integer, kappa_, potfun_, matkeys__ ] := With[
	{
		mu = QuantityMagnitude[ matpars[ matkeys, "mu" ] ],
		shift = 1
	},
	If[
		( List[ matkeys ] )[[ 1 ]] == "p",
		- ( 1 / 2 ) * ( ( 1 / mx ) * D[ f[ x, y ], { x, 2 }, { y, 0 } ] + ( 1 / my ) * D[ f[ x, y ], { x, 0 }, { y, 2 } ] ) + ( potfun[ nbn, kappa, matkeys ] + shift ) * f[ x, y ],
		- ( f''[r] * ( 1 / ( 2 * mu ))) - ( f'[r] * ( 1 / ( 2 * r * mu ))) + ( f[ r ] * ( potfun[ nbn, kappa, matkeys ] + ( angm^2 / ( 2 * r^2 * mu )) + shift))
	]
]

(* NDEigensystem Wrapper *)
(* solveNDE[...] := Module[*)
(*   {*)

(*   },*)
(*   makeDE[ angm, nbn, kappa, potfun, matkeys ]*)

(* Use an Association to write NDE functions??? *)
(* This assoc will produce one solution for one input (no tables or anything yet) *)
(* solveNDE = Association[*)
(*   "solve" -> Module[*)
(*     {},*)

(*   ]*)
(* ];*)

(* Write a function that will run **inside** the MPDS that will take parameters appropriately and return a result for NDEigensystem *)
(* dsnde[externalparams_]:=Module[ [> External params are e.g. kappa, D, nmax, mmax (for TMDCs/Xenes) etc <]*)
(*   {*)
(*   },*)
(*   diffeq = gendiffeq[...];*)
(*   transdiffeq = Simplify[ diffeq /. mattransrules ];*)
(*   {{ev, ef},stats} = med[Unevaluated[NDEigensystem[*)
(*     {*)
(*       transdiffeq,*)
(*       DirichletCondition[ [> 2D vs. 1D rules here <]]*)
(*     },*)
(*     u[ [> coords from diff eq...? <] ],*)
(*     [> coords from diff eq somehow <],*)
(*     nmax,*)
(*     ndeopts*)
(*   ]]];*)
(*   nef [> 2D/1D function for normalizing EFs <];*)
(*   {*)
(*     {params},*)
(*     {evs},*)
(*     {nefs}*)
(*   }*)

