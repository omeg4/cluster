SetDirectory["/home/mbrunetti/cluster/diss/2dp/v1"]

(* Physical Constants *)
ce0 = Quantity[ "ElectricConstant" ];
ckk = 1/( 4 * Pi * e0 );
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

tssf[x___]:=ToString@StringForm[x]
med[f_]:=({#["Result"],KeyTake[#,{"AbsoluteTiming","MessagesText","Timing"}]}&[KeyTake[EvaluationData[f],{"Result","AbsoluteTiming","MessagesText","Timing"}]])
med[f_,x___]:=({#["Result"],KeyTake[#,{"AbsoluteTiming","MessagesText","Timing"}]}&[KeyTake[EvaluationData[f[x]],{"Result","AbsoluteTiming","MessagesText","Timing"}]])
bg[expr_, clr_RGBColor] := Framed[expr, FrameStyle -> None, Background -> clr]
bg[clr_RGBColor, expr_] := Framed[expr, FrameStyle -> None, Background -> clr]
fulldate:=DateString[{"Month","/","Day"," @ ","Time"," | "}]
time:=DateString[{"Time"}]

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

(* /*{{{*/ Assoc. of KE operators + (shift * f) for an/iso exciton diff. eq. *)
diffeqkeop = <|
	"iso" -> - ( f''[r] * ( 1 / ( 2 * mu ))) - ( f'[r] * ( 1 / ( 2 * r * mu ))) + ( f[ r ] * ( angm^2 / ( 2 * r^2 * mu ))) + ( shift * f[ r ] ),
	"ani" -> - ( 1 / 2 ) * ( ( 1 / mx ) * D[ f[ x, y ] , { x, 2 }, { y, 0 } ] + ( 1 / my ) * D[ f[ x, y ], { x, 0 }, { y, 2 } ] ) + ( shift * f[ x, y ] )
|>
(* /*}}}*/*)

(* /*{{{*/ Assoc. of transformation rules for the diff eqs. *)
diffeqtransrules = <|
	"1dp" -> (Simplify[ # /. f -> (g[ ArcTan[ # / c ] ] &) /. r-> ( c * Tan[s] ), 0 < s < Pi/2]&),
	"2dc" -> (Simplify[ # /.{f -> (g[ ArcTan[ #1 / c ], ArcTan[ #2 / c ] ] &) } /. { x -> ( c * Tan[ \[Xi] ]), y -> ( c * Tan[ \[Psi] ]) }, { -( Pi/2 ) < \[Xi] < ( Pi/2 ), -( Pi/2 ) < \[Psi] < ( Pi/2 ) } ]&),
	"2dp" -> (Simplify[ Simplify[
						# /.
				{ f -> ( u[ Sqrt[ #1^2 + #2^2 ], ArcTan[ #2 / #1 ] ] &) } /.
				{ x -> r * Cos[ \[Theta] ], y -> r * Sin[ \[Theta] ] }, (* These two rules transform 2D cartesian to 2D polar *)
				{ -pio2 < \[Theta] < pio2, 0 <= r < Infinity } (* Apparently simplifying here with -pio2<th<pio2 is necessary because otherwise the expression doesn't simplify properly and NDE breaks... *)
			] /.
		{ u -> ( g[ ArcTan[ #1 / c ], #2 ] &) } /. { r -> ( c * Tan[ s ] ) },
		{ 0 <= s < Pi/2 }
	]&)
|>
(* /*}}}*/*)

(* /*{{{*/ Assoc. of BC for diff eqs. *)
diffeqbc = <|
	"1dp" -> DirichletCondition[ g[ s ] == 0, s == Pi/2 ],
	"2dc" -> DirichletCondition[ g[ \[Xi], \[Psi] ] == 0, Abs[ \[Xi] ] == finpi || Abs[ \[Psi] ] == finpi ],
	"2dp" -> {
		DirichletCondition[ g[ \[Psi], \[Theta] ] == 0, \[Psi] == finpi && 0 < \[Theta] <= 2*Pi ],
		PeriodicBoundaryCondition[g[\[Psi],\[Theta]], \[Theta]==0, TranslationTransform[{0,2*Pi}]]
	}
|>
(* /*}}}*/*)

(* /*{{{*/ Assoc. of NDE-ready differential equations *)
fulldiffeqs = <|
	"1dp" -> <|
		"rk" -> diffeqtransrules[ "1dp" ][ diffeqkeop[ "iso" ] + pots[ "iso", "rk" ] * f[ r ] ],
		"cl" -> diffeqtransrules[ "1dp" ][ diffeqkeop[ "iso" ] + pots[ "iso", "cl" ] * f[ r ] ]
	|>,
	"2dp" -> <|
		"rk" -> diffeqtransrules[ "2dp" ][ diffeqkeop[ "ani" ] + pots[ "ani", "rk" ] * f[ x, y ] ],
		"cl" -> diffeqtransrules[ "2dp" ][ diffeqkeop[ "ani" ] + pots[ "ani", "cl" ] * f[ x, y ] ]
	|>,
	"2dc" -> <|
		"rk" -> diffeqtransrules[ "2dc" ][ diffeqkeop[ "ani" ] + pots[ "ani", "rk" ] * f[ x, y ] ],
		"cl" -> diffeqtransrules[ "2dc" ][ diffeqkeop[ "ani" ] + pots[ "ani", "cl" ] * f[ x, y ] ]
	|>
|>
(* /*}}}*/*)

(* /*{{{*/ Assoc. of an/iso RK/C potentials *)
pots = <|
	"iso" -> <|
		"rk" -> - ( Pi / ( 2 * kappa * rho0 )) * ( StruveH[ 0, Sqrt[ r^2 + z^2 ] / rho0 ] - BesselY[ 0, Sqrt[ r^2 + z^2 ] / rho0 ] ),
		"cl" -> - ( 1 / ( kappa * Sqrt[ r^2 + z^2 ] ) )
	|>,
	"ani" -> <|
		"rk" -> - ( Pi / ( 2 * kappa * rho0 )) * ( StruveH[ 0, Sqrt[ x^2 + y^2 + z^2 ] / rho0 ] - BesselY[ 0, Sqrt[ x^2 + y^2 + z^2 ] / rho0 ] ),
		"cl" -> - ( 1 / ( kappa * Sqrt[ x^2 + y^2 + z^2 ] ) )
	|>
|>
(* /*}}}*/ *)

(* /*{{{*/ Repurpose 2D Polar to fit new approach *)
d2pndetest[ params_List, calcparams_List ]:=Module[
	{
		mx = params[[1]],
		my = params[[2]],
		chi = params[[3]],
		l = params[[4]],
		kappa = params[[5]],
		z = params[[6]],
		rho0 = ( 2 * Pi * params[[3]] / params[[5]] ),
		tde = ( diffeqtransrules[ "2dp" ][ diffeqkeop[ "ani" ] + pots[ "ani", "rk" ]*f[ x, y ] ] ),
		tdenum = With[ { mx = params[[1]], my = params[[2]], chi = params[[3]], l = params[[4]], kappa = params[[5]], rho0 = ( 2 * Pi * params[[3]] / params[[5]] ), z = params[[6]], c = calcparams[[2]], shift = 1 }, Evaluate[ ( diffeqtransrules[ "2dp" ][ diffeqkeop[ "ani" ] + pots[ "ani", "rk" ]*f[ x, y ] ] ) ] ],
		eps = calcparams[[1]], maxi = 10^12,
		c = calcparams[[2]],
		shift = 1,
		nmax = 12
	},
	{evs,efs}=NDEigensystem[
		{
			tdenum,
			DirichletCondition[g[ s, \[Theta] ] == 0, ( s == (pio2) ) && 0 < \[Theta] <= 2*Pi],
			PeriodicBoundaryCondition[ g[ s, \[Theta] ], \[Theta]==0, TranslationTransform[{0,2*Pi}]]
		},
	 g[ s, \[Theta] ],
	 { s, \[Theta] } \[Element] Rectangle[ { 0, 0 }, { pio2, 2*Pi } ],
	 nmax,
	 Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->eps}}},"Eigensystem"->{"Arnoldi","MaxIterations"->maxi},"VectorNormalization"->None}
	];
	{
		{ mx, my, chi, l, kappa, rho0 }, { tde, tdenum }, { eps, c },
		UnitConvert[Quantity[ evs-shift, "Hartrees" ], "Millielectronvolts" ]
		(* Head/@efs*) (* NOTE: NOT SAVING EFS FOR NOW, JUST WANT TO SEE WHAT HAPPENS TO EIGENENERGIES *)
	}
]
(* /*}}}*/*)

Table[
	Module[
		{ res,stats },
		{ res,stats } = med[Unevaluated[d2pndetest[ { 0.06, 0.9, 13, 9, 1, 0 }, { eps, cscale } ]]];
		{ Export[ "2dpol_eps-`1`_cscale-`2`_res.m", res ], Export[ "2dpol_eps-`1`_cscale-`2`_stats.m", stats ] }
	],
	{
		eps,
		{ 10^-3, 9*10^-4, 8*10^-4, 7*10^-4, 6*10^-4, 5*10^-4, 4*10^-4, 3*10^-4, 2*10^-4, 10^-4 }
	},
	{
	cscale,
		{ 1, 5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25, 30 }
	}
];

Quit[]
