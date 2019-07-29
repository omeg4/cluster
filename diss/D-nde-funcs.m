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

(* /*{{{*/ Biiiiiiig function that does everything I need in 1D up to optical props. (no polaritons yet) *)
d1pnde[ params_List ] := Module[
	{
		mu = params[[1]],
		chi = params[[2]],
		l = params[[3]],
		kappa = params[[4]],
		rho0 = (2*Pi*params[[2]]/params[[4]]),
		z = 0,
		c = 10,
		shift = 1,
		tde = (diffeqtransrules["1dp"][diffeqkeop["iso"] + pots["iso", "rk"]*f[r]]),
		tdenum = With[{mu = params[[1]], chi = params[[2]], l = params[[3]], kappa = params[[4]], rho0 = (2*Pi*params[[2]]/params[[4]]), z = 0, c = 10, shift = 1}, Evaluate@(diffeqtransrules["1dp"][diffeqkeop["iso"] + pots["iso", "rk"]*f[r]])],
		nmax = 3,
		res, evs, efs, esform, evform, efform, tefs, ndestats,
		memi, memf, memx,
		out, norms, normstats, ntefs, normcheck, normcheckstats, rmax = 10^6,
		rsq, xbhrad, f0, alpha, afac, alcon,
		lc, ldbr, leff, Mcav, Mex, hxsq, hcsq, gamcav, gamcavreff, gampol,
		exk, eck, ephotpct, elpup, elp, eup, mlpup, exres, vsav, vsavmod,
		vflat, EfromT, TfromE, U0, Cs, tcpol, tcpol2, critdensfromT,
		normdens, tcdirect, tcpostsolve, ns, tc0, ncrit,
		ndeopts = {Method -> {"PDEDiscretization" -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> 10^-3}}, "Eigensystem" -> {"Arnoldi", "MaxIterations" -> 10^12}, "VectorNormalization" -> None}},
		nintopts = {Method -> "GlobalAdaptive", MinRecursion -> 10, MaxRecursion -> 10^6}
  },
 (memi = N@UnitConvert[Quantity[MemoryInUse[], "Bytes"], "Gigabytes"]);
 (* Solve the Schrodinger equation *)
 res = med[Unevaluated[Table[
     NDEigensystem[
      {
       Evaluate[(tdenum /. angm -> mind)],
       DirichletCondition[g[s] == 0, s == Pi/2]
       },
      g[s],
      {s, 0, Pi/2},
      nmax - mind,
      Evaluate@FilterRules[{ndeopts}, Options[NDEigensystem]]
      ],
     {mind, 0, nmax - 1}
     ]]];
 {evs = UnitConvert[Quantity[res[[1, ;; , 1]] - shift, "Hartrees"], 
    "Millielectronvolts"], (efs = res[[1, ;; , 2]]), 
  ndestats = res[[2]]};
 (* Transform and normalize the EFs *)
 tefs = Table[
   Function[{r}, Evaluate[Head[efs[[i, j]]][ArcTan[r/c]]]], {i, 
    nmax}, {j, 1, nmax + 1 - i}];
 {norms, normstats} = med[Unevaluated[Table[
     NIntegrate[
      r*(tefs[[i, j]][r])^2,
      {r, 0, rmax},
      Evaluate@FilterRules[{nintopts}, Options[NIntegrate]]
      ],
     {i, nmax}, {j, 1, nmax + 1 - i}
     ]]];
 ntefs = Table[
   Function[{r}, Evaluate[tefs[[i, j]][r]/Sqrt[norms[[i, j]]]]],
   {i, nmax}, {j, 1, nmax + 1 - i}
   ];
 {normcheck, normcheckstats} = med[Unevaluated[Table[
     Table[
      NIntegrate[r*ntefs[[i, j]][r]*ntefs[[i, jp]][r], {r, 0, rmax}, 
       Evaluate@FilterRules[{nintopts}, Options[NIntegrate]]],
      {j, 1, nmax + 1 - i}, {jp, 1, nmax + 1 - i}
      ],
     {i, nmax}
     ]]];
 (* Reorganize the EVs and ntEFs to make other calculations more \
straightforward *)
 esform = <|Table[
    tssf["`1`,`2`", i, i - j] -> {evs[[i - j + 1, j]], 
      ntefs[[i - j + 1, j]]},
    {i, nmax}, {j, i, 1, -1}
    ]|>;
 evform = Table[evs[[i - j + 1, j]], {i, nmax}, {j, i, 1, -1}];
 efform = Table[ntefs[[i - j + 1, j]], {i, nmax}, {j, i, 1, -1}];
 rsq = Table[
   UnitConvert[
    Quantity[
     Sqrt[NIntegrate[(r^3)*(efform[[n, m, 2]][r])^2, {r, 0, rmax}, 
       Evaluate@FilterRules[{nintopts}, Options[NIntegrate]]]], 
     "BohrRadius"], "Nanometers"], {n, nmax}, {m, n}];
 xbhrad = 
  Table[UnitConvert[
    Quantity[
     NMaximize[{Abs[r*efform[[n, m, 2]][r]], 0 <= r <= 3000}, r, 
      Method -> "SimulatedAnnealing"], "BohrRadius"], 
    "Nanometers"], {n, nmax}, {m, n}];
 (* EXCITON OPTICAL PROPERTIES *)
 f0 = Table[(* THIS IS f0 for ONE transition *)
   2*QuantityMagnitude[
     evform[[n, 1]] - 
      evform[[1, 1]]]*(NIntegrate[
       r*efform[[1, 0]][r]*efform[[n, 1]][r], {r, 0, rmax}, 
       Evaluate@FilterRules[{nintopts}, Options[NIntegrate]]])^2,
   {n, 2, nmax}];
 alcon = (Pi*cee^2)/(2*ce0*Sqrt[kappa]*mu*ccc);
 alpha = Table[ (* ALPHA AND AFAC HAVE A FACTOR OF TWO SO WERE CONSIDERING THE SYMMETRIC TRANSITIONS TO m = +/- 1 *)
   Function[{nx, gamma}, 
    Evaluate[
     2*UnitConvert[
       alcon*(Quantity[nx, "Meters"^-2]/l)*
        f0[[i]]*2/Quantity[gamma, "Seconds"^-1]], "Meters"^-1]],
   {i, nmax - 1}];
 afac = Table[
   Function[{nx, gamma}, 
    Evaluate[
     1 - Exp[-2*alcon*(Quantity[nx, "Meters"^-2])*
        f0[[i]]*2/Quantity[gamma, "Seconds"^-1]]]],
   {i, nmax - 1}];
	out = Dataset[<|
		"evs" -> evform,
		"efs" -> efform,
		"rsq" -> rsq,
		"xbohr" -> xbhrad,
		"f0" -> f0,
		"alpha" -> alpha,
		"afac" -> afac,
		"miscres" -> <|
		"efs" -> efs,
		"tefs" -> tefs,
		"norms" -> norms,
		"normcheck" -> normcheck
		|>,
		"params" -> <|
			"inputs" -> <|"mu" -> mu, "chi" -> chi, "l" -> l, "kappa" -> kappa, "D" -> z, "cscale" -> c, "shift" -> shift|>,
			"calcs" -> <|"ndeopts" -> ndeopts, "nintopts" -> nintopts|>,
			"tde" -> tde,
			"tdenum" -> tdenum
		|>,
		"stats" -> <|
			"nde" -> ndestats,
			"norm" -> normstats,
			"normcheck" -> normcheckstats,
			"meminit" -> memi,
			"memfinl" -> memf,
			"memmax" -> memx
		|>
	|>]
]
(* /*}}}*/*)

(* /*{{{*/ Update ndesscale for new approach *)
ndesscale[nmax_, mx_, my_, pot_, kappa_, chi_, d_, eps_, c_, maxi_, vn_]:=Module[
	{
		tds,
		tdstrans,
		shift=10,
		stranseps=N[(Pi/2)-$MachineEpsilon],
		evs,efs
	},
	tds = (-(1/2)*((1/mx)*D[f[x, y], {x, 2}, {y, 0}] + (1/my)* D[f[x, y], {x, 0}, {y, 2}]) + pot[kappa][d, chi][x, y]*f[x, y] + shift*f[x, y]);
	tdstrans = Simplify[tds /. {f -> (u[ArcTan[#1/c], ArcTan[#2/c]] &)} /. {x -> (c*Tan[\[Xi]]), y -> (c*Tan[\[Psi]])}, {(Pi/2) > \[Xi] > -(Pi/2), (Pi/2) > \[Psi] > -(Pi/2)}];
	{evs,efs}=NDEigensystem[
		{
			tdstrans,
			DirichletCondition[u[\[Xi], \[Psi]] == 0,Abs[\[Xi]] == (stranseps) || Abs[\[Psi]] == (stranseps)]
		},
	 u[\[Xi], \[Psi]],
	 {\[Xi], \[Psi]} \[Element] Rectangle[{-stranseps, -stranseps}, {stranseps, stranseps}],
	 nmax,
	 Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->eps}}},"Eigensystem"->{"Arnoldi","MaxIterations"->maxi},"VectorNormalization"->vn}
	];
	{
		evs-shift,
		Table[
			Function[{x,y},
				Evaluate[
					Head[efs[[i]]][ArcTan[x/c],ArcTan[y/c]]
				]
			],
			{i,nmax}
		]
	}
]
(* /*}}}*/*)

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

