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
				{ 0 < \[Theta] < 2 * Pi, 0 <= r < Infinity }
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

(* /*{{{*/ 1D Polar NDE solver PLUS optical properties of exciton *)
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

Table[
	Export[
		tssf[ "1dp_l-`1`_chi-`2`_kappa-`3`_mu-`4`.m", l, chi, kappa, mu ],
		d1pnde[ { mu, QuantityMagnitude[Quantity[chi,"Angstroms"],"BohrRadius"], QuantityMagnitude[Quantity[l,"Angstroms"],"BohrRadius"], kappa } ]
	],
	{ l, 4, 8, 0.5 }
	{ chi, 4, 10, 0.5 }
	{ kappa, 1, 7, 0.5 }
	{ mu, 0.1, 1.5, 0.05 }
];

Quit[]
