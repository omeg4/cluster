(* Set up conversion factors and general constants *)
e = 1;
H2J = 4.35974417 * (10^-18);
T2S = 2.41884326505 * (10^-17);
B2nm = 0.052917721092;
qbhr = Quantity["BohrRadius"];
qnm = Quantity["Nanometers"];
qmev = Quantity["Millielectronvolts"];
qev = Quantity["Electronvolts"];
qhart = Quantity["Hartrees"];
HH = Quantity["PlanckConstant"];
HB = Quantity["ReducedPlanckConstant"]
CC = Quantity["SpeedOfLight"];
H2eV = 27.21138602;
na = (5 * 10^15) * (10^-18);
damp = (10^13) * (T2S);
lBN = (1 / B2nm) * 0.333;
lphos = (1 / B2nm) * 0.541;


(* Set up phosphorene parameters *)
mey = { 1.2, 1.3, 0.7527, 1.12 };
mex = { 0.17, 0.1, 0.1990, 0.17 };
mhy = { 5.0, 2.8, 5.3525, 6.35 };
mhx = { 0.1, 0.2, 0.1678, 0.15 };
masses = { mey, mex, mhy, mhx };
mus = Table[
	{
		mex[[i]]*mhx[[i]]/(mex[[i]] + mhx[[i]]),
		mey[[i]]*mhy[[i]]/(mey[[i]] + mhy[[i]])
	},
	{i, Length[mey]}
];
chiphos = QuantityMagnitude[Quantity[0.41,"Nanometers"],"BohrRadius"];
(* rho = 2 * Pi * chiphos / kappa;*)

(* Define the potentials *)
VCoul[k_][d_,chi_][x_,y_]:= -1 / (k*Sqrt[x^2 + y^2 + d^2])
VKeld[k_][d_,chi_][x_,y_] := - (1/(4*chi)) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/(2*Pi*chi/k)] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/(2*Pi*chi/k)])
VCho[eps_][d_,chi_][x_,y_]:= (1/(2*eps*(d^3)))*(x^2+y^2) - (1/(eps*d))
VKho[kappa_][d_,chi_][x_,y_]:=-Pi/(4*kappa*((2*Pi*chi/(2*kappa))^2)*d)*(StruveH[-1,d/((2*Pi*chi/(2*kappa)))]-BesselY[-1,d/((2*Pi*chi/(2*kappa)))])
hoCpsi[m_,eps_,d_][n_][x_]:=With[{a=Sqrt[1/Sqrt[(2*m*(1/(2*eps*d^3)))]]},Sqrt[1/(Sqrt[Pi]*a)]*1/Sqrt[(2^n)*Factorial[n]]*Exp[-(x^2)/(2*(a^2))]*HermiteH[n,x/a]]
hoCen[mx_,my_,eps_,d_][n_,m_]:=With[{gamma=1/(2*eps*d^3),v0=(1/(eps*d))},Sqrt[2*gamma/mx]*(n+1/2)+Sqrt[2*gamma/my]*(m+1/2)-v0]

(* Quick function that returns 0 for direct ("d" = -1 for direct) but Nhbn*lbn + lphos if "d" >= 0 *)
rightd[d_]:=If[d==-1,0,d*lBN + lphos]

CompPhos2[{mx_,my_,nhbn_,pot_,kappa_},eps_,nmax_]:=Module[
	{
		chi = chiphos,
		rho = chiphos/kappa,
		shift = 10,
		d = rightd[nhbn],
		ev,
		ef,
		norms,
		tds,
		tdstrans
	},
	s=1/eps;
	tds=- (1/2) * ( (1/mx)*D[f[x,y],{x,2},{y,0}] + (1/my)*D[f[x,y],{x,0},{y,2}]) + pot[kappa][d,chi][x,y] * f[x,y] + shift*f[x,y];
	tdstrans=Simplify[tds /. {
			f -> (u[ArcTan[#1],ArcTan[#2]]&)
		} /. {
			x -> (Tan[\[Xi]]),
			y->(Tan[\[Psi]])
		} , {
			(Pi/2)>\[Xi]>-(Pi/2),
			(Pi/2)>\[Psi]>-(Pi/2)
		}];
	{ev, ef} = NDEigensystem[
			{
				tdstrans,
				DirichletCondition[u[\[Xi],\[Psi]]==0,Abs[\[Xi]]==(Pi/2-eps)||Abs[\[Psi]]==(Pi/2-eps)]
			},
			u[ \[Xi], \[Psi]],
			{\[Xi], \[Psi]} \[Element] Rectangle[{-Pi/2+eps,-Pi/2+eps},{Pi/2-eps,Pi/2-eps}],
			nmax,
			Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->eps}}},"Eigensystem"->{"Arnoldi","MaxIterations"->10^7}}
	];
	{
		{{mx, my}, {nhbn, d}, pot, kappa, chi, eps},
		{
			UnitConvert[Quantity[ev - shift,"Hartrees"],"Millielectronvolts"],
			(* NormalizeEF[(#/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),s]&/@ef*)
			(# / Sqrt[ xysubNInt[ #, #, 1, s] ])&/@(ef/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})
		}
	}
]

(* freezing this equation because i want to test something re: conjugates *)
(* NormalizeEF[EF_,s_] := Module[*)
(*   {norm},*)
(*   norm = Quiet[NIntegrate[(EF)^2,{x,-s,s},{y,-s,s}, MinRecursion -> 5,*)
(*     MaxRecursion -> 20],{NIntegrate::slwcon,NIntegrate::eincr,General::stop}];*)
(*   EF/Sqrt[norm]*)
(* ]*)

NormalizeEF[EF_,s_] := Module[
  {norm},
  norm = Quiet[NIntegrate[Conjugate[EF]*EF,{x,-s,s},{y,-s,s}, MinRecursion -> 5,
    MaxRecursion -> 20],{NIntegrate::slwcon,NIntegrate::eincr,General::stop}];
  EF/Sqrt[norm]
]

xysubNInt[efi_,eff_,op_,s_]:=Chop[N[Quiet[NIntegrate[
	Conjugate[eff]*op*efi,
	{x,-s,s}, {y,-s,s},
	MinRecursion -> 5, MaxRecursion -> 20
	],{NIntegrate::slwcon,NIntegrate::eincr,General::stop}]]
]

(* {{{ *)
makeproc[eps_:(10^-3),nmax_:10,mutab_:mus,ntab_:Table[i,{i,0,10}],kaptab_:{1,4.89}]:=Module[
	{
		mun=Dimensions[mutab] [[1]],
		nn=Length[ntab],
		kn=Length[kaptab],
		DRare,IRare,DRtime,IRtime,
		s = 1/eps,
		muth=Function[{rare, th},	rare[[ 1, 1, 1 ]]*Cos[th]^2 + rare[[ 1, 1, 2 ]]*Sin[th]^2],
		dpmeth=Function[{rare, i, j, th}, (rare[[ 3, i, j ]]*Cos[th]^2 + rare[[ 4, i, j ]]*Sin[th]^2)^2],
		etrfunc=Function[{rare, ni, nf}, rare[[ 2, nf ]] - rare[[ 2, ni ]] ]
	},
	DRtime = AbsoluteTiming[
		DRare = Flatten[
			Table[
				Module[
					{
						raw
					},
					raw=CompPhos2[{mutab[[ mui, 1 ]], mutab[[ mui, 2 ]], -1, pot, kaptab[[ki]]}, eps, nmax];
					Parallelize@{
						raw[[1]],
						raw[[2, 1]],
						Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], x, s ], {i, nmax}, {j, i, nmax}],
						Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], y, s ], {i, nmax}, {j, i, nmax}],
						Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], 1, s ], {i, nmax}, {j, i, nmax}],
						Dimensions[raw]
					}
				],
				{mui, mun}, {pot, {VKeld, VCoul}}, {ki, 2}
			],
			{{2},{3},{1}}
		]
	] [[1]];
	IRtime = AbsoluteTiming[
		IRare = Flatten[
			Table[
				Module[
					{
						raw
					},
					raw=CompPhos2[{mutab[[ mui, 1 ]], mutab[[ mui, 2 ]], ntab[[ni]], pot, 4.89}, eps, nmax];
					Parallelize@{
						raw[[1]],
						raw[[2, 1]],
						Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], x, s ], {i, nmax}, {j, i, nmax}],
						Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], y, s ], {i, nmax}, {j, i, nmax}],
						Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], 1, s ], {i, nmax}, {j, i, nmax}],
						Dimensions[raw]
					}
				],
				{ni, nn}, {mui, mun}, {pot, {VKeld, VCoul}}
			],
			{{3},{1},{2}}
		]
	] [[1]];
	assoctime = AbsoluteTiming[ assoc = Join[
		Association[ (* First level: generic stats and parameters *)
			"lbn" -> UnitConvert[Quantity[lBN,"BohrRadius"],"Nanometers"],
			"lphos" -> UnitConvert[Quantity[lphos, "BohrRadius"], "Nanometers"],
			"chi2d" -> Quantity[0.41, "Nanometers"],
			"dims" -> {
				ToString@StringForm["Nmu = ``", mun],
				ToString@StringForm["Nkappa = ``", kn],
				ToString@StringForm["Nhbn = ``", nn],
				ToString@StringForm["Nmax = ``", nmax]
			},
			"stats" -> {
				ToString@StringForm["DRtime = ``", DRtime],
				ToString@StringForm["IRtime = ``", IRtime]
			}
		],
		Map[
			<|Table[ #["kornkey"] [[korn]] ->
				<|Table[ ToString@StringForm["mu``",mu] ->
					<|
						"mux" -> #["rare"] [[ korn, mu, 1, 1, 1]],
						"muy" -> #["rare"] [[ korn, mu, 1, 1, 2]],
						"params" -> #["rare"] [[ korn, mu, 1 ]],
						"raredims" -> Dimensions[#["rare"]],
						"muth" -> Interpolation[
							Table[
								{
									th,
									muth[ #["rare"] [[korn, mu]], th ]
								},
								{th, 0, 2 * Pi, Pi/100}
							]
						],
						"evs" -> #["rare"] [[ korn, mu, 2 ]],
						"etr" -> Table[ etrfunc[ #["rare"] [[korn, mu]], i, j ], {i, nmax}, {j, i, nmax}],
						"dpmeth" -> Table[
							Interpolation[
								Table[
									{
										th,
										dpmeth[ #["rare"] [[korn, mu]], i, j, th ]
									},
									{th, 0, 2 * Pi, Pi/100}
								]
							],
							{i, nmax}, {j, 1, nmax - i + 1}
						],
						"f0th" -> Table[
							Interpolation[
								Table[
									{
										th,
										2 * muth[ #["rare"] [[korn, mu]], th ] * dpmeth[ #["rare"] [[korn, mu]], i, j, th ] * QuantityMagnitude[etrfunc[ #["rare"] [[korn, mu]], i, j ], "Hartrees"]
									},
									{th, 0, 2 * Pi, Pi/100}
								]
							],
							{i, nmax}, {j, 1, nmax - i + 1}
						],
						"norm" -> #["rare"] [[korn, mu, 5]]
					|>,
					{mu,mun}
				]|>,
			{korn, #["korniter"]}
			]|>&,
			<|
				"K" -> <|
					"D" -> <|
						"rare" -> DRare[[1]],
						"korn" -> ToString[k],
						"korniter" -> kn,
						"kornkey" -> {"FS", "Enc"}
					|>,
					"I" -> <|
						"rare" -> IRare[[1]],
						"korn" -> "Nhbn",
						"korniter" -> nn,
						"kornkey" -> Table[ToString@StringForm["n``", n], {n, nn}]
					|>
				|>,
				"C" -> <|
					"D" -> <|
						"rare" -> DRare[[2]],
						"korn" -> ToString[k],
						"korniter" -> kn,
						"kornkey" -> {"FS", "Enc"}
					|>,
					"I" -> <|
						"rare" -> IRare[[2]],
						"korn" -> "Nhbn",
						"korniter" -> nn,
						"kornkey" -> Table[ToString@StringForm["n``", n], {n, nn}]
					|>
				|>
			|>,
			{2}
		]
	]] [[1]];
	assoc["stats"] = Append[ assoc["stats"], ToString@StringForm[ "assoctime = ``", assoctime ] ];
	assoc
] (* function (Module) ends here *)
(* }}} *)

(* functions with options *)

CompPhos2opts[{mx_,my_,nhbn_,pot_,kappa_},eps_,nmax_,wp_,ndeopts:OptionsPattern[],niopts:OptionsPattern[]]:=Module[
	{
		chi = chiphos,
		rho = chiphos/kappa,
		shift = 10,
		d = rightd[nhbn],
		s=1/eps,
		ev,
		ef,
		norms,
		tds,
		tdstrans,
		stranseps=SetPrecision[(Pi/2)-eps, wp]
	},
	tds=(-(1/2) * (
			(1/mx)*D[f[x,y],{x,2},{y,0}] +
			(1/my)*D[f[x,y],{x,0},{y,2}]
		) +
		pot[kappa][d,chi][x,y] * f[x,y] +
		shift*f[x,y]
	);
	tdstrans=Simplify[tds /. {
			f -> (u[ArcTan[#1],ArcTan[#2]]&)
		} /. {
			x -> (Tan[\[Xi]]),
			y->(Tan[\[Psi]])
		} , {
			(Pi/2)>\[Xi]>-(Pi/2),
			(Pi/2)>\[Psi]>-(Pi/2)
		}];
	ndetime = AbsoluteTiming[{ev, ef} = NDEigensystem[
			{(* N[tdstrans,wprec],*)
				tdstrans,
				DirichletCondition[u[\[Xi],\[Psi]]==0,Abs[\[Xi]]==(stranseps)||Abs[\[Psi]]==(stranseps)]
			},
			u[ \[Xi], \[Psi]],
			{\[Xi], \[Psi]} \[Element] Rectangle[{-stranseps,-stranseps},{stranseps,stranseps}],
			nmax,
			(* Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->eps}}},"Eigensystem"->{"Arnoldi","MaxIterations"->10^7}}*)
			Evaluate@FilterRules[{ndeopts},Options[NDEigensystem]]
		]
	][[1]];
	{
		{{mx, my}, {nhbn, d}, pot, kappa, chi, eps, ndetime},
		{
			UnitConvert[Quantity[ev - shift,"Hartrees"],"Millielectronvolts"],
			Head[#]&/@ef
		}
	}
]

xysubNIntopts[efi_,eff_,op_,s_,nopts:OptionsPattern[]]:=Chop[
	N[
		Quiet[
			NIntegrate[
				Conjugate[eff]*op*efi,
				{x,-s,s}, {y,-s,s},
				Evaluate@FilterRules[{nopts},Options[NIntegrate]]
			],
		{NIntegrate::slwcon,NIntegrate::eincr,General::stop}]
	]
]

headNInt[efi_, eff_, op_, s_, niopts:OptionsPattern[]]:=NIntegrate[
	Conjugate[ eff[ ArcTan[x], ArcTan[y] ] ] * op * efi[ ArcTan[x], ArcTan[y] ],
	{x,-s,s},{y,-s,s},
	Evaluate@FilterRules[ {niopts}, Options[ NIntegrate ] ]
]

(* FUNCTION FOR GETTING METHODS AND OPTIONS *)
 getList[name_String] := Module[{options, idx}, options = Names[name <> "`*"];
   options = ToExpression /@ options;
   options = {#, Options[#]} & /@ options;
   idx = Range[Length[options]];
   options = {#[[1]], TableForm[#[[2]]]} & /@ options;
   options = Insert[options[[#]], #, 1] & /@ idx;
   options = Insert[options, {"#", "Option", "Options to this option"}, 1]
];

(*
	USE IT LIKE THIS

	Grid[getList["NIntegrate"], Frame -> All, Alignment -> Left, FrameStyle -> Directive[Thickness[.005], Gray]]

*)

(* This function only searches for Method options *)
getList2[name_String] := Module[{options, idx,z1,z2},
   options = Names[name <> "`*"];
   options = ToExpression /@ options;
   options = Flatten[Last@Reap@Do[z1 = Options[options[[i]]];
        If[z1 != {},
         z2 = Cases[z1, Rule["Method", x_] :> Method -> x];
         If[Length[z2] != 0 , Sow[{options[[i]], z2}]]
         ],
        {i, Length[options]}
        ], 1];
   (* rest for formatting*)
   idx = Range[Length[options]];
   options = {#[[1]], TableForm[#[[2]]]} & /@ options;
   options = Insert[options[[#]], #, 1] & /@ idx;
   options = Insert[options, {"#", "Option", "Options to this option"}, 1]
];

(*
	USAGE

	r = getList2["NDSolve"];
	Grid[r, Frame -> All, Alignment -> Left]

*)

(* New Series of functions *)
ndef[nmax_,mx_,my_,pot_,kappa_,chi_,d_,eps_,maxi_,vn_]:=Module[
	{
		tds,
		tdstrans,
		shift=10,
		stranseps=N[(Pi/2)-$MachineEpsilon],
		evs,efs
	},
	tds = (-(1/2)*((1/mx)*D[f[x, y], {x, 2}, {y, 0}] + (1/my)* D[f[x, y], {x, 0}, {y, 2}]) + pot[kappa][d, chi][x, y]*f[x, y] + shift*f[x, y]);
	tdstrans = Simplify[tds /. {f -> (u[ArcTan[#1], ArcTan[#2]] &)} /. {x -> (Tan[\[Xi]]), y -> (Tan[\[Psi]])}, {(Pi/2) > \[Xi] > -(Pi/2), (Pi/2) > \[Psi] > -(Pi/2)}];
	EvaluationData[
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
			Head/@efs
		}
	]
]

ndefmep[nmax_,mx_,my_,pot_,kappa_,chi_,d_,eps_,maxi_,vn_]:=Module[
	{
		tds,
		tdstrans,
		shift=10,
		stranseps=N[(Pi/2)-$MachineEpsilon],
		evs,efs
	},
	tds = (-(1/2)*((1/mx)*D[f[x, y], {x, 2}, {y, 0}] + (1/my)* D[f[x, y], {x, 0}, {y, 2}]) + pot[kappa][d, chi][x, y]*f[x, y] + shift*f[x, y]);
	tdstrans = Simplify[tds /. {f -> (u[ArcTan[#1], ArcTan[#2]] &)} /. {x -> (Tan[\[Xi]]), y -> (Tan[\[Psi]])}, {(Pi/2) > \[Xi] > -(Pi/2), (Pi/2) > \[Psi] > -(Pi/2)}];
	EvaluationData[
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
			Head/@efs
		}
	]
]

normraw[ef_,s_,opts:OptionsPattern[]]:=Module[
	{nevadata, norm, normef},
	nevadata=EvaluationData[
		norm=NIntegrate[
			Conjugate[ ef[ ArcTan[x], ArcTan[y] ] ] * ef[ ArcTan[x], ArcTan[y] ],
			{x,-s,s},{y,-s,s},
			Evaluate@FilterRules[{opts},Options[NIntegrate]]
		]
	];
	normef=Function[{x, y}, ef[ ArcTan[x], ArcTan[y] ] / Sqrt[ norm ]];
	{
		normef,
		nevadata
	}
]

apnint[ef1_, ef2_, s_, op_, opts:OptionsPattern[]]:=EvaluationData[
	NIntegrate[
		Conjugate[ ef1[ x, y ] ] * op * ef2[ x, y ],
		{x,-s,s},{y,-s,s},
		Evaluate@FilterRules[{opts},Options[NIntegrate]]
	]
]

triassnint[ efs_, s_, op_, opts:OptionsPattern[]]:=Association[
	Table[
		$Messages={};
		ToString@StringForm["i``",i]->Association[
			Table[
				ToString@StringForm["j``",j]->apnint[
					efs[[ j ]], efs[[ i ]], s, op, opts
				],
				{j,i,Length@efs}
			]
		],
		{i,Length@efs}
	]
]

(* rewrite again *)
(* start with "modified eval data" function *)
med[f_]:=({#["Result"],KeyTake[#,{"AbsoluteTiming","MessagesText","Timing"}]}&[KeyTake[EvaluationData[f],{"Result","AbsoluteTiming","MessagesText","Timing"}]])
med[f_,x___]:=({#["Result"],KeyTake[#,{"AbsoluteTiming","MessagesText","Timing"}]}&[KeyTake[EvaluationData[f[x]],{"Result","AbsoluteTiming","MessagesText","Timing"}]])


ndes[nmax_, mx_, my_, pot_, kappa_, chi_, d_, eps_, maxi_, vn_]:=Module[
	{
		tds,
		tdstrans,
		shift=10,
		stranseps=N[(Pi/2)-$MachineEpsilon],
		evs,efs
	},
	tds = (-(1/2)*((1/mx)*D[f[x, y], {x, 2}, {y, 0}] + (1/my)* D[f[x, y], {x, 0}, {y, 2}]) + pot[kappa][d, chi][x, y]*f[x, y] + shift*f[x, y]);
	tdstrans = Simplify[tds /. {f -> (u[ArcTan[#1], ArcTan[#2]] &)} /. {x -> (Tan[\[Xi]]), y -> (Tan[\[Psi]])}, {(Pi/2) > \[Xi] > -(Pi/2), (Pi/2) > \[Psi] > -(Pi/2)}];
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
		Head/@efs
	}
]

ndespolar[nmax_, mx_, my_, pot_, kappa_, chi_, d_, eps_, maxi_, vn_]:=Module[
	{
		tdspt,
		shift=10,
		pio2 = Pi/2,
		stranseps=N[(Pi/2)-$MachineEpsilon],
		evs,efs
	},
	tdspt=FullSimplify[
	 FullSimplify[
			Simplify[
				 ((-1/(2*mx))*D[f[x, y], {x, 2}, {y, 0}]) + ((-1/(2*my))*
						D[f[x, y], {x, 0}, {y, 2}]) + 
					VKeld[kappa][d, chi][x, y]*f[x, y] + 10*f[x, y]
					] /. {
					f -> (u[Sqrt[#1^2 + #2^2], ArcTan[#2/#1]] &)} /. {
				x -> r*Cos[\[Theta]], y -> r*Sin[\[Theta]]
				},
			{-pio2 < \[Theta] < pio2, 0 <= r < Infinity}
		] /. {u -> (g[ArcTan[#1], #2] &)} /. {r -> Tan[\[Psi]]},
		{0 <= \[Psi] < Pi/2}
	];
	{evs,efs}=NDEigensystem[
		{
			tdspt,
			DirichletCondition[g[\[Psi], \[Theta]] == 0, (\[Psi] == (pio2)) && 0 < \[Theta] <= 2*Pi],
			(* DirichletCondition[(Cos[\[Psi]]^2)*D[g[\[Psi], \[Theta]],\[Psi]] == 0, (\[Psi] == (pio2)) && 0 < \[Theta] <= 2*Pi],*)
			PeriodicBoundaryCondition[g[\[Psi],\[Theta]], \[Theta]==0, TranslationTransform[{0,2*Pi}]]
		},
	 g[\[Psi], \[Theta]],
	 {\[Psi], \[Theta]} \[Element] Rectangle[{0, 0}, {pio2, 2*Pi}],
	 nmax,
	 Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->eps}}},"Eigensystem"->{"Arnoldi","MaxIterations"->maxi},"VectorNormalization"->vn}
	];
	{
		evs-shift,
		Head/@efs
	}
]

ning[ef1_, ef2_, op_, s_]:=NIntegrate[
	Conjugate[ ef1 ] * op * ef2,
	{x,-s,s},{y,-s,s},
	Method->"LocalAdaptive"
]

nestedgridWheaders[data_,{innerrow_, innercolumn_,innergropts : OptionsPattern[]},{outerrow_, outercolumn_,outergropts : OptionsPattern[]}]:=Grid[
  Join[
   {PadLeft[outerrow, Length[outerrow] + 1, " "]} // Transpose,
   Join[
    {outercolumn},
    Map[
     Grid[
       Join[
        {PadLeft[innerrow, Length[innerrow] + 1, " "]} // Transpose,
        Join[
         {innercolumn},
         #
         ],
        2
        ],
       Dividers -> {{False, True}, {False, True}},
       ItemStyle -> 
        Directive[FontSize -> 14, Black, FontFamily -> "Arial"],
       Evaluate@FilterRules[{innergropts}, Options[Grid]]
       ] &,
     data,
     {2}
     ]
    ],
   2
   ],
  Dividers -> {{False, True}, {False, True}},
  ItemStyle -> Directive[FontSize -> 14, Black, FontFamily -> "Arial"],
  Evaluate@FilterRules[{outergropts}, Options[Grid]]
  ]
