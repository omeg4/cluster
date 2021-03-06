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

gridWheaders[data_, rowheads_, columnheads_, gropts : OptionsPattern[]] := Grid[
  Join[
   {PadLeft[rowheads, Length[rowheads] + 1, " "]} // Transpose,
   Join[
    {columnheads},
    data
    ],
   2
   ],
  Dividers -> {{False, True}, {False, True}},
  ItemStyle -> Directive[FontSize -> 20, Black, FontFamily -> "Arial"],
  Evaluate@FilterRules[{gropts}, Options[Grid]]
]


newvarplt[data_,bounds_,frmlbl_,lgndlbl_,colors_,size_:{1600,900},lblsize_:30,plopts:OptionsPattern[]]:=
	{
	Framed@LineLegend[colors,lgndlbl,LegendLayout->{"Row",1},LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],LegendFunction->"Frame"],
	Plot[
		data,
		bounds,
		PlotTheme->"Detailed",
		FrameLabel->frmlbl,
		GridLinesStyle->Directive[Gray,Thin],
		LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
		ImageSize->size,
		AspectRatio->size[[2]]/size[[1]],
		PlotStyle->colors,
		PlotLegends->None,
		Evaluate@FilterRules[{plopts},Options[Plot]]
	]
	}

newlineplt[data_,frmlbl_,lgndlbl_,colors_,size_:{1600,900},lblsize_:30,plopts:OptionsPattern[]]:=
{
Framed@LineLegend[colors,lgndlbl,LegendLayout->{"Row",1},LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],LegendFunction->"Frame",LegendMarkerSize->{75,25}],
ListLinePlot[
data,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
FrameLabel->frmlbl,
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
Evaluate@FilterRules[{plopts},Options[ListLinePlot]]
]
}

newplt[data_,frmlbl_,lgndlbl_,colors_,{markers_,msize_:24,lmsize_:30},size_:{1600,900},lblsize_:30,plopts:OptionsPattern[]]:=
{
Framed@PointLegend[colors,lgndlbl,LegendMarkers->Table[Style[ToString@markers[[i]],lmsize],{i,Length[markers]}],LegendMarkerSize->lmsize,LegendLayout->{"Row",1},LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],LegendFunction->"Frame"],
ListPlot[
data,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
FrameLabel->frmlbl,
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
PlotMarkers->Table[{markers[[i]],msize},{i,Length[markers]}],
Evaluate@FilterRules[{plopts},Options[ListPlot]]
]
}

insetvarplt[data_,bounds_,colors_,ipos_,size_:{800,450},lblsize_:24,plopts:OptionsPattern[]]:=
Inset[
Framed[
Plot[
data,
bounds,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
PlotLegends->None,
Evaluate@FilterRules[{plopts},Options[Plot]]
],
Background->White,FrameStyle->None
],
Scaled[{ipos[[1]],ipos[[2]]}],Scaled[{ipos[[3]],ipos[[4]]}]
]

insetlineplt[data_,colors_,ipos_,size_:{800,450},lblsize_:24,plopts:OptionsPattern[]]:=
Inset[
Framed[
ListLinePlot[
data,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
Evaluate@FilterRules[{plopts},Options[ListLinePlot]]
],
Background->White,FrameStyle->None
],
Scaled[{ipos[[1]],ipos[[2]]}],Scaled[{ipos[[3]],ipos[[4]]}]
]

insetlistplt[data_,colors_,markers_,ipos_,size_:{800,450},lblsize_:24,plopts:OptionsPattern[]]:=
	Inset[
		Framed[
				ListPlot[
				data,
				PlotTheme->"Detailed",
				GridLinesStyle->Directive[Gray,Thin],
				LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
				ImageSize->size,
				AspectRatio->size[[2]]/size[[1]],
				PlotStyle->colors,
				PlotMarkers->Table[{markers[[i]],lblsize},{i,Length[markers]}],
				Evaluate@FilterRules[{plopts},Options[ListPlot]]
			],
		Background->White,FrameStyle->None
		],
	Scaled[{ipos[[1]],ipos[[2]]}],Scaled[{ipos[[3]],ipos[[4]]}]
]

SetDirectory["$(pwd)"];
day:=DateList[TimeZone->"America/New_York"][[3]];
month:=DateList[TimeZone->"America/New_York"][[2]];
hour:=DateList[TimeZone->"America/New_York"][[4]];
minute:=DateList[TimeZone->"America/New_York"][[5]];
gettime:=DateList[TimeZone->"America/New_York"];
ee=Quantity["ElementaryCharge"];
cc=Quantity["SpeedOfLight"];
e0=Quantity["ElectricConstant"];
niopts={Method->"GlobalAdaptive",MinRecursion->1000,MaxRecursion->10^6};
s=10^4;
c=10;
nmax=12;
eps=10^-3;
maxiter=10^15;
mindmax = 4;
dindmax = 8;
(* This chunk is just re-doing my long Keldysh calculations (need to check that the optical results are the same, otherwise I'm kinda boned *)
Table[
	Module[
		{
			rev,ref,ndestats,nef,norms,normstats,onc,oncstats,fox,foxstats,foy,foystats,
			paramsassoc,alphax,alphay,afax,afay,
			starttime,
			mx=mus[[muind,1]],
			my=mus[[muind,2]],
			pot=VKeld,
			kappa=4.89,
		},
		starttime=gettime;
		ToString[StringForm["At `3`h`4`m on `5`/`6`, starting calculation for (RK, hbn-enc) mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Solve the schrodinger equation *)
		{{rev,ref},ndestats}=med[Unevaluated[ndesscale[nmax,mx,my,pot,kappa,chiphos,rightd[d],eps,c,maxiter,None]]];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, NDE calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* adding a third list element to the "NDE results" where I can store computational parameters *)
		paramsassoc=Association[
			"mu index" -> muind,
			"d" -> d,
			"dval" -> rightd[d],
			"mux" -> mx,
			"muy" -> my,
			"pot" -> pot,
			"kappa" -> kappa,
			"chi" -> chiphos,
			"nde eps" -> eps,
			"c rescale" -> c,
			"maxiter" -> maxiter,
			"int box s" -> s,
			"int opts" -> ToString[niopts]
		];
		(* save the data *)
		Export[ToString@StringForm["ndeevs_rk_he_mu`1`_d`2`.m",muind,(d+1)],rev];
		Export[ToString@StringForm["ndeefs_rk_he_mu`1`_d`2`.m",muind,(d+1)],ref];
		Export[ToString@StringForm["ndeparams_rk_he_mu`1`_d`2`.m",muind,(d+1)],paramsassoc];
		Export[ToString@StringForm["ndestats_rk_he_mu`1`_d`2`.m",muind,(d+1)],ndestats];
		(* Calculate the normalization constants *)
		{norms,normstats}=Flatten[Table[med[Unevaluated[NIntegrate[
							(Conjugate[#]*#)&[ref[[n]][x,y]],
						{x,-s,s},{y,-s,s},
						Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax}],{{2}}];
		Export[ToString@StringForm["normconsts_rk_he_mu`1`_d`2`.m",muind,(d+1)],norms];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, norms calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Make the normalized EFs *)
		nef = Table[
			Function[{x,y},
				Evaluate[
					ref[[ n ]][ x,y ] / Sqrt[ norms[[ n ]] ]
				]
			],
			{n, nmax}
		];
		(* Save the normalized EFs because why not? *)
		Export[ToString@StringForm["ndenefs_rk_he_mu`1`_d`2`.m",muind,(d+1)],nef];
		(* check orthonormality *)
		{onc,oncstats}=Flatten[Table[med[Unevaluated[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["onc_rk_he_mu`1`_d`2`.m",muind,(d+1)],onc];
		Export[ToString@StringForm["oncstats_rk_he_mu`1`_d`2`.m",muind,(d+1)],oncstats];
		(* calculate f0x *)
		{fox,foxstats}=Flatten[Table[med[Unevaluated[2 * mx * (rev[[m]] - rev[[n]]) * (Abs[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * x * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]^2)]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["fox_rk_he_mu`1`_d`2`.m",muind,(d+1)],fox];
		Export[ToString@StringForm["foxstats_rk_he_mu`1`_d`2`.m",muind,(d+1)],foxstats];
		(* Calculate f0y *)
		{foy,foystats}=Flatten[Table[med[Unevaluated[2 * my * (rev[[m]] - rev[[n]]) * (Abs[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * y * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]^2)]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["foy_rk_he_mu`1`_d`2`.m",muind,(d+1)],foy];
		Export[ToString@StringForm["foystats_rk_he_mu`1`_d`2`.m",muind,(d+1)],foystats];
		(* Calculate alpha x and y *)
		{alphax,alphay}=Table[
			With[
				{
					mu=Quantity[xory[[1]],"ElectronMass"],
					f0xy=xory[[2]],
					l=If[rightd[d]==0,Quantity[lphos,"BohrRadius"],Quantity[2*lphos,"BohrRadius"]]
				},
				Table[
					Function[{nx,gamma},Evaluate[UnitConvert[(( Pi*(ee^2) )/( 2 * e0 * Sqrt[kappa] * mu * cc )) * ( ( Quantity[nx,"Meters"^-2] ) / ( l ) ) * f0xy[[ni,nf]]*(2/Quantity[gamma,"Seconds"^-1]),"Meters"^-1]]],
					{ni,nmax},{nf,1,nmax-ni}
				]
			],{xory,{{mx,fox},{my,foy}}}
		];
		Export[ToString@StringForm["alphax_rk_he_mu`1`_d`2`.m",muind,(d+1)],alphax];
		Export[ToString@StringForm["alphay_rk_he_mu`1`_d`2`.m",muind,(d+1)],alphay];
		(* Calculate afac x and y *)
		{afax,afay}=Table[
			With[
				{
					mu=Quantity[xory[[1]],"ElectronMass"],
					f0xy=xory[[2]]
				},
				Table[
					Function[{nx,gamma},Evaluate[1-Exp[-(( Pi*(ee^2) )/( 2 * e0 * Sqrt[kappa] * mu * cc )) * ( ( Quantity[nx,"Meters"^-2] ) ) * f0xy[[ni,nf]] * (2/Quantity[gamma,"Seconds"^-1])]]],
					{ni,nmax},{nf,1,nmax-ni}
				]
			],{xory,{{mx,fox},{my,foy}}}
		];
		Export[ToString@StringForm["afax_rk_he_mu`1`_d`2`.m",muind,(d+1)],afax];
		Export[ToString@StringForm["afay_rk_he_mu`1`_d`2`.m",muind,(d+1)],afay];
		ToString[StringForm["`3`h`4`m on `5`/`6`: loop complete. Total runtime for this leg is `7`",muind,d,hour,minute,month,day,DateDifference[starttime,gettime,{"Hours","Minutes","Seconds"}]]]>>>mylog.txt;
		ToString[StringForm["-------------------------------------------",muind,d,hour,minute,month,day]]>>>mylog.txt;
	],
	{muind,mindmax},{d,-1,dindmax}
];
ToString[StringForm["============================================================================================",muind,d,hour,minute,month,day]]>>>mylog.txt;
ToString[StringForm["`3`h`4`m on `5`/`6`: RK, h-BN encap done!!",muind,d,hour,minute,month,day]]>>>mylog.txt;
ToString[StringForm["============================================================================================",muind,d,hour,minute,month,day]]>>>mylog.txt;
(* Now do a quick version of the above for ONLY direct excitons in MLBP with kappa = 1, (3.8 + 1)/2, (4.89 + 1)/2 *)
Table[
	Module[
		{
			rev,ref,ndestats,nef,norms,normstats,onc,oncstats,fox,foxstats,foy,foystats,
			paramsassoc,alphax,alphay,afax,afay,
			starttime,
			mx=mus[[muind,1]],
			my=mus[[muind,2]],
			pot=VKeld,
			kappa=kapiter[[1]],
			d=-1,
			envname=ToString[kapiter[[2]]]
		},
		starttime=gettime;
		ToString[StringForm["At `3`h`4`m on `5`/`6`, starting calculation for (RK, `7`) mu=`1`, d=`2`",muind,d,hour,minute,month,day,envname]]>>>mylog.txt;
		(* Solve the schrodinger equation *)
		{{rev,ref},ndestats}=med[Unevaluated[ndesscale[nmax,mx,my,pot,kappa,chiphos,rightd[d],eps,c,maxiter,None]]];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, NDE calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* adding a third list element to the "NDE results" where I can store computational parameters *)
		paramsassoc=Association[
			"mu index" -> muind,
			"d" -> d,
			"dval" -> rightd[d],
			"mux" -> mx,
			"muy" -> my,
			"pot" -> pot,
			"kappa" -> kappa,
			"chi" -> chiphos,
			"nde eps" -> eps,
			"c rescale" -> c,
			"maxiter" -> maxiter,
			"int box s" -> s,
			"int opts" -> ToString[niopts]
		];
		(* save the data *)
		Export[ToString@StringForm["ndeevs_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],rev];
		Export[ToString@StringForm["ndeefs_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],ref];
		Export[ToString@StringForm["ndeparams_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],paramsassoc];
		Export[ToString@StringForm["ndestats_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],ndestats];
		(* Calculate the normalization constants *)
		{norms,normstats}=Flatten[Table[med[Unevaluated[NIntegrate[
							(Conjugate[#]*#)&[ref[[n]][x,y]],
						{x,-s,s},{y,-s,s},
						Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax}],{{2}}];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, norms calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Make the normalized EFs *)
		nef = Table[
			Function[{x,y},
				Evaluate[
					ref[[ n ]][ x,y ] / Sqrt[ norms[[ n ]] ]
				]
			],
			{n, nmax}
		];
		(* Save the normalized EFs because why not? *)
		Export[ToString@StringForm["ndenefs_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],nef];
		(* check orthonormality *)
		{onc,oncstats}=Flatten[Table[med[Unevaluated[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["onc_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],onc];
		Export[ToString@StringForm["oncstats_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],oncstats];
		(* calculate f0x *)
		{fox,foxstats}=Flatten[Table[med[Unevaluated[2 * mx * (rev[[m]] - rev[[n]]) * (Abs[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * x * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]^2)]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["fox_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],fox];
		Export[ToString@StringForm["foxstats_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],foxstats];
		(* Calculate f0y *)
		{foy,foystats}=Flatten[Table[med[Unevaluated[2 * my * (rev[[m]] - rev[[n]]) * (Abs[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * y * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]^2)]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["foy_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],foy];
		Export[ToString@StringForm["foystats_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],foystats];
		(* Calculate alpha x and y *)
		{alphax,alphay}=Table[
			With[
				{
					mu=Quantity[xory[[1]],"ElectronMass"],
					f0xy=xory[[2]],
					l=If[rightd[d]==0,Quantity[lphos,"BohrRadius"],Quantity[2*lphos,"BohrRadius"]]
				},
				Table[
					Function[{nx,gamma},Evaluate[UnitConvert[(( Pi*(ee^2) )/( 2 * e0 * Sqrt[kappa] * mu * cc )) * ( ( Quantity[nx,"Meters"^-2] ) / ( l ) ) * f0xy[[ni,nf]]*(2/Quantity[gamma,"Seconds"^-1]),"Meters"^-1]]],
					{ni,nmax},{nf,1,nmax-ni}
				]
			],{xory,{{mx,fox},{my,foy}}}
		];
		Export[ToString@StringForm["alphax_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],alphax];
		Export[ToString@StringForm["alphay_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],alphay];
		(* Calculate afac x and y *)
		{afax,afay}=Table[
			With[
				{
					mu=Quantity[xory[[1]],"ElectronMass"],
					f0xy=xory[[2]]
				},
				Table[
					Function[{nx,gamma},Evaluate[1-Exp[-(( Pi*(ee^2) )/( 2 * e0 * Sqrt[kappa] * mu * cc )) * ( ( Quantity[nx,"Meters"^-2] ) ) * f0xy[[ni,nf]] * (2/Quantity[gamma,"Seconds"^-1])]]],
					{ni,nmax},{nf,1,nmax-ni}
				]
			],{xory,{{mx,fox},{my,foy}}}
		];
		Export[ToString@StringForm["afax_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],afax];
		Export[ToString@StringForm["afay_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],afay];
		ToString[StringForm["`3`h`4`m on `5`/`6`: loop complete. Total runtime for this leg is `7`",muind,d,hour,minute,month,day,DateDifference[starttime,gettime,{"Hours","Minutes","Seconds"}]]]>>>mylog.txt;
		ToString[StringForm["-------------------------------------------",muind,d,hour,minute,month,day]]>>>mylog.txt;
	],
	{muind,mindmax},{kapiter,{{1,ToString["fs"]},{((3.8+1)/2),ToString["ss"]},{((4.89+1)/2),ToString["hs"]}}}
];
ToString[StringForm["============================================================================================",muind,d,hour,minute,month,day]]>>>mylog.txt;
ToString[StringForm["`3`h`4`m on `5`/`6`: Direct X, RK, different kappa done!!",muind,d,hour,minute,month,day]]>>>mylog.txt;
ToString[StringForm["============================================================================================",muind,d,hour,minute,month,day]]>>>mylog.txt;
(* Do the same as the first big table here for the coulomb potential. Wasn't planning on doing direct excitons here but its a relatively minor inclusion considering how long this will run, so whatever. *)
Table[
	Module[
		{
			rev,ref,ndestats,nef,norms,normstats,onc,oncstats,fox,foxstats,foy,foystats,
			paramsassoc,alphax,alphay,afax,afay,
			starttime,
			mx=mus[[muind,1]],
			my=mus[[muind,2]],
			pot=VCoul,
			kappa=4.89
		},
		starttime=gettime;
		ToString[StringForm["At `3`h`4`m on `5`/`6`, starting calculation for (coul, hbn-enc) mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Solve the schrodinger equation *)
		{{rev,ref},ndestats}=med[Unevaluated[ndesscale[nmax,mx,my,pot,kappa,chiphos,rightd[d],eps,c,maxiter,None]]];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, NDE calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* adding a third list element to the "NDE results" where I can store computational parameters *)
		paramsassoc=Association[
			"mu index" -> muind,
			"d" -> d,
			"dval" -> rightd[d],
			"mux" -> mx,
			"muy" -> my,
			"pot" -> pot,
			"kappa" -> kappa,
			"chi" -> chiphos,
			"nde eps" -> eps,
			"c rescale" -> c,
			"maxiter" -> maxiter,
			"int box s" -> s,
			"int opts" -> ToString[niopts]
		];
		(* save the data *)
		Export[ToString@StringForm["ndeevs_c_he_mu`1`_d`2`.m",muind,(d+1)],rev];
		Export[ToString@StringForm["ndeefs_c_he_mu`1`_d`2`.m",muind,(d+1)],ref];
		Export[ToString@StringForm["ndeparams_c_he_mu`1`_d`2`.m",muind,(d+1)],paramsassoc];
		Export[ToString@StringForm["ndestats_c_he_mu`1`_d`2`.m",muind,(d+1)],ndestats];
		(* Calculate the normalization constants *)
		{norms,normstats}=Flatten[Table[med[Unevaluated[NIntegrate[
							(Conjugate[#]*#)&[ref[[n]][x,y]],
						{x,-s,s},{y,-s,s},
						Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax}],{{2}}];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, norms calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Make the normalized EFs *)
		nef = Table[
			Function[{x,y},
				Evaluate[
					ref[[ n ]][ x,y ] / Sqrt[ norms[[ n ]] ]
				]
			],
			{n, nmax}
		];
		(* Save the normalized EFs because why not? *)
		Export[ToString@StringForm["ndenefs_c_he_mu`1`_d`2`.m",muind,(d+1)],nef];
		(* check orthonormality *)
		{onc,oncstats}=Flatten[Table[med[Unevaluated[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["onc_c_he_mu`1`_d`2`.m",muind,(d+1)],onc];
		Export[ToString@StringForm["oncstats_c_he_mu`1`_d`2`.m",muind,(d+1)],oncstats];
		(* calculate f0x *)
		{fox,foxstats}=Flatten[Table[med[Unevaluated[2 * mx * (rev[[m]] - rev[[n]]) * (Abs[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * x * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]^2)]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["fox_c_he_mu`1`_d`2`.m",muind,(d+1)],fox];
		Export[ToString@StringForm["foxstats_c_he_mu`1`_d`2`.m",muind,(d+1)],foxstats];
		(* Calculate f0y *)
		{foy,foystats}=Flatten[Table[med[Unevaluated[2 * my * (rev[[m]] - rev[[n]]) * (Abs[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * y * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]^2)]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["foy_c_he_mu`1`_d`2`.m",muind,(d+1)],foy];
		Export[ToString@StringForm["foystats_c_he_mu`1`_d`2`.m",muind,(d+1)],foystats];
		(* Calculate alpha x and y *)
		{alphax,alphay}=Table[
			With[
				{
					mu=Quantity[xory[[1]],"ElectronMass"],
					f0xy=xory[[2]],
					l=If[rightd[d]==0,Quantity[lphos,"BohrRadius"],Quantity[2*lphos,"BohrRadius"]]
				},
				Table[
					Function[{nx,gamma},Evaluate[UnitConvert[(( Pi*(ee^2) )/( 2 * e0 * Sqrt[kappa] * mu * cc )) * ( ( Quantity[nx,"Meters"^-2] ) / ( l ) ) * f0xy[[ni,nf]]*(2/Quantity[gamma,"Seconds"^-1]),"Meters"^-1]]],
					{ni,nmax},{nf,1,nmax-ni}
				]
			],{xory,{{mx,fox},{my,foy}}}
		];
		Export[ToString@StringForm["alphax_c_he_mu`1`_d`2`.m",muind,(d+1)],alphax];
		Export[ToString@StringForm["alphay_c_he_mu`1`_d`2`.m",muind,(d+1)],alphay];
		(* Calculate afac x and y *)
		{afax,afay}=Table[
			With[
				{
					mu=Quantity[xory[[1]],"ElectronMass"],
					f0xy=xory[[2]]
				},
				Table[
					Function[{nx,gamma},Evaluate[1-Exp[-(( Pi*(ee^2) )/( 2 * e0 * Sqrt[kappa] * mu * cc )) * ( ( Quantity[nx,"Meters"^-2] ) ) * f0xy[[ni,nf]] * (2/Quantity[gamma,"Seconds"^-1])]]],
					{ni,nmax},{nf,1,nmax-ni}
				]
			],{xory,{{mx,fox},{my,foy}}}
		];
		Export[ToString@StringForm["afax_c_he_mu`1`_d`2`.m",muind,(d+1)],afax];
		Export[ToString@StringForm["afay_c_he_mu`1`_d`2`.m",muind,(d+1)],afay];
		ToString[StringForm["`3`h`4`m on `5`/`6`: loop complete. Total runtime for this leg is `7`",muind,d,hour,minute,month,day,DateDifference[starttime,gettime,{"Hours","Minutes","Seconds"}]]]>>>mylog.txt;
		ToString[StringForm["-------------------------------------------",muind,d,hour,minute,month,day]]>>>mylog.txt;
	],
	{muind,mindmax},{d,-1,dindmax}
];
Quit[]
