(* Set up conversion factors and general constants *)
e = 1
H2J = 4.35974417 * (10^-18)
T2S = 2.41884326505 * (10^-17)
B2nm = 0.052917721092
H2eV = 27.21138602
na = (5 * 10^15) * (10^-18)
damp = (10^13) * (T2S)
lBN = (1 / B2nm) * 0.333
lphos = (1 / B2nm) * 0.85


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
	tdstrans=Simplify[tds/.{f->(u[ArcTan[#1],ArcTan[#2]]&)}/.{x->(Tan[\[Xi]]),y->(Tan[\[Psi]])},{(Pi/2)>\[Xi]>-(Pi/2),(Pi/2)>\[Psi]>-(Pi/2)}];
	{ev, ef} = NDEigensystem[
			{
				tdstrans,
				DirichletCondition[u[\[Xi],\[Psi]]==0,Abs[\[Xi]]==(Pi/2-eps)||Abs[\[Psi]]==(Pi/2-eps)]
			},
			u[\[Xi],\[Psi]],
			{\[Xi],\[Psi]}\[Element]Rectangle[{-Pi/2+eps,-Pi/2+eps},{Pi/2-eps,Pi/2-eps}],
			nmax,
			Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->eps}}},"Eigensystem"->{"Arnoldi","MaxIterations"->10^7}}
	];
	{
		{{mx,my},{nhbn,d},pot,kappa,chi,eps},
		{
		UnitConvert[Quantity[ev - shift,"Hartrees"],"Millielectronvolts"],
		NormalizeEF[(#/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),s]&/@ef
		}
	}
]

NormalizeEF[EF_,s_] := Module[
  {norm},
  norm = NIntegrate[(EF)^2,{x,-s,s},{y,-s,s}, MinRecursion -> 5,
    MaxRecursion -> 20];
  EF/Sqrt[norm]
]

xysubNInt[efi_,eff_,op_,s_]:=Chop[N[NIntegrate[
	Conjugate[eff]*op*efi,
	{x,-s,s},{y,-s,s},
	MinRecursion->20,MaxRecusion->200,WorkingPrecision->100,AccuracyGoal->6,PrecisionGoal->6
	]^2]]

makeproc[eps_:(10^-3),nmax_:10,mutab_:mus,ntab_:Table[i,{i,0,10}],kaptab_:{1,4.89}]:=Module[
	{
		mun=Dimensions[mutab][[1]],
		nn=Length[ntab],
		kn=Length[kaptab],
		DRare,IRare,DRtime,IRtime,
		s = 1/eps,
		muth=Function[{rare, th},	rare[[ 1, 1, 1 ]]*Cos[th]^2 + rare[[ 1, 1, 2 ]]*Sin[th]^2],
		dpmeth=Function[{rare, i, th}, (rare[[ 3, i ]]*Cos[th]^2 + rare[[ 4, i ]]*Sin[th]^2)^2],
		etrfunc=Function[{rare, ni, nf}, rare[[ 2, nf ]] - rare[[ 2, ni ]] ]
	},
	DRtime=AbsoluteTiming[DRare=Module[
		{
			raw
		},
		Flatten[ParallelTable[Quiet[
			raw=CompPhos2[{mutab[[ mui, 1 ]], mutab[[ mui, 2 ]], -1, pot, kaptab[[ki]]}, eps, nmax];
			{
				raw[[1]],
				raw[[2,1]],
				Table[xysubNInt[ raw[[ 2, 2, 1 ]], raw[[ 2, 2, i ]], x ], {i, nmax}],
				Table[xysubNInt[ raw[[ 2, 2, 1 ]], raw[[ 2, 2, i ]], y ], {i, nmax}],
				Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], 1 ], {i, nmax}, {j, i, nmax}],
				Dimensions[raw]
			}],
			{mui, mun}, {pot, {VKeld, VCoul}}, {ki, 2}
		],{{2},{3},{1}}]
	]][[1]];
	IRtime=AbsoluteTiming[IRare=Module[
		{
			raw
		},
		Flatten[ParallelTable[Quiet[
			raw=CompPhos2[{mutab[[ mui, 1 ]], mutab[[ mui, 2 ]], ntab[[ni]], pot, 4.89}, eps, nmax];
			{
				raw[[1]],
				raw[[2,1]],
				Table[xysubNInt[ raw[[ 2, 2, 1 ]], raw[[ 2, 2, i ]], x ], {i, nmax}],
				Table[xysubNInt[ raw[[ 2, 2, 1 ]], raw[[ 2, 2, i ]], y ], {i, nmax}],
				Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], 1 ], {i, nmax}, {j, i, nmax}],
				Dimensions[raw]
			}],
			{ni, nn}, {mui, mun}, {pot, {VKeld, VCoul}}
		],{{3},{1},{2}}]
	]][[1]];
	Join[
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
			<|Table[ #["kornkey"][[korn]] ->
				<|Table[ ToString@StringForm["mu``",mu] ->
					<|
						"mux" -> #["rare"][korn, mu, 1, 1, 1],
						"muy" -> #["rare"][korn, mu, 1, 1, 2],
						"muth" -> Interpolation[ Table[ muth[ #["rare"][korn, mu], th ], {th, 0, Pi, Pi/20} ] ],
						"evs" -> #["rare"][korn, mu, 2],
						"etr" -> Table[ etrfunc[ #["rare"][korn, mu], i, j ], {i, nmax}, {j, i, nmax}],
						"f0th" -> Table[
								Interpolation[Table[
									2 * muth[ #["rare"][korn, mu], th ] * dpmeth[ #["rare"][korn, mu], j, th ] * etrfunc[ #["rare"][korn, mu], 1, j ],
									{th, 0, Pi, Pi/20}
								]],
							{j,nmax}
						],
						"norm" -> #["rare"][korn, mu, 5]
					|>,
					{mu,mun}
				]|>,
			#["korniter"]
			]|>&,
			<|
				"K" -> <|
					"D" -> <|
						"rare" -> DRare[[1]],
						"korn" -> ToString[k],
						"korniter" -> {korn, kn},
						"kornkey" -> {"FS", "Enc"}
					|>,
					"I" -> <|
						"rare" -> IRare[[1]],
						"korn" -> "Nhbn",
						"korniter" -> {korn, nn},
						"kornkey" -> Table[ToString@StringForm["n``", n], {n, nn}]
					|>
				|>,
				"C" -> <|
					"D" -> <|
						"rare" -> DRare[[2]],
						"korn" -> ToString[k],
						"korniter" -> {korn, kn},
						"kornkey" -> {"FS", "Enc"}
					|>,
					"I" -> <|
						"rare" -> IRare[[2]],
						"korn" -> "Nhbn",
						"korniter" -> {korn, nn},
						"kornkey" -> Table[ToString@StringForm["n``", n], {n, nn}]
					|>
				|>
			|>,
			{2}
		]
	]
] (* function (Module) ends here *)
