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

(* Quick function that returns 0 for direct ("dee" = -1 for direct) but Nhbn*lbn + lphos if "dee" >= 0 *)
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

makeproc[eps_:(10^-3),nmax_:10,mutab_:mus,ntab_:Table[i,{i,0,10}],kaptab_:{1,4.89}]:=Module[
	{
		mun=Dimensions[mutab][[1]],
		nn=Length[ntab],
		kn=Length[kaptab],
		DRare,IRare,
		s = 1/eps,
		muth=Function[{rare, mui, th},	rare[[ 1, 1, 1 ]]*Cos[th]^2 + rare[[ 1, 1, 2 ]]*Sin[th]^2],
		dpmeth=Function[{rare, th}, (rare[[ 3, 1 ]]*Cos[th]^2 + rare[[ 3, 2 ]]*Sin[th]^2)^2],
		etrfunc=Function[{rare, ni, nf}, rare[[ 2, nf ]] - rare[[ 2, ni ]] ]
	},
	(* Need to re-write this whole section, running into RAM issues!
	{
		rawRKD=Flatten[ParallelTable[
			CompPhos2[{mutab[[mui,1]],mutab[[mui,2]],-1,VKeld,kaptab[[ki]]},eps,nmax],
			{mui,mun},{ki,kn}
		],{{2},{1}}],
		rawCD=Flatten[ParallelTable[
			CompPhos2[{mutab[[mui,1]],mutab[[mui,2]],-1,VCoul,kaptab[[ki]]},eps,nmax],
			{mui,mun},{ki,kn}
		],{{2},{1}}],
		rawRKI=ParallelTable[
			CompPhos2[{mutab[[mui,1]],mutab[[mui,2]],ntab[[ni]],VKeld,4.89},eps,nmax],
			{ni,nn},{mui,mun}
		],
		rawCI=ParallelTable[
			CompPhos2[{mutab[[mui,1]],mutab[[mui,2]],ntab[[ni]],VCoul,4.89},eps,nmax],
			{ni,nn},{mui,mun}
		]
	};
	{
		dpmeRKD=ParallelTable[
			xysubNInt[
				rawRKD[[ki,mui,2,2,1]],
				rawRKD[[ki,mui,2,2,i]],
				#
			]&/@{x,y},
			{mui,mun},{ki,kn}
		],
		dpmeCD=ParallelTable[
			xysubNInt[
				rawCD[[ki,mui,2,2,1]],
				rawCD[[ki,mui,2,2,i]],
				#
			]&/@{x,y},
			{mui,mun},{ki,kn}
		],
		dpmeRKI=ParallelTable[
			xysubNInt[
				rawRKI[[ni,mui,2,2,1]],
				rawRKI[[ni,mui,2,2,i]],
				#
			]&/@{x,y},
			{ni,nn},{mui,mun}
		],
		dpmeCI=ParallelTable[
			xysubNInt[
				rawCI[[ni,mui,2,2,1]],
				rawCI[[ni,mui,2,2,i]],
				#
			]&/@{x,y},
			{ni,nn},{mui,mun}
		]
	};
	*)
	DRare=Module[
		{
			raw
		},
		Flatten[ParallelTable[
			raw=CompPhos2[{mutab[[ mui, 1 ]], mutab[[ mui, 2 ]], -1, pot, kaptab[[ki]]}, eps, nmax];
			{
				raw[[1]],
				raw[[2,1]],
				Table[xysubNInt[ raw[[ 2, 2, 1 ]], raw[[ 2, 2, i ]], x ], {i, nmax}],
				Table[xysubNInt[ raw[[ 2, 2, 1 ]], raw[[ 2, 2, i ]], y ], {i, nmax}],
				Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], 1 ], {i, nmax}, {j, i, nmax}]
			},
			{mui, mun}, {pot, {VKeld, VCoul}}, {ki, 2}
		],{{2},{3},{1}}]
	];
	IRare=Module[
		{
			raw
		},
		Flatten[ParallelTable[
			raw=CompPhos2[{mutab[[ mui, 1 ]], mutab[[ mui, 2 ]], ntab[[ni]], pot, 4.89}, eps, nmax];
			{
				raw[[1]],
				raw[[2,1]],
				Table[xysubNInt[ raw[[ 2, 2, 1 ]], raw[[ 2, 2, i ]], x ], {i, nmax}],
				Table[xysubNInt[ raw[[ 2, 2, 1 ]], raw[[ 2, 2, i ]], y ], {i, nmax}],
				Table[xysubNInt[ raw[[ 2, 2, i ]], raw[[ 2, 2, j ]], 1 ], {i, nmax}, {j, i, nmax}]
			},
			{ni, nn}, {mui, mun}, {pot, {VKeld, VCoul}}
		],{{3},{1},{2}}]
	];
	Join[
		Association[
			"lbn"->UnitConvert[Quantity[lBN,"BohrRadius"],"Nanometers"],
			"lphos"->UnitConvert[Quantity[lphos,"BohrRadius"],"Nanometers"],
			"chi2d"->Quantity[0.41,"Nanometers"]
		],
		Map[
			Association[
				"evi"->Function[{korn,mui,i},
					#["rare"][[ korn, mui, 2, i ]]
				],
				"etr"->Function[{korn,mui,i,j},
					etrfunc[ #["rare"][[ korn, mui ]], i, j]
				],
				"f0th"->Function[{korn,mui,j,theta},
					2 * muth[ mui, theta ] * dpmeth[ #["rare"][[ korn, mui ]], theta ] * UnitConvert[etrfunc[ #["rare"][[ korn, mui ]], 1, j],"Hartrees"]
				]
			]&,
			Association[
				"K"->Association[
					"D"->Association[
						"rare"->DRare[[1]]
					],
					"I"->Association[
						"rare"->IRare[[1]]
					]
				],
				"C"->Association[
					"D"->Association[
						"rare"->DRare[[2]]
					],
					"I"->Association[
						"rare"->IRare[[2]],
					]
				]
			], (* outer association with "C"/"K" and "D"/"I" ends here *)
			{2}
		] (* map ends here *)
	] (* join ends here *)
] (* function (Module) ends here *)
