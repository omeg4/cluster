(* ::Package:: *)

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



(* ::Package:: *)

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
chiphos = 0.41 / B2nm;
rho = 2 * Pi * chiphos / kappa;
(* ::Package:: *)

(* Set up the NDEigensystem function wrapper *)
VCoul[k_][d_,chi_][x_,y_]:= -1 / (k*Sqrt[x^2 + y^2 + d^2])
VKeld[k_][d_,chi_][x_,y_] := - (1/(4*chi)) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/(2*Pi*chi/k)] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/(2*Pi*chi/k)])
VCho[eps_][d_,chi_][x_,y_]:= (1/(2*eps*(d^3)))*(x^2+y^2) - (1/(eps*d))
VKho[kappa_][d_,chi_][x_,y_]:=-Pi/(4*kappa*((2*Pi*chi/(2*kappa))^2)*d)*(StruveH[-1,d/((2*Pi*chi/(2*kappa)))]-BesselY[-1,d/((2*Pi*chi/(2*kappa)))])
hoCpsi[m_,eps_,d_][n_][x_]:=With[{a=Sqrt[1/Sqrt[(2*m*(1/(2*eps*d^3)))]]},Sqrt[1/(Sqrt[Pi]*a)]*1/Sqrt[(2^n)*Factorial[n]]*Exp[-(x^2)/(2*(a^2))]*HermiteH[n,x/a]]
hoCen[mx_,my_,eps_,d_][n_,m_]:=With[{gamma=1/(2*eps*d^3),v0=(1/(eps*d))},Sqrt[2*gamma/mx]*(n+1/2)+Sqrt[2*gamma/my]*(m+1/2)-v0]

(* Quick function that returns 0 for direct ("dee" = -1 for direct) but Nhbn*lbn + lphos if "dee" >= 0 *)
rightd[d_]:=If[d==-1,0,d*lBN + lphos]

CompPhos[nmax_,params_,dee_,pot_,eps_]:=Module[
	{
		mx = params[[1]],
		my = params[[2]],
		chi = params[[3]],
		kappa = params[[4]],
		shift = 10,
		d = rightd[dee],
		ev,
		ef,
		norms,
		tds,
		tdstrans
	},
	tds=- (1/2) * ( (1/mx)*D[f[x,y],{x,2},{y,0}] + (1/my)*D[f[x,y],{x,0},{y,2}]) + pot[d,chi][x,y] * f[x,y] + shift*f[x,y];
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
		{mx,my,chi,kappa,d,eps},
		{UnitConvert[Quantity[ev - shift,"Hartrees"],"Millielectronvolts"], ef}
	}
]

NormalizeEF[EF_,s_] := Module[
  {norm},
  norm = NIntegrate[(EF)^2,{x,-s,s},{y,-s,s}, MinRecursion -> 5,
    MaxRecursion -> 20];
  EF/Sqrt[norm]
]

ProcessPhos[data_]:=Module[
	{
		params = data[[1]],
		mx=data[[1]][[1]],
		my=data[[1]][[2]],
		evs = data[[2]][[1]],
		efs = data[[2]][[2]],
		s=1/data[[1]][[6]],
		norms,
		efn,
		normtable,
		xtable,
		ytable,
		xyavgtable,
		xyredtable
	},
	(* normalize eigenfunctions first *)
	norms=Table[
			NIntegrate[
				(efs[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})^2,
				{x,-s,s},{y,-s,s},
				MinRecursion->20,MaxRecursion->200,WorkingPrecision->100
			],
		{i,Length[efs]}
	];
	efn=Table[efs[[i]]/Sqrt[norms[[i]]],{i,Length[efs]}];
	normtable=Table[
		N[NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100],5]//Chop,
		{i,10},{j,i,10}
		];
	xtable=Table[
		N[2*mx*Abs[QuantityMagnitude@UnitConvert[evs[[i]]-evs[[j]],"Hartrees"]]*(NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*x*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100])^2,5]//Chop,
		{i,10},{j,i,10}
		];
	ytable=Table[
		N[2*my*Abs[QuantityMagnitude@UnitConvert[evs[[i]]-evs[[j]],"Hartrees"]]*(NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*y*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100])^2,5]//Chop,
		{i,10},{j,i,10}
		];
	xyavgtable=Table[
		N[2*((mx+my)/2)*Abs[QuantityMagnitude@UnitConvert[evs[[i]]-evs[[j]],"Hartrees"]]*(NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*Sqrt[x^2+y^2]*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100])^2,5]//Chop,
		{i,10},{j,i,10}
		];
	xyredtable=Table[
		N[2*((mx*my)/(mx+my))*Abs[QuantityMagnitude@UnitConvert[evs[[i]]-evs[[j]],"Hartrees"]]*(NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*Sqrt[x^2+y^2]*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100])^2,5]//Chop,
		{i,10},{j,i,10}
		];
		{params,evs,normtable,xtable,ytable,xyavgtable,xyredtable}
]

ProcessPhosInd[data_]:=Table[
	Module[
	{
		params = data[[i]] [[d]] [[1]],
		mx=data[[i]] [[d]] [[1]] [[1]],
		my=data[[i]] [[d]] [[1]] [[2]],
		evs = data[[i]] [[d]] [[2]] [[1]],
		efs = data[[i]] [[d]] [[2]] [[2]],
		s=1/data[[i]] [[d]] [[1]] [[6]],
		norms,
		efn,
		normtable,
		xtable,
		ytable,
		xyavgtable,
		xyredtable
	},
	(* normalize eigenfunctions first *)
	norms=Table[
			NIntegrate[
				(efs[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})^2,
				{x,-s,s},{y,-s,s},
				MinRecursion->20,MaxRecursion->200,WorkingPrecision->100
			],
		{i,Length[efs]}
	];
	efn=Table[efs[[i]]/Sqrt[norms[[i]]],{i,Length[efs]}];
	normtable=Table[
		N[NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100],5]//Chop,
		{i,10},{j,i,10}
		];
	(* xtable=Table[*)
	(*     N[2*mx*Abs[QuantityMagnitude@UnitConvert[evs[[i]]-evs[[j]],"Hartrees"]]*(NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*x*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100])^2,5]//Chop,*)
	(*     {i,10},{j,i,10}*)
	(*     ];*)
	(* ytable=Table[*)
	(*     N[2*my*Abs[QuantityMagnitude@UnitConvert[evs[[i]]-evs[[j]],"Hartrees"]]*(NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*y*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100])^2,5]//Chop,*)
	(*     {i,10},{j,i,10}*)
	(*     ];*)
	(* xyavgtable=Table[*)
	(*     N[2*((mx+my)/2)*Abs[QuantityMagnitude@UnitConvert[evs[[i]]-evs[[j]],"Hartrees"]]*(NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*Sqrt[x^2+y^2]*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100])^2,5]//Chop,*)
	(*     {i,10},{j,i,10}*)
	(*     ];*)
	(* xyredtable=Table[*)
	(*     N[2*((mx*my)/(mx+my))*Abs[QuantityMagnitude@UnitConvert[evs[[i]]-evs[[j]],"Hartrees"]]*(NIntegrate[(efn[[i]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]})*Sqrt[x^2+y^2]*(Conjugate[efn[[j]]]/.{\[Xi]->ArcTan[x],\[Psi]->ArcTan[y]}),{x,-s,s},{y,-s,s}, MinRecursion->20,MaxRecursion->200,WorkingPrecision->100])^2,5]//Chop,*)
	(*     {i,10},{j,i,10}*)
	(*     ];*)
		{params,evs,normtable(*,xtable,ytable,xyavgtable,xyredtable*)}
	],
	{i,Dimensions[data][[1]]},{d,Dimensions[data][[2]]}
]

FPparams[proc_]:=proc[[1]]
FPev[proc_,n_]:=proc[[2]][[n]]
FPnorms[proc_,ni_,nf_]:=proc[[3]][[ni]][[nf]]
FPf0x[proc_,ni_,nf_]:=proc[[4]][[ni]][[nf]]
FPf0y[proc_,ni_,nf_]:=proc[[5]][[ni]][[nf]]
FPf0xya[proc_,ni_,nf_]:=proc[[6]][[ni]][[nf]]
FPf0xyr[proc_,ni_,nf_]:=proc[[7]][[ni]][[nf]]
SetDirectory["/home/mbrunetti/cluster/phos/results/VicCxact-mr-2018-11-29"]
kappa = 4.89;
Export["inps.txt",{VicCxact-mr,mus,chiphos,kappa,{50,50},VCoul[4.89],10^-4}];
result=Table[
	CompPhos[10,{mus[[i]][[1]],mus[[i]][[2]],chiphos,kappa},nhbn,VCoul[4.89],10^-4],
	{i,1},{nhbn,50,50}
	];
Export["suite.m",result];
Export["proc.m",ProcessPhosInd[result]];
Quit[]
