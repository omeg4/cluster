(* ::Package:: *)

(* Set up the NDEigensystem function wrapper *)
Vcoul[k_][d_][x_,y_]:= -1 / (k*Sqrt[x^2 + y^2 + d^2])
Vkeld[k_][d_,chi_][x_,y_] := - (1/(4*chi)) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/(2*Pi*chi/k)] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/(2*Pi*chi/k)])
VCho[eps_][d_][x_,y_]:= (1/(2*eps*(d^3)))*(x^2+y^2) - (1/(eps*d))
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
	{i,4},{d,Length[data [[1]]]}
]

FPparams[proc_]:=proc[[1]]
FPev[proc_,n_]:=proc[[2]][[n]]
FPnorms[proc_,ni_,nf_]:=proc[[3]][[ni]][[nf]]
FPf0x[proc_,ni_,nf_]:=proc[[4]][[ni]][[nf]]
FPf0y[proc_,ni_,nf_]:=proc[[5]][[ni]][[nf]]
FPf0xya[proc_,ni_,nf_]:=proc[[6]][[ni]][[nf]]
FPf0xyr[proc_,ni_,nf_]:=proc[[7]][[ni]][[nf]]
