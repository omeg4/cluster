(* ::Package:: *)

(* Set up the NDEigensystem function wrapper *)
Coulomb[k_,d_][x_,y_]:= -1 / (k*Sqrt[x^2 + y^2 + d^2])

Keldysh[k_,d_,chi_][x_,y_] := - (Pi / (chi)) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/(chi / (2 * k))] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/(chi / (2 * k))])

Harmosc[mx_,my_,eps_,d_,v0_][x_,y_]:= (1/(2*eps*(d^3)))*(x^2+y^2) - (1/(eps*d))

psi1d[a_][n_][x_]:=Sqrt[1/(Sqrt[Pi]*a)]*1/Sqrt[(2^n)*Factorial[n]]*Exp[-(x^2)/(2*a)]*HermiteH[n,x/a]

V[pot_] := If[pot=="K",
	Keldysh[kappa,d,chi],
	Coulomb[kappa,d]
]

omega[reg_] := If[reg=="D",
	Disk[{0,0},s],
	Rectangle[{-s,-s},{s,s}]
]

CompPhos[nmax_,params_,dee_,pot_]:=Module[
	{
		mx = params[[1]],
		my = params[[2]],
		chi = params[[3]],
		kappa = params[[4]],
		shift = 10,
		d = dee,
		ev,
		ef,
		norms,
		tds,
		tdstrans,
		eps=10^-3,
		MC=10^-4
	},
	tds=- (1/2) * ( (1/mx)*D[f[x,y],{x,2},{y,0}] + (1/my)*D[f[x,y],{x,0},{y,2}]) + pot[x,y] * f[x,y] + shift*f[x,y];
	tdstrans=Simplify[tds/.{f->(u[ArcTan[#1],ArcTan[#2]]&)}/.{x->(Tan[\[Xi]]),y->(Tan[\[Psi]])},{(Pi/2)>\[Xi]>-(Pi/2),(Pi/2)>\[Psi]>-(Pi/2)}];
	{ev, ef} = NDEigensystem[
			{
				tdstrans,
				DirichletCondition[u[\[Xi],\[Psi]]==0,Abs[\[Xi]]==(Pi/2-eps)||Abs[\[Psi]]==(Pi/2-eps)]
			},
			u[\[Xi],\[Psi]],
			{\[Xi],\[Psi]}\[Element]Rectangle[{-Pi/2+eps,-Pi/2+eps},{Pi/2-eps,Pi/2-eps}],
			nmax,
			Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->MC}}},"Eigensystem"->{"Arnoldi",MaxIterations->10^5}}
	];
	{
		{mx,my,chi,kappa,d},
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
		evs = data[[2]][[1]],
		efs = data[[2]][[2]]
	},
	{params,
		Table[
			{
				evs[[i]],
				Table[
					{
						N[NIntegrate[efs[[i]]*Conjugate[efs[[j]]],{x,-s,s},{y,-s,s}, MinRecursion->10,MaxRecursion->100,WorkingPrecision->100],5],

						N[NIntegrate[efs[[i]]*x*Conjugate[efs[[j]]],{x,-s,s},{y,-s,s}, MinRecursion->10,MaxRecursion->100,WorkingPrecision->100],5],
						N[NIntegrate[efs[[i]]*y*Conjugate[efs[[j]]],{x,-s,s},{y,-s,s}, MinRecursion->10,MaxRecursion->100,WorkingPrecision->100],5]
					},
					{j,Length[evs]}
				]
			},
			{i,Length[evs]}
		]
	}
]
