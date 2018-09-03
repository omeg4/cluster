(* Set up the NDEigensystem function wrapper *)
Coulomb[k_,d_][x_,y_]:= -1 / (k*Sqrt[x^2 + y^2 + d^2])

Keldysh[k_,d_,chi_][x_,y_] := - (Pi / (chi)) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/(chi / (2 * k))] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/(chi / (2 * k))])

V[pot_] := If[pot=="K",
	Keldysh[kappa,d,chi],
	Coulomb[kappa,d]
]

omega[reg_] := If[reg=="D",
	Disk[{0,0},s],
	Rectangle[{-s,-s},{s,s}]
]

CompPhos[nmax_,params_,dee_,pot_,reg_,s_,ns_]:=Module[
	{
		mx = params[[1]],
		my = params[[2]],
		chi = params[[3]],
		kappa = params[[4]],
		shift = 10,
		d = dee,
		ev,
		ef,
		norms
	},
	{ev, ef} = NDEigensystem[
			{
				- (1/2) * ( (1/mx)*D[f[x,y],{x,2},{y,0}] + (1/my)*D[f[x,y],{x,0},{y,2}]) + pot[x,y] * f[x,y] + shift*f[x,y],
				DirichletCondition[f[x,y] == 0, True]
			},
			f[x,y],
			{x,y} \[Element] reg,
			nmax,
			Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->s/ns}}},"Eigensystem"->{"Arnoldi",MaxIterations->10^8}}
	];
	{
		{mx,my,chi,kappa,d,s,ns},
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
						N[NIntegrate[efs[[i]]*Conjugate[efs[[j]]],{x,-s,s},{y,-s,s}, MinRecursion->10,MaxRecursion->100,WorkingPrecision->100],10],

						N[NIntegrate[efs[[i]]*x*Conjugate[efs[[j]]],{x,-s,s},{y,-s,s}, MinRecursion->10,MaxRecursion->100,WorkingPrecision->100],10],
						N[NIntegrate[efs[[i]]*y*Conjugate[efs[[j]]],{x,-s,s},{y,-s,s}, MinRecursion->10,MaxRecursion->100,WorkingPrecision->100],10]
					},
					{j,Length[evs]}
				]
			},
			{i,Length[evs]}
		]
	}
]
