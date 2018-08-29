(* Set up the NDEigensystem function wrapper *)
Coulomb[k_,d_][x_,y_]:= -1 / (k*Sqrt[x^2 + y^2 + d^2])

Keldysh[k_,d_,chi_][x_,y_] := - (Pi / (chi)) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/(chi / (2 * k))] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/(chi / (2 * k))])

CompPhos[nmax_,params_,d_,pot_,reg_,s_,ns_]:=Module[
	{
		mx = params[[1]],
		my = params[[2]],
		chi = params[[3]],
		kappa = params[[4]],
		shift = 10,
		V = If[pot=="K",
			Keldysh[params[[4]],d,params[[3]]],
			Coulomb[params[[4]],d]
		],
		omega = If[reg=="D",
			Disk[{0,0},s],
			Rectangle[{-s,-s},{s,s}]
		],
		ev,
		ef,
		norms
	},
	{ev, ef} = NDEigensystem[
			{
				- (1/2) * ( (1/mx)*D[f[x,y],{x,2},{y,0}] + (1/my)*D[f[x,y],{x,0},{y,2}]) + V[x,y] * f[x,y] + shift*f[x,y],
				DirichletCondition[f[x,y] == 0, True]
			},
			f[x,y],
			{x,y} \[Element] omega,
			nmax,
			Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->s/ns}}},"Eigensystem"->{"Arnoldi",MaxIterations->10^6}}
	];
	{
		{mx,my,chi,kappa,d,s},
		{ev - shift, ef}
	}
]

NormalizeEF[EF_,s_] := Module[
  {norm},
  norm = NIntegrate[(EF)^2,{x,-s,s},{y,-s,s}, MinRecursion -> 5,
    MaxRecursion -> 20];
  EF/Sqrt[norm]
]

(* ProcessPhos[]:=Module[*)
(*     {*)
(*         data = Import["suite.m"],*)
(*         params = data[[1]],*)
(*         dtab*)
(*     },*)
(*     Nd = Dimensions[data][[1]];*)
(*     {params,*)
(*         Table[*)
(*             {*)
(*                 Quantity[1000*H2eV*data*)
(*     }*)
(* ]*)
