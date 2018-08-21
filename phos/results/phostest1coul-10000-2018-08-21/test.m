(* Set up conversion factors and general constants *)
kappa = 4.89
e = 1
H2J = 4.35974417 * (10^-18)
T2S = 2.41884326505 * (10^-17)
B2nm = 0.052917721092
H2eV = 27.21138602
na = (5 * 10^15) * (10^-18)
damp = (10^13) * (T2S)
lBN = (1 / B2nm) * 0.333
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
(* Set up the NDEigensystem function wrapper *)
Coulomb[k_,d_][x_,y_]:= -1 / (k*Sqrt[x^2 + y^2 + d^2])

Keldysh[k_,d_,chi_][x_,y_] := - (Pi / (chi)) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/(chi / (2 * k))] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/(chi / (2 * k))])

CompPhos[nmax_,params_,d_,pot_,s_]:=Module[
	{
		mx = params[[1]],
		my = params[[2]],
		chi = params[[3]],
		kappa = params[[4]],
		shift = 10,
		omega = 10,
		V = If[pot=="K",
			Keldysh[params[[4]],d,params[[3]]],
			Coulomb[params[[4]],d]
		],
		ev,
		ef,
		norms
	},
	{ev, ef} = NDEigensystem[
			{
				- (1/2) * ( (1/mx)*D[f[x,y],{x,2},{y,0}] + (1/my)*D[f[x,y],{x,0},{y,2}]) - V[x,y] * f[x,y] + shift*f[x,y],
				DirichletCondition[f[x,y] == 0, True]
			},
			f[x,y],
			{x,y} \[Element] Rectangle[{-s,-s},{s,s}],
			nmax,
			Method->{"SpatialDiscretization"->{"FiniteElement"},"Eigensystem"->{"Arnoldi",MaxIterations->10^5}}
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
SetDirectory["/home/mbrunetti/cluster/phos/results/phostest1coul-10000-2018-08-21"]
Export["inps.txt",{phostest1coul,"C",1,1,1,1,0,10000}];
result=CompPhos[4,{1,1,1,1},0,"C",10000];
Export["suite.m",result]
Quit[]
