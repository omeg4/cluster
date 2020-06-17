(* Set up conversion factors and general constants *)
kappa = 4.89;
e = 1;
H2J = 4.35974417 * (10^-18);
T2S = 2.41884326505 * (10^-17);
B2nm = 0.052917721092;
H2eV = 27.21138602;
na = (5 * 10^15) * (10^-18);
damp = (10^13) * (T2S);
lBN = (1 / B2nm) * 0.333;
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
CompPhosKeld[mux_,muy_,kap_,rh_,dee_,nmax_,s_,ns_]:=Module[
	{
		kappa = kap,
		rho = rh,
		d = dee,
		mx = mux,
		my = muy,
		maxcell = s/ns,
		shift = 10,
		omega = 10,
		ev,
		ef,
		norms
	},
	{ev, ef} = NDEigensystem[
			{
				- (1/2) * ( (1/mx)*D[f[x,y],{x,2},{y,0}] + (1/my)*D[f[x,y],{x,0},{y,2}]) + ( (Pi / ( 2 * kappa * rho )) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/rho] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/rho] )) * f[x,y] + shift*f[x,y],
				DirichletCondition[f[x,y] == 0, True]
			},
			f[x,y],
			{x,y} \[Element] Rectangle[{-s,-s},{s,s}],
			nmax,
			Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->maxcell}}},"Eigensystem"->{"Arnoldi",MaxIterations->10^3}}
	];
	{ev - shift, ef}
]
s = 500;
ns = 5000;
result=CompPhosKeld[mus[[1]][[1]],mus[[1]][[2]],4.89,rho,40*lBN,10,s,ns]//AbsoluteTiming;
Save[ToString@StringForm['/home/mbrunetti/phos/`1`-`2`_keld_s`3`_ns`4`,DateObject[][[1]][[2]],DateObject[][[1]][[3]],s,ns],result];
Export[ToString@StringForm['/home/mbrunetti/phos/`1`-`2`_keld_s`3`_ns`4`,DateObject[][[1]][[2]],DateObject[][[1]][[3]],s,ns],result];
Quit[]
