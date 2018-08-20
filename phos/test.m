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
CompPhosKeld[mux_,muy_,kap_,rh_,dee_,nmax_,s_,ns_]:=Module[
	{
		kappa = kap,
		rho = rh,(* Don't I need (chi_2D * 2pi / kappa) ??? *)
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
				- (1/2) * ( (1/mx)*D[f[x,y],{x,2},{y,0}] + (1/my)*D[f[x,y],{x,0},{y,2}]) - ( (Pi / ( 2 * kappa * rho )) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/rho] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/rho] )) * f[x,y] + shift*f[x,y],
				DirichletCondition[f[x,y] == 0, True]
			},
			f[x,y],
			{x,y} \[Element] Rectangle[{-s,-s},{s,s}],
			nmax,
			Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->maxcell}}},"Eigensystem"->{"Arnoldi",MaxIterations->10^5}}
	];
	{ev - shift, ef}
]
s = 10000
ns = 10000
projdir = ToString[2018-04-26_sns10K_2]
filename = ToString[sns10K_2]
result=Table[
	{
		mux,
		muy,
		Table[{nhbn,CompPhosKeld[mux,muy,4.89,rho,nhbn*lBN,10,s,ns]},{nhbn,10}]
	},
	{mux,{0.062963,0.0910365}},{muy,{0.659901,0.967742}}
]//AbsoluteTiming;
Export[
	"/home/mbrunetti/phos/results/2018-04-26_sns10K_2/2018-04-26_sns10K_2.m",
	result
];
getEV[mux_,muy_,d_,n_]:=result[[mux]][[muy]][[3]][[d]][[2]][[1]][[n]]
getEF[mux_,muy_,d_,n_]:=result[[mux]][[muy]][[3]][[d]][[2]][[2]][[n]]
getmux[mux_]:=result[[mux]][[1]][[1]]
getmuy[muy_]:=result[[1]][[muy]][[2]]
pickmu[mux_,muy_,xy_]:=result[[mux]][[muy]][[xy]]
Lphos=1.117/B2nm
getf0[mux_,muy_,d_,n_,xy_]:=2*pickmu[mux,muy,xy]*(getEV[mux,muy,d,n]-getEV[mux,muy,d,1])*NIntegrate[getEF[mux,muy,d,n]*If[xy==1,x,y]*getEF[mux,muy,d,1],{x,-s,s},{y,-s,s},MaxRecursion->100,WorkingPrecision->100]
(* How to organize processed data? *)
phosprocess={
	chiphos*B2nm*Quantity["Nanometers"],
	Table[
	{
		result[[mux]][[muy]][[1]]*Quantity["ElectronMass"],
		result[[mux]][[muy]][[2]]*Quantity["ElectronMass"],
		{
			"orthonormcheck",
			Table[
				NIntegrate[getEF[mux,muy,d,ni]*getEF[mux,muy,d,nj],{x,-s,s},{y,-s,s},MaxRecursion->100,WorkingPrecision->100],
				{ni,10},{nj,10}
			]	
		},
		{
			"f0 matrix with Etr",
			Table[
				{
					(getEV[mux,muy,d,n]-getEV[mux,muy,d,1])*H2eV*1000,
					getf0[mux,muy,d,n,1],
					getf0[mux,muy,d,n,2]
				},
				{n,2,10}
			]
		}
	},
	{mux,2},{muy,2},{d,10}
	]
};
Export[ToString@StringForm["/home/mbrunetti/phos/results/`1`/processed.m",projdir],phosprocess]
Quit[]
