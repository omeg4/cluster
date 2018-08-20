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
