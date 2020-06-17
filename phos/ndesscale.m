(* Define the potentials *)
VCoul[k_][d_,chi_][x_,y_]:= -1 / (k*Sqrt[x^2 + y^2 + d^2])
VKeld[k_][d_,chi_][x_,y_] := - (1/(4*chi)) * (StruveH[0,Sqrt[x^2 + y^2 + d^2]/(2*Pi*chi/k)] - BesselY[0,Sqrt[x^2 + y^2 + d^2]/(2*Pi*chi/k)])

(*
function ndesscale solves the schrodinger equation describing the electron and hole, returning the eigenfunctions and eigenenergies of the exciton
	inputs:
		nmax: number of eigenstates to calculate
		mx, my: anisotropic exciton effective masses
		pot: a function Head which describes the electron-hole potential. note that only the function head (e.g. VKeld) is passed as an input argument to ndesscale.
		kappa: (\[epsilon]_1 + \[epsilon]_2)/2
		chi: the 2d polarizability
		d: vertical separation between electron and hole. d=0 for direct excitons, d>0 for indirect.
		eps: maximum distance between two finite mesh elements.
		c: transformation axes scaling factor. c \approx 10 seemed to work very well for the level of detail i needed. it would be interesting to explore optimizing this number more.
		maxi: max number of iterations for converging the eigenfunctions and eigenvalues. i found diminishing returns past about 10^12.
		vn: just keep it None.
	outputs:
		eigenvalues (List of numbers): starting with the ground state energy.
		eigenfunctions (List of InterpolatingFunction's): again starting with the ground state energy. the List elements are implicit functions of [x,y].
*)
ndesscale[nmax_, mx_, my_, pot_, kappa_, chi_, d_, eps_, c_, maxi_, vn_]:=Module[
	{
		tds,
		tdstrans,
		shift=10,
		stranseps=N[(Pi/2)-$MachineEpsilon],
		evs,efs
	},
	tds = (-(1/2)*((1/mx)*D[f[x, y], {x, 2}, {y, 0}] + (1/my)* D[f[x, y], {x, 0}, {y, 2}]) + pot[kappa][d, chi][x, y]*f[x, y] + shift*f[x, y]);
	tdstrans = Simplify[tds /. {f -> (u[ArcTan[#1/c], ArcTan[#2/c]] &)} /. {x -> (c*Tan[\[Xi]]), y -> (c*Tan[\[Psi]])}, {(Pi/2) > \[Xi] > -(Pi/2), (Pi/2) > \[Psi] > -(Pi/2)}];
	{evs,efs}=NDEigensystem[
		{
			tdstrans,
			DirichletCondition[u[\[Xi], \[Psi]] == 0,Abs[\[Xi]] == (stranseps) || Abs[\[Psi]] == (stranseps)]
		},
	 u[\[Xi], \[Psi]],
	 {\[Xi], \[Psi]} \[Element] Rectangle[{-stranseps, -stranseps}, {stranseps, stranseps}],
	 nmax,
	 Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->eps}}},"Eigensystem"->{"Arnoldi","MaxIterations"->maxi},"VectorNormalization"->vn}
	];
	{
		evs-shift,
		Table[
			Function[{x,y},
				Evaluate[
					Head[efs[[i]]][ArcTan[x/c],ArcTan[y/c]]
				]
			],
			{i,nmax}
		]
	}
]
