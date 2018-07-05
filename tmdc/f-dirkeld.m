DirKeld[maxm_, matpars_, kappa_] := With[
  {
   (* unpack material parameters *)
   me = matpars[[1]],
   mh = matpars[[2]],
   chi2d = UnitConvert[Quantity[matpars[[3]],"Angstroms"],"BohrRadius"]//QuantityMagnitude,
   \[Kappa] = kappa,
   (* Define all the units *)
   Uhart = Quantity["Hartrees"],
   Ujoul = Quantity["Joules"],
   Uev = Quantity["Electronvolts"],
   Umev = Quantity["Millielectronvolts"],
   Um = Quantity["Meters"],
   Unm = Quantity["Nanometers"],
   Uang = Quantity["Angstroms"],
   Ubohr = Quantity["BohrRadius"],
   Uvolt = Quantity["Volts"],
   Usec = Quantity["Seconds"]
   },
  Module[
   {
    maxcell = 10^-4,
    maxiter = 10^7,
	\[Rho] = 2 * Pi * chi2d / \[Kappa],
	\[Mu] = (me * mh) / (me + mh),
    e = 1,
    shift = 10,
    mmax = maxm,
    radialeqs,
    solnmat
    },
   radialEqKeld = -(1/(\[Mu]*2)) f''[r] - (1/(2*\[Mu]*r))* f'[r] - (((((Pi*(e^2))/(2*\[Kappa]*\[Rho]))*(StruveH[0,r/\[Rho]] - BesselY[0, r/\[Rho]]))) - (m^2/(2*\[Mu]* r^2)))* f[r];
   radial\[Xi]Keld[m_] = Simplify[radialEqKeld /. f -> (\[Psi][ArcTan[#]] &) /. r -> (Tan[\[Xi]]), Pi/2 > \[Xi] > 0];
   solnmat = {};
   evTab = {};
   efTab = {};
   bigarray = {};
   evout = {};
   efout = {};
   Do[
    ev = {};
    {ev, ef} = 
     NDEigensystem[{radial\[Xi]Keld[mind] + shift \[Psi][\[Xi]], 
       DirichletCondition[\[Psi][\[Xi]] == 0, \[Xi] == 
         Pi/2]}, \[Psi][\[Xi]], {\[Xi], 0, Pi/2}, mmax - mind + 1, 
      Method -> {"SpatialDiscretization" -> {"FiniteElement", \
{"MeshOptions" -> {"MaxCellMeasure" -> maxcell}}}, 
        "Eigensystem" -> {"Arnoldi", MaxIterations -> maxiter}}];
    evTab = Append[evTab, ev - shift];
    efTab = Append[efTab, ef /. \[Xi] -> ArcTan[r]],
    {mind, 0, mmax}
    ];
   {
	{Egap, \[Mu], Ep, 0},
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i, 1, -1}],
	Table[NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i, 1, -1}]
	}
	]
]
