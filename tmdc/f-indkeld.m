IndKeld[maxm_, matpars_, kappa_, dee_] := With[
  {
   (* unpack material parameters *)
   me = matpars[[1]],
   mh = matpars[[2]],
   chi2d = UnitConvert[Quantity[matpars[[3]],"Angstroms"],"BohrRadius"]//QuantityMagnitude,
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
    d = dee,
    maxcell = 10^-4,
    maxiter = 10^7,
    \[Kappa] = kappa,
    \[Rho] = 2 * Pi * chi2d / kappa,
    \[Mu] = (me * mh) / (me + mh),
	e = 1,
    shift = 10,
    mmax = maxm,
    radialeqs,
    solnmat
    },
   radialEqKInd = -(1/(\[Mu]*2)) f''[r] - (1/(2 \[Mu]*r))* f'[r] - (((Pi/(2*\[Kappa]*\[Rho])*(StruveH[0, Sqrt[r^2 + d^2]/\[Rho]] - BesselY[0, Sqrt[r^2 + d^2]/\[Rho]]))) - (m^2/(2*\[Mu]* r^2)))* f[r];
   radial\[Xi]KInd[m_] = Simplify[radialEqKInd /. f -> (\[Psi][ArcTan[#]] &) /. r -> (Tan[\[Xi]]), Pi/2 > \[Xi] > 0];
   solnmat = {}; evTab = {}; efTab = {};
   bigarray = {};
   Do[
    ev = {};
	{ev, ef} =
     NDEigensystem[{radial\[Xi]KInd[mind] + shift \[Psi][\[Xi]], 
       DirichletCondition[\[Psi][\[Xi]] == 0, \[Xi] == 
         Pi/2]}, \[Psi][\[Xi]], {\[Xi], 0, Pi/2}, mmax - mind + 1, 
      Method -> {"SpatialDiscretization" -> {"FiniteElement", \
{"MeshOptions" -> {"MaxCellMeasure" -> maxcell}}}, 
        "Eigensystem" -> {"Arnoldi", MaxIterations -> maxiter}}]; 
    evTab = Append[evTab, ev - shift]; efTab = Append[efTab/.{\[Xi]->ArcTan[r]}, ef],
    {mind, 0, mmax}
    ];
	{
	{Egap, \[Mu], Ep, d},
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i, 1, -1}],
    Table[NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i, 1, -1}]
	}
   ]
]
