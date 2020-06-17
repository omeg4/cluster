(* ::Package:: *)

MBDirKeld[maxm_, matpars_, kappa_,guess_] := With[
  {
   (* unpack material parameters *)
   a = UnitConvert[Quantity[matpars[[1]],"Angstroms"],"BohrRadius"]//QuantityMagnitude,
   t = UnitConvert[Quantity[matpars[[2]],"Electronvolts"],"Hartrees"]//QuantityMagnitude,
   delta = UnitConvert[Quantity[matpars[[3]],"Electronvolts"],"Hartrees"]//QuantityMagnitude,
   lambda = UnitConvert[Quantity[matpars[[4]],"Electronvolts"],"Hartrees"]//QuantityMagnitude,
   chi2d = UnitConvert[Quantity[matpars[[5]],"Angstroms"],"BohrRadius"]//QuantityMagnitude,
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
    e = 1,
    shift = 10,
    mmax = maxm,
    radialeqs,
    solnmat,
    dstau,
    \[Mu]
    },
    dstau=delta - (lambda/2);
    \[Mu]=dstau/(2*(a^2)*(t^2));
    Vr=(((((Pi*(e^2))/(2*\[Kappa]*\[Rho]))*(StruveH[0,r/\[Rho]] - BesselY[0, r/\[Rho]]))) - (m^2/(2*\[Mu]* r^2)));
    radialEqKeld = -((2*(a^2)*(t^2))/(guess+Vr))*f''[r]-D[((2*(a^2)*(t^2))/(guess+Vr)),r]* f'[r] +(dstau-Vr)* f[r];
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
	{dstau, \[Mu], 0, 0},
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i, 1, -1}],
	Table[efTab[[i - j + 1]][[j]], {i, Dimensions[efTab][[1]]}, {j, i, 1, -1}]
	}
	]
]


ConvMB[maxm_,matpars_,kappa_,guess_,tol_]:=Module[
	{
	EB = guess,
	LEB = guess+100,
	EBlist={guess+100,guess},
	result,
	LGSE,GSE=guess,dst
	},
	While[Abs[LEB-EB] > tol,
	LGSE = GSE;
	LEB=EB;
	result=MBDirKeld[maxm,matpars,kappa,GSE];
	dst=result[[1]][[1]];
	GSE=result[[2]][[1]][[1]];
	EB=H2meV[dst-Abs@GSE];
	EBlist=Append[EBlist,EB];
	Print[{H2meV[dst],H2meV[GSE],EB}]
	]
	]


H2meV[in_]:=UnitConvert[Quantity[in,"Hartrees"],"Millielectronvolts"]//QuantityMagnitude
meV2H[in_]:=UnitConvert[Quantity[in,"Millielectronvolts"],"Hartrees"]//QuantityMagnitude
