IndCoul[maxm_, matpars_, ingap_, Eperp_, kappa_, dee_, pm1_] := With[
  {
   (* unpack material parameters *)
   d0 = matpars[[1]],
   vF = matpars[[2]],
   ds = matpars[[3]],
   epsrel = matpars[[4]] - 1,
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
    maxiter = 4*10^5,
    \[Kappa] = kappa,
    \[Rho] = QuantityMagnitude[(ds*Unm/Ubohr)*epsrel/(2*kappa)],
    (* for \[Mu], 
    convert input to SI units and then take the number (should come \
out in kg) and convert to units of Subscript[m, 0] *)
    Egap = Abs[
      pm1*(ingap*Umev/Ujoul*Ujoul) - (d0*Unm/Um)*(Eperp*Uev/Uang*
          Um)],
    \[Mu] = 
     QuantityMagnitude[
      Abs[pm1*(ingap*Umev/Ujoul*Ujoul) - (d0*Unm/Um)*(Eperp*Uev/Uang*
            Um)]/(2*(vF*Um/Usec)^2), "ElectronMass"],
    e = 1,
    shift = 10,
    mmax = maxm,
    radialeqs,
    solnmat,
    Ep = Eperp
    },
   radialEqKInd = -(1/(\[Mu]*2))*f''[r] - (1/(2 \[Mu]*r))* f'[r] - ((1/(kappa*Sqrt[r^2 + d^2])) - (m^2/(2*\[Mu]* r^2)))* f[r];
   radial\[Xi]KInd[m_] = Simplify[radialEqKInd /. f -> (\[Psi][ArcTan[#]] &) /. r -> (Tan[\[Xi]]), Pi/2 > \[Xi] > 0];
   solnmat = {}; evTab = {}; efTab = {};
   bigarray = {};
   Do[
    ev = {};
	{ev, ef} =
     NDEigensystem[{radial\[Xi]KInd[mind] + shift \[Psi][\[Xi]],
       DirichletCondition[\[Psi][\[Xi]] == 0, \[Xi] ==
         Pi/2]}, \[Psi][\[Xi]], {\[Xi], 0, Pi/2}, mmax - mind + 1,
      Method -> {"SpatialDiscretization" -> {"FiniteElement",{"MeshOptions" -> {"MaxCellMeasure" -> maxcell}}},
        "Eigensystem" -> {"Arnoldi", MaxIterations -> maxiter}}];
    evTab = Append[evTab, ev - shift]; efTab = Append[efTab, ef/.{\[Xi]->ArcTan[r]}],
    {mind, 0, mmax}
    ];
	{
	{Egap, \[Mu], Ep, d},
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i, 1, -1}],
    Table[NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i, 1, -1}]
	}
   ]
]

