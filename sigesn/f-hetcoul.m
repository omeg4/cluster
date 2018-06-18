HetCoul[maxm_, matpars1_, matpars2_, Eperp_, kappa_, dee_, pm1_] := With[
  {
   (* unpack material parameters *)
   ingap1 = matpars1[[1]],
   d01 = matpars1[[2]],
   vF1 = matpars1[[3]],
   ds1 = matpars1[[4]],
   epsrel1 = matpars1[[5]],
   ingap2 = matpars2[[1]],
   d02 = matpars2[[2]],
   vF2 = matpars2[[3]],
   ds2 = matpars2[[4]],
   epsrel2 = matpars2[[5]],
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
    maxcell = 10^-5,
    maxiter = 10^7,
    \[Kappa] = kappa,
    \[Rho] = (QuantityMagnitude[(ds1*Unm/Ubohr)*epsrel1/(2*kappa)]+QuantityMagnitude[(ds2*Unm/Ubohr)*epsrel2/(2*kappa)])/2,
    (* for \[Mu], 
    convert input to SI units and then take the number (should come \
out in kg) and convert to units of Subscript[m, 0] *)
    Egap1 = Abs[
      pm1*(ingap1*Umev/Ujoul*Ujoul) - (d01*Unm/Um)*(Eperp*Uev/Uang*
          Um)],
    Egap2 = Abs[
      pm1*(ingap2*Umev/Ujoul*Ujoul) - (d02*Unm/Um)*(Eperp*Uev/Uang*
          Um)],
	mstar1 = QuantityMagnitude[
      Abs[pm1*(ingap1*Umev/Ujoul*Ujoul) - (d01*Unm/Um)*(Eperp*Uev/Uang*
            Um)]/((vF1*Um/Usec)^2), "ElectronMass"],
	mstar2 = QuantityMagnitude[
      Abs[pm1*(ingap2*Umev/Ujoul*Ujoul) - (d02*Unm/Um)*(Eperp*Uev/Uang*
            Um)]/((vF2*Um/Usec)^2), "ElectronMass"],
    e = 1,
    shift = 10,
    mmax = maxm,
    radialeqs,
    solnmat,
    Ep = Eperp
    },
	\[Mu] = (mstar1 * mstar2)/(mstar1 + mstar2);
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
	{{Egap1,Egap2}, \[Mu], Ep, d},
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i, 1, -1}],
    Table[NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i, 1, -1}]
	}
   ]
]


