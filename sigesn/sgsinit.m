(* Generic values for absorption: na, \[CapitalGamma], etc. *)
B2nm = 0.052917721092;
H2J = 4.35974417*(10^-18);
H2eV = 27.21138602;
T2S = 2.41884326505*(10^-17);
lBN = (1/B2nm)*0.333;
na = (5*10^15)*(10^-18)*(B2nm^2) ;
damp = (10^13)*(T2S);
damp2 = (10^12)*T2S;
damp3 = (10^11)*(T2S);
l = (3.1 * 10^-10)*(10^9)/(B2nm) ;
sipar2 = {1.9, 0.046, 6.5*10^5, 0.4, 11.9};
gepar2 = {33, 0.0676, 6.2*10^5, 0.45, 16};
snpar2 = {101, 0.085, 5.5*10^5, 0.5, 24};
labels = {"B Exciton (Small Gap)", "A Exciton (Large Gap)"};
Uhart = Quantity["Hartrees"];
Ujoul = Quantity["Joules"];
Uev = Quantity["Electronvolts"];
Umev = Quantity["Millielectronvolts"];
Um = Quantity["Meters"];
Unm = Quantity["Nanometers"];
Uang = Quantity["Angstroms"];
Ubohr = Quantity["BohrRadius"];
Uvolt = Quantity["Volts"];
Usec = Quantity["Seconds"];
Uemass = Quantity["ElectronMass"];
DirKeld[maxm_, matpars_, ingap_, Eperp_, kappa_, pm1_] := With[
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
    maxcell = 10^-4,
    maxiter = 4*10^5,
    \[Kappa] = kappa,
    \[Rho] = QuantityMagnitude[(ds*Unm/Ubohr)*epsrel/(2*kappa)],
    (* for \[Mu], convert input to SI units and then take the number (should come out in kg) and convert to units of Subscript[m, 0] *)
    Egap = Abs[pm1*(ingap*Umev/Ujoul*Ujoul) - (d0*Unm/Um)*(Eperp*Uev/Uang*Um)],
	\[Mu] = QuantityMagnitude[Abs[pm1*(ingap*Umev/Ujoul*Ujoul) - (d0*Unm/Um)*(Eperp*Uev/Uang*Um)]/(2*(vF*Um/Usec)^2), "ElectronMass"],
    e = 1,
    shift = 10,
    mmax = maxm,
    radialeqs,
    solnmat,
    Ep = Eperp
    (*radialEqKeld,radialEqDir,radialEqInd,radial\[Xi]Keld,
    radial\[Xi]Dir,radial\[Xi]Ind*)
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
	{Egap, \[Mu], Ep},
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i, 1, -1}]*H2eV,
     Table[
      NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i, 1, -1}]
	}
	]
]
IndKeld[maxm_, matpars_, ingap_, Eperp_, kappa_, dee_, pm1_] := With[
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
   radialEqKInd = -(1/(\[Mu]*2)) f''[r] - (1/(2 \[Mu]*r))* 
      f'[r] - (((Pi/(2*\[Kappa]*\[Rho])*(StruveH[0, 
              Sqrt[r^2 + d^2]/\[Rho]] - 
             BesselY[0, Sqrt[r^2 + d^2]/\[Rho]]))) - (m^2/(2*\[Mu]* 
            r^2)))* f[r]; 
   radial\[Xi]KInd[m_] = 
    Simplify[
     radialEqKInd /. f -> (\[Psi][ArcTan[#]] &) /. r -> (Tan[\[Xi]]), 
     Pi/2 > \[Xi] > 0]; solnmat = {}; evTab = {}; efTab = {}; 
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
    evTab = Append[evTab, ev - shift]; efTab = Append[efTab, ef],
    {mind, 0, mmax}
    ]; {{Egap, \[Mu]},
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i,
        1, -1}]*H2eV,
    Table[efTab[[i - j + 1]][[j]], {i, Dimensions[efTab][[1]]}, {j, i,
       1, -1}]}
   ]
  ]
NormalizeEF[EF_, rmax_] := Module[
  {norm},
  norm = NIntegrate[r*(EF)^2, {r, 0, 10^6}, MinRecursion -> 5,
    MaxRecursion -> 20];
  EF/Sqrt[norm]
  ]
SiGeSuite[mmax_, matpars_, eztab_, kappa_] :=
  Module[
   {
    (* Unpack material parameters and get units straight *)
    Egap = matpars[[1]],
    (* Express bandgaps in Joules *)
    d0 = matpars[[2]]*Unm/Um,
	(* Express buckling constant in meters *)
    vF = matpars[[3]],
	(* Fermi velocity is already in m/s *)
    ds = matpars[[4]]*Unm/Um, (* Monolayer in thickness in meters *)
	erel = matpars[[5]] - 1,
    \[Kappa] = kappa,
	min,
	max
	MUTABMIN={};
	MUTABMAX={};
	EVTABMIN={};
	EVTABMAX={};
	ETRTABMIN={};
	ETRTABMAX={};
	F0TAB={};
    },
    Table[
	Export[ToString@StringForm["calcs_e`1`.m",Eperp],
	{
		Eperp,
		min={"min", DirKeld[mmax, matpars[[2 ;;]], Egap, Eperp, \[Kappa], 1]},
		max={"max", DirKeld[mmax, matpars[[2 ;;]], Egap, Eperp, \[Kappa], -1]}
	}],
	{Eperp, eztab[[1]], eztab[[2]], eztab[[3]]}
	];
	Export[ToString@StringForm["tabledone.txt","Table done, files saved. Starting analysis of data"];
	resultfiles = FileNames["*calcs*"];
	(* What to do with these separated files?
		Analyze each and do:
		List of {{Eperps},{quantities}}
		Quantities would be EVs and Etrs and also f0, alpha, abs
		*)
	Table[
		{
			thismin=Import[ToString@resultfiles[[i]]],
			MUTABMIN=Append[MUTABMIN, {thismin[[1]],thismin[[2]][[3]]}],
			EVTABMIN=Append[EVTAB, {thismin[[1]],thismin[[2]][[4]]}],
			ETRTABMIN=Append[ETRTABMIN, {thismin[[1]],Table[thismin[[2]][[4]][[j]][[j]]-thismin[[2]][[4]][[j+1]][[j+1]],{j,Length[thismin]-1}]}]
		},
		{i, Length[resultfiles]}
	];
	Export["fin-mutab-min.m",MUTABMIN//Transpose];
	Export["fin-evtab-min.m",EVTABMIN//Transpose];
	Export["fin-etrtab-min.m",EVTABMIN//Transpose];
]
getindeb[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[2]][[1]][[1]]
getindmu[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[1]][[2]]
getindeperp[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[1]]
getindgap[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[1]][[1]]
getindd[matrix_, type_, dind_] := matrix[[2]][[type]][[dind]][[1]]
getindEVs[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[2]]
getindEFs[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[3]]
calcindf0[matrix_, type_, dind_, eperp_] := 
 2*getindmu[matrix, type, dind, 
   eperp]*((getindEVs[matrix, type, dind, eperp][[2]][[2]] - 
      getindEVs[matrix, type, dind, eperp][[1]][[1]])/
    H2eV)*((1/2)*
     NIntegrate[(r^2)*(getindEFs[matrix, type, dind, eperp][[2]][[
          2]] /. \[Xi] -> 
          ArcTan[r])*(getindEFs[matrix, type, dind, eperp][[1]][[
          1]] /. \[Xi] -> ArcTan[r]), {r, 0, 10^6}, 
      MaxRecursion -> 20])^2
calcindalpha[matrix_, pars_, type_, dind_, eperp_] := 
 2*((2*Pi)/(Sqrt[4.89]*(137)))*(na/(pars[[4]]*
      getindmu[matrix, type, dind, eperp]))*
  calcindf0[matrix, type, dind, eperp]*(2/damp)
(* Why wasn't I multiplying by 2l in the argument of the exponential???? \
that worries me. Need to dig though all my old calculations for \
TMDC's now.... *)

calcindafac[matrix_, pars_, type_, dind_, eperp_] := 
 1 - Exp[-calcindalpha[matrix, pars, type, dind, eperp]*2*pars[[4]]]
getdireb[matrix_, type_, eperp_] := 
 matrix[[2]][[eperp]][[type]][[4]][[1]][[1]]
getdirmu[matrix_, type_, eperp_] := matrix[[2]][[eperp]][[type]][[3]]
getdirgap[matrix_, type_, eperp_] := 
 matrix[[2]][[eperp]][[type]][[2]]
getdireperp[matrix_, eperp_] := matrix[[2]][[eperp]][[1]]
getdirEVs[matrix_, type_, eperp_] := 
 matrix[[2]][[eperp]][[type]][[4]]
getdirEFs[matrix_, type_, eperp_] := 
 matrix[[2]][[eperp]][[type]][[5]]
calcdirf0[matrix_, type_, eperp_] := 
 2*getdirmu[matrix, type, 
   eperp]*((getdirEVs[matrix, type, eperp][[2]][[2]] - 
      getdirEVs[matrix, type, eperp][[1]][[1]])/
    H2eV)*((1/2)*
     NIntegrate[(r^2)*getdirEFs[matrix, type, eperp][[2]][[2]]*
       getdirEFs[matrix, type, eperp][[1]][[1]], {r, 0, 10^6}, 
      MaxRecursion -> 20])^2
(* I think I need to add a factor of 2 here because direct excitons \
live in one monolayer *)

calcdiralpha[matrix_, pars_, type_, eperp_] := 
 2*((4*Pi)/(Sqrt[4.89]*(137)))*(na/((pars[[4]]/B2nm)*
      getdirmu[matrix, type, eperp]))*
  calcdirf0[matrix, type, eperp]*(2/damp)
calcdirafac[matrix_, pars_, type_, eperp_] := 
 1 - Exp[-calcdiralpha[matrix, pars, type, eperp]*(pars[[4]]/B2nm)]
getTindeb[matrix_, type_, dind_, eperp_] := 
 matrix[[type]][[2]][[dind]][[2]][[eperp]][[4]][[1]][[1]][[1]] // 
  QuantityMagnitude
getTindmu[matrix_, type_, dind_, eperp_] := 
 matrix[[type]][[2]][[dind]][[2]][[eperp]][[3]] // QuantityMagnitude
getTindeperp[matrix_, type_, dind_, eperp_] := 
 matrix[[type]][[2]][[dind]][[2]][[eperp]][[1]] // QuantityMagnitude
getTindgap[matrix_, type_, dind_, eperp_] := 
 matrix[[type]][[2]][[dind]][[2]][[eperp]][[2]] // QuantityMagnitude
getTindd[matrix_, type_, dind_] := 
 matrix[[type]][[2]][[dind]][[1]] // QuantityMagnitude
getTindEVs[matrix_, type_, dind_, 
  eperp_] := (1/1000)*
   Table[matrix[[type]][[2]][[dind]][[2]][[eperp]][[4]][[i]][[j]][[
     1]], {i, 3}, {j, i}] // QuantityMagnitude
getTindEFs[matrix_, type_, dind_, eperp_] := 
 Table[matrix[[type]][[2]][[dind]][[2]][[eperp]][[4]][[i]][[j]][[
   2]], {i, 3}, {j, i}]
calcTindf0[matrix_, type_, dind_, eperp_, f_] := 
 2*getTindmu[matrix, type, dind, 
   eperp]*((getTindEVs[matrix, type, dind, eperp][[f]][[2]] - 
      getTindEVs[matrix, type, dind, eperp][[1]][[1]])/
    H2eV)*((1/2)*
     NIntegrate[(r^2)*getTindEFs[matrix, type, dind, eperp][[f]][[2]]*
       getTindEFs[matrix, type, dind, eperp][[1]][[1]], {r, 0, 10^6}, 
      MaxRecursion -> 20])^2
calcTindalpha[matrix_, pars_, type_, dind_, eperp_, f_] := 
 2*((2*Pi)/(Sqrt[4.89]*(137)))*(na/((pars[[4]]/B2nm)*
      getTindmu[matrix, type, dind, eperp]))*
  calcTindf0[matrix, type, dind, eperp, f]*(2/damp)
(* Why wasn't I multiplying by 2l in the argument of the exponential???? \
that worries me. Need to dig though all my old calculations for \
TMDC's now.... *)

calcTindafac[matrix_, pars_, type_, dind_, eperp_, f_] := 
 1 - Exp[-calcTindalpha[matrix, pars, type, dind, eperp, 
      f]*2*(pars[[4]]/B2nm)]
normindEFs[mat_, rmax_] := With[
  {
   names = {"\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(min\)]\), \
\[Xi]\[Sigma]=-1", 
     "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(min\)]\), \[Xi]\
\[Sigma]=1", 
     "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(max\)]\), \[Xi]\
\[Sigma]=-1", 
     "\!\(\*SubscriptBox[\(\[CapitalDelta]\), \(max\)]\), \[Xi]\
\[Sigma]=1"}
   },
  Table[
   {
    names[[type]],
    Table[
     {
      d,
      Table[
       {
        getindeperp[mat, type, d, ep]*
         Quantity["Volts"]/Quantity["Angstroms"],
        getindgap[mat, type, d, ep],
        getindmu[mat, type, d, ep]*Quantity["ElectronMass"],
        Table[
         {
          getindEVs[mat, type, d, ep][[i]][[j]]*1000*
           Quantity["Millielectronvolts"],
          (NIntegrate[
             r*(getindEFs[mat, type, d, ep][[i]][[j]] /. \[Xi] -> 
                  ArcTan[r])^2, {r, 0, rmax}, MinRecursion -> 5, 
             MaxRecursion -> 20]^(-1/
             2))*(getindEFs[mat, type, d, ep][[i]][[j]] /. \[Xi] -> 
              ArcTan[r])
          },
         {i, 3}, {j, i}]
        }, {ep, 17}]
      },
     {d, 10}]
    },
   {type, 1, 4}]
  ]
SGSIndSuite[pars_, ezrange_, drange_, name_, projdir_] := Module[
  {make,
   normed,
   processed,
   ei = ezrange[[1]],
   ef = ezrange[[2]],
   estep = ezrange[[3]],
   di = drange[[1]],
   df = drange[[2]],
   dstep = drange[[3]]
   },
  make = {
    Table[{dd, 
      Table[{Eperp, 
        IndKeld[2, pars[[2 ;;]], pars[[1]], Eperp, 4.89, 
         dd*lBN + (pars[[4]]/B2nm), 1]}, {Eperp, ei, ef, 
        estep}]}, {dd, di, df, dstep}],
    Table[{dd, 
      Table[{Eperp, 
        IndKeld[2, pars[[2 ;;]], pars[[1]], Eperp, 4.89, 
         dd*lBN + (pars[[4]]/B2nm), -1]}, {Eperp, ei, ef, 
        estep}]}, {dd, di, df, dstep}]
    };
  normed = normindEFs[make, 10^5];
  processed = {
    pars,
    Table[
     {
      labels[[type]],
      Table[
       {
        Nindeperp[normed, type, 1, eperp],
        Nindgap[normed, type, 1, eperp],
        Nindmu[normed, type, 1, eperp],
        Table[
         {
          NindD[normed, type, d],
          Table[
           {
            NindEV[normed, type, d, eperp, n, l],
            Sqrt[NindNINT[normed, type, d, eperp, n, l, n, l, 3]]
            },
           {n, 3}, {l, 0, n - 1}
           ],
          Table[
           {
            (NindEV[normed, type, d, eperp, nf, 1] - 
              NindEV[normed, type, d, eperp, 1, 0]),
            Nindf0[normed, type, d, eperp, 1, 0, nf, 1],
            Nindalpha[normed, type, d, eperp, 1, 0, nf, 1, pars],
            Nindafac[normed, type, d, eperp, 1, 0, nf, 1]
            },
           {nf, {2, 3}}
           ]
          },
         {d, di, df, dstep}
         ]
        },
       {eperp, 2, ((ef - ei)/estep) + 1}
       ]
      },
     {type, 2}
     ]
    };
  Export[projdir <> 
    ToString[
     StringForm["`1`-`2`-`3`_`4`_normed.m", DateObject[][[1]][[1]], 
      DateObject[][[1]][[2]], DateObject[][[1]][[3]], name]], 
   normed];
  Export[projdir <> 
    ToString[
     StringForm["`1`-`2`-`3`_`4`_processed.m", DateObject[][[1]][[1]],
       DateObject[][[1]][[2]], DateObject[][[1]][[3]], name]], 
   processed];
  Export[projdir <> 
    ToString[
     StringForm["`1`-`2`-`3`_`4`_info.txt", DateObject[][[1]][[1]], 
      DateObject[][[1]][[2]], DateObject[][[1]][[3]], name]],
   ToString[
    StringForm[
      "input parameters: ``", pars];
    ]
	]
]
