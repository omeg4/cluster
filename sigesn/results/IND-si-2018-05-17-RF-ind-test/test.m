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
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i,
        1, -1}]*H2eV,
    Table[NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i,
       1, -1}]}
   ]
]
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
   radialEqKInd = -(1/(\[Mu]*2))*f''[r] - (1/(2 \[Mu]*r))* f'[r] - ((1/(kappa*r)) - (m^2/(2*\[Mu]* r^2)))* f[r];
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
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i,
        1, -1}]*H2eV,
    Table[NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i,
       1, -1}]}
   ]
]

SGSIndSuite[pars_, kappa_, ezrange_, drange_,funcname_] := Module[
	{
		make,
		normed,
		processed,
		ei = ezrange[[1]],
		ef = ezrange[[2]],
		estep = ezrange[[3]],
		en = (ezrange[[2]]-ezrange[[1]])/ezrange[[3]],
		di = drange[[1]],
		df = drange[[2]],
		dstep = drange[[3]],
		dn = (drange[[2]]-drange[[1]])/drange[[3]],
		min,
		max,
		MUTABMIN={},
		EGTABMIN={},
		EVTABMIN={},
		F0TABMIN={},
		ABSTABMIN={},
		AFACTABMIN={},
		MUTABMAX={},
		EGTABMAX={},
		EVTABMAX={},
		F0TABMAX={},
		ABSTABMAX={},
		AFACTABMAX={}
   },
	nd=0;
	ne=0;
	Table[
	{
		nd = nd + 1, ne = ne + 1,
		Export[
			ToString@StringForm["calcs_d`1`_e`2`.m",IntegerString[nd,10,3],IntegerString[ne,10,3]],
			{
				{"min",funcname[3, pars[[2;;]], pars[[1]], ep, kappa, d, 1]},
				{"max",funcname[3, pars[[2;;]], pars[[1]], ep, kappa, d, -1]}
			}
		]
	},
	{d,di,df,dstep},{ep,ei,ef,estep}
	];
	Export[ToString@StringForm["tabledone.txt","Table done, starting analysis"]];
	resultfiles = FileNames["*calcs*"];
	Table[
		Module[
			{
			thismin=Import[ToString@resultfiles[[i]]][[1]][[2]],
			thismax=Import[ToString@resultfiles[[i]]][[2]][[2]]
			},
			Module[
				{
					ep=thismin[[1]][[3]],
					d=thismin[[1]][[4]],
					minegap=thismin[[1]][[1]],
					minmu=thismin[[1]][[2]],
					minevs=thismin[[2]],
					minefs=thismin[[3]],
					maxegap=thismax[[1]][[1]],
					maxmu=thismax[[1]][[2]],
					maxevs=thismax[[2]],
					maxefs=thismax[[3]],
				},
			MUTABMIN=Append[MUTABMIN, {d, ep,minmu}];
			EGTABMIN=Append[EGTABMIN, {d, ep,minegap}];
			EVTABMIN=Append[EVTABMIN, {d, ep,minevs}];
			F0TABMIN=Append[F0TABMIN, {d, ep, {f0fromfile[thismin,2],f0fromfile[thismin,3]}}];
			ABSTABMIN=Append[ABSTABMIN, {d, ep, {absfromfile[thismin,2,minmu,matpars[[4]],kappa,na,damp]/B2nm,absfromfile[thismin,3,minmu,matpars[[4]],kappa,na,damp]/B2nm}}];
			AFACTABMIN=Append[AFACTABMIN, {d, ep, {Exp[-absfromfile[thismin,2,minmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm],Exp[-absfromfile[thismin,2,minmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm]}}];
			MUTABMAX=Append[MUTABMAX, {d, ep,maxmu}];
			EGTABMAX=Append[EGTABMAX, {d, ep,maxegap}];
			EVTABMAX=Append[EVTABMAX, {d, ep,maxevs}];
			F0TABMAX=Append[F0TABMAX, {d, ep, {f0fromfile[thismax,2],f0fromfile[thismax,3]}}];
			ABSTABMAX=Append[ABSTABMAX, {d, ep, {absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]/B2nm,absfromfile[thismax,3,maxmu,matpars[[4]],kappa,na,damp]/B2nm}}];
			AFACTABMAX=Append[AFACTABMAX, {d, ep, {Exp[-absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm],Exp[-absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm]}}];
			]
		]
		{i, Length[resultfiles]}
	];
	Export["fin-mutab-min.m",MUTABMIN//Transpose];
	Export["fin-egtab-min.m",EGTABMIN//Transpose];
	Export["fin-evtab-min.m",EVTABMIN//Transpose];
	Export["fin-f0tab-min.m",F0TABMIN//Transpose];
	Export["fin-abstab-min.m",ABSTABMIN//Transpose];
	Export["fin-afactab-min.m",AFACTABMIN//Transpose];
	Export["fin-mutab-max.m",MUTABMAX//Transpose];
	Export["fin-egtab-max.m",EGTABMAX//Transpose];
	Export["fin-evtab-max.m",EVTABMAX//Transpose];
	Export["fin-f0tab-max.m",F0TABMAX//Transpose];
	Export["fin-abstab-max.m",ABSTABMAX//Transpose];
	Export["fin-afactab-max.m",AFACTABMAX//Transpose];
]
Nintfromfile[ef1_,ef2_]:=NIntegrate[r*ef1*ef2, {r,0,10^6},MinRecursion->10,MaxRecursion->50]

f0fromfile[filedata_,nf_]:=Module[{mu=filedata[[1]][[2]],gsev=filedata[[2]][[1]][[1]],exef=filedata[[2]][[nf]][[2]],gsef=filedata[[3]][[1]][[1]],exef=filedata[[3]][[nf]][[2]]},2*mu*(QuantityMagnitude@UnitConvert[exev-gsev,"Hartrees"])*Nintfromfile[gsef,exef]]

absfromfile[filedata_,nf_,mu_,hl_,kappa_,na_,damp_]:=2*((4*\[Pi])/(Sqrt[kappa]*(137)))*(na/((hl/B2nm)*mu))*f0fromfile[filedata,nf]*(2/damp)
params={1.9,0.046,650000,0.4,11.9};
etab={0,2.5,0.5};
dtab={1,5,};
Export["diag1.txt","Params, Dtab, and Etab initialized"]
SGSIndSuite[params,4.89,etab,dtab,IndCoul];
Quit[]
