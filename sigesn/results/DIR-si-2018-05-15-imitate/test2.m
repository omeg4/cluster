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
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i, 1, -1}]*H2eV*1000*Quantity["Millielectronvolts"],
	Table[NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i, 1, -1}]
	}
	]
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
	max,
    },
	ne=0;
    Table[{ne = ne + 1,
	Export[ToString@StringForm["calcs_e`1`.m",IntegerString[ne,10,3]],
	{
		{"min", DirKeld[mmax, matpars[[2 ;;]], Egap, Eperp, \[Kappa], 1]},
		{"max", DirKeld[mmax, matpars[[2 ;;]], Egap, Eperp, \[Kappa], -1]}
	}]},
	{Eperp, eztab[[1]], eztab[[2]], eztab[[3]]}
	];
	Export["tabledone.txt","Table done, files saved. Starting analysis of data"];
]
NormalizeEF[EF_, rmax_] := Module[
  {norm},
  norm = NIntegrate[r*(EF)^2, {r, 0, 10^6}, MinRecursion -> 5,
    MaxRecursion -> 20];
  EF/Sqrt[norm]
]
BuildData[]:=Module[
		{
				resultfiles = FileNames["*results*"],
				inps = Import["inp.m"],
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
	Table[
		Module[
			{
			thismin=Import[resultfiles[[i]]][[1]][[2]],
			thismax=Import[resultfiles[[i]]][[2]][[2]],
			kappa=inps[[2]],
			matpars=inps[[1]]
			},
			Module[
				{
					ep=thismin[[1]][[3]],
					minegap=thismin[[1]][[1]],
					minmu=thismin[[1]][[2]],
					minevs=thismin[[2]],
					minefs=thismin[[3]],
					maxegap=thismax[[1]][[1]],
					maxmu=thismax[[1]][[2]],
					maxevs=thismax[[2]],
					maxefs=thismax[[3]]
				},
			MUTABMIN=Append[MUTABMIN, {ep,minmu}];
			EGTABMIN=Append[EGTABMIN, {ep,minegap}];
			EVTABMIN=Append[EVTABMIN, {ep,minevs}];
			F0TABMIN=Append[F0TABMIN, {ep, {f0fromfile[thismin,2],f0fromfile[thismin,3]}}];
			ABSTABMIN=Append[ABSTABMIN, {ep, {absfromfile[thismin,2,minmu,matpars[[4]],kappa,na,damp]/B2nm,absfromfile[thismin,3,minmu,matpars[[4]],kappa,na,damp]/B2nm}}];
			AFACTABMIN=Append[AFACTABMIN, {ep, {1 - Exp[-absfromfile[thismin,2,minmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm], 1 - Exp[-absfromfile[thismin,2,minmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm]}}];
			MUTABMAX=Append[MUTABMAX, {ep,maxmu}];
			EGTABMAX=Append[EGTABMAX, {ep,maxegap}];
			EVTABMAX=Append[EVTABMAX, {ep,maxevs}];
			F0TABMAX=Append[F0TABMAX, {ep, {f0fromfile[thismax,2],f0fromfile[thismax,3]}}];
			ABSTABMAX=Append[ABSTABMAX, {ep, {absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]/B2nm,absfromfile[thismax,3,maxmu,matpars[[4]],kappa,na,damp]/B2nm}}];
			AFACTABMAX=Append[AFACTABMAX, {ep, {1 - Exp[-absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm], 1 - Exp[-absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm]}}];
			]
		],
		{i, Length[resultfiles]}
	];
	Export["fin-mutab.m",{MUTABMIN//Transpose,MUTABMAX//Transpose}];
	Export["fin-egtab.m",{EGTABMIN//Transpose,EGTABMAX//Transpose}];
	Export["fin-evtab.m",{EVTABMIN//Transpose,EVTABMAX//Transpose}];
	Export["fin-f0tab.m",{F0TABMIN//Transpose,F0TABMAX//Transpose}];
	Export["fin-abstab.m",{ABSTABMIN//Transpose,ABSTABMAX//Transpose}];
	Export["fin-afactab.m",{AFACTABMIN//Transpose,AFACTABMAX//Transpose}];
]

Nintfromfile[rn_,ef1_,ef2_]:=NIntegrate[(r^rn)*ef1*ef2, {r,0,10^6},MinRecursion->10,MaxRecursion->50]

f0fromfile[filedata_,nf_]:=Module[{mu=filedata[[1]][[2]],gsev=filedata[[2]][[1]][[1]],exev=filedata[[2]][[nf]][[2]],gsef=filedata[[3]][[1]][[1]],exef=filedata[[3]][[nf]][[2]]},2*mu*(QuantityMagnitude@UnitConvert[exev-gsev,"Hartrees"])*((1/4)*Nintfromfile[2,gsef,exef]^2)]

absfromfile[filedata_,nf_,mu_,hl_,kappa_,na_,damp_]:=2*((4*\[Pi])/(Sqrt[kappa]*(137)))*(na/((hl/B2nm)*mu))*f0fromfile[filedata,nf]*(2/damp)
SetDirectory["/home/mbrunetti/cluster/sigesn/results/DIR-si-2018-05-15-imitate"]
BuildData[];
Quit[]
