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
    Table[evTab[[i - j + 1]][[j]], {i, Dimensions[evTab][[1]]}, {j, i, 1, -1}],
    Table[NormalizeEF[efTab[[i - j + 1]][[j]], 10^6], {i, Dimensions[efTab][[1]]}, {j, i, 1, -1}]
	}
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

SGSIndSuite[pars_, kappa_, ezrange_, drange_,funcname_] := Module[
	{
		make,
		normed,
		processed,
		ei = ezrange[[1]],
		ef = ezrange[[2]],
		estep = ezrange[[3]],
		en = 1 + (ezrange[[2]]-ezrange[[1]])/ezrange[[3]],
		di = drange[[1]],
		df = drange[[2]],
		dstep = drange[[3]],
		dn = 1 + (drange[[2]]-drange[[1]])/drange[[3]],
		min,
		max
	},
	Table[
	{
		funcname[3, pars[[2;;]], pars[[1]], ep, kappa, (d*lBN)+(pars[[4]]/B2nm), 1],
		funcname[3, pars[[2;;]], pars[[1]], ep, kappa, (d*lBN)+(pars[[4]]/B2nm), -1]
	},
	{d,di,df,dstep},{ep,ei,ef,estep}
	]
]
NormalizeEF[EF_, rmax_] := Module[
  {norm},
  norm = NIntegrate[r*(EF)^2, {r, 0, 10^6}, MinRecursion -> 5,
    MaxRecursion -> 20];
  EF/Sqrt[norm]
]
ProcessSuite[]:=Module[
	{
		data,
		params,
		kappa,
		etab,
		dtab,
		labels={"min","max"}
	},
	data=Import["suite.m"];
	{params,kappa,etab,dtab}=Import["inp.m"];
	Nd=Dimensions[ data ][[1]];
	Ne=Dimensions[ data ][[2]];
	{
	Flatten[{params,kappa}],
	Table[(*Big outer table, 2 elements, one for each type (min/max) *)
		{
			labels[[type]],
			Table[(* This is the table of D's, with Nd elements *)
				{
					dfromsuite[ data[[Dn]][[1]][[type]] ],
					Table[
					{
						Quantity[eperpfromsuite[ data[[Dn]][[En]][[type]] ],"Volts"/"Angstroms"],
						Quantity[egapfromsuite[ data[[Dn]][[En]][[type]] ], "Millielectronvolts"],
						Quantity[mufromsuite[ data[[Dn]][[En]][[type]] ],"ElectronMass"],
						Table[
							{
								Quantity[1000*H2eV*EVfromsuite[ data[[Dn]][[En]][[type]],n,l ],"Millielectronvolts"],
								Quantity[B2nm*r2fromsuite[ data[[Dn]][[En]][[type]],n,l ],"BohrRadius"]
							},
							{n,3},{l,0,n-1}
						],
						Table[
							{
								Quantity[1000*H2eV*Etrfromsuite[ data[[Dn]][[En]][[type]], 1, nf, 0, 1 ],"Millielectronvolts"],
								f0fromsuite[ data[[Dn]][[En]][[type]], 1, nf, 0, 1 ],
								UnitConvert[Quantity[absfromsuite[ data[[Dn]][[En]][[type]], params, kappa, 1, nf, 0, 1],1/"BohrRadius"],1/"Meters"],
								afacfromsuite[ data[[Dn]][[En]][[type]], params, kappa, 1, nf, 0, 1 ]
							},
							{nf,{2,3}}
						]
					},
					{En,Ne}
					]
				},
				{Dn,Nd}
			]
		},
		{type,2}
	]
	}
]

eperpfromsuite[onerun_]:=onerun[[1]][[3]]//QuantityMagnitude
mufromsuite[onerun_]:=onerun[[1]][[2]]//QuantityMagnitude
egapfromsuite[onerun_]:=onerun[[1]][[1]]//QuantityMagnitude
dfromsuite[onerun_]:=onerun[[1]][[4]]
EVfromsuite[onerun_,n_,l_]:=onerun[[2]][[n]][[l+1]]
EFfromsuite[onerun_,n_,l_]:=onerun[[3]][[n]][[l+1]]

Etrfromsuite[onerun_,ni_,nf_,li_,lf_]:=(EVfromsuite[ onerun, nf, lf ] - EVfromsuite[ onerun, ni, li])
normcheckfromsuite[onerun_,n_,l_]:=NIntegrate[r*(onerun[[3]][[n]][[l]]^2),{r,0,10^6},MinRecursion->10,MaxRecursion->50]
r2fromsuite[onerun_,n_,l_]:=Sqrt@NIntegrate[(r^3)*(onerun[[3]][[n]][[l+1]]^2),{r,0,10^6},MinRecursion->10,MaxRecursion->50]
f0fromsuite[onerun_,ni_,nf_,li_,lf_]:=2*mufromsuite[onerun]*Etrfromsuite[ onerun, ni, nf, li, lf ]*((1/4)*Nintfromfile[2, EFfromsuite[onerun, ni, li], EFfromsuite[onerun, nf, lf]]^2)
absfromsuite[onerun_,params_,kappa_,ni_,nf_,li_,lf_]:=2*((4*\[Pi])/(Sqrt[kappa]*(137)))*(na/((params[[4]]/B2nm)*mufromsuite[onerun]))*f0fromsuite[onerun,ni,nf,li,lf]*(2/damp)
afacfromsuite[onerun_,params_,kappa_,ni_,nf_,li_,lf_]:=1-Exp[-absfromsuite[onerun,params,kappa,ni,nf,li,lf]*(params[[4]]/B2nm)]

FPparams[proc_]:=				proc[[1]]
FPnhbn[proc_,nd_]:=				proc[[2]][[1]][[2]][[nd]][[1]]
FPeperp[proc_,ne_]:=			proc[[2]][[1]][[2]][[1]][[2]][[ne]][[1]]
FPegap[proc_,mm_,ne_]:=			proc[[2]][[mm]][[2]][[1]][[2]][[ne]][[2]]
FPmu[proc_,mm_,ne_]:=			proc[[2]][[mm]][[2]][[1]][[2]][[ne]][[3]]
FPev[proc_,mm_,nd_,ne_,n_,l_]:=	proc[[2]][[mm]][[2]][[nd]][[2]][[ne]][[4]][[n]][[l+1]][[1]]
FPr2[proc_,mm_,nd_,ne_,n_,l_]:=	proc[[2]][[mm]][[2]][[nd]][[2]][[ne]][[4]][[n]][[l+1]][[2]]
FPetr[proc_,mm_,nd_,ne_,nf_]:=	proc[[2]][[mm]][[2]][[nd]][[2]][[ne]][[5]][[nf-1]][[1]]
FPf0[proc_,mm_,nd_,ne_,nf_]:=	proc[[2]][[mm]][[2]][[nd]][[2]][[ne]][[5]][[nf-1]][[2]]
FPabs[proc_,mm_,nd_,ne_,nf_]:=	proc[[2]][[mm]][[2]][[nd]][[2]][[ne]][[5]][[nf-1]][[3]]
FPafac[proc_,mm_,nd_,ne_,nf_]:=	proc[[2]][[mm]][[2]][[nd]][[2]][[ne]][[5]][[nf-1]][[4]]
FPDdim[proc_]:=					Dimensions[ proc[[2]][[1]][[2]] ][[1]]
FPEdim[proc_]:=					Dimensions[ proc[[2]][[1]][[2]][[1]][[2]] ][[1]]

mkmuplt[prc_]:=ListPlot[
	Table[{FPeperp[prc, ne], FPmu[prc, type, ne]},{type,2},{ne,FPEdim[prc]}],
	PlotLabel->"\[Mu] vs. Eperp",
	LabelStyle->Directive[24,Black],
	FrameLabel->{"Eperp [V/\[Angstrom]]","\[Mu] [m0]"},
	PlotTheme->"Detailed",
	ImageSize->{1600,900},
	PlotLegends->SwatchLegend[{"Min","Max"}]
]

mkegplt[prc_]:=ListPlot[
	Table[{FPeperp[prc, ne], FPegap[prc, type, ne]},{type,2},{ne,FPEdim[prc]}],
	PlotLabel->"Bandgap vs. Eperp",
	LabelStyle->Directive[24,Black],
	FrameLabel->{"Eperp [V/\[Angstrom]]","Egap [meV]"},
	PlotTheme->"Detailed",
	ImageSize->{1600,900},
	PlotLegends->SwatchLegend[{"Min","Max"}]
]

mkebplt[prc_]:=ListPlot[
	Flatten[Table[{FPeperp[prc, ne], -FPev[prc, type, nd, ne, 1, 0]},{type,2},{nd,FPDdim[prc]},{ne,FPEdim[prc]}],1],
	PlotLabel->"Binding Energy vs. Eperp",
	LabelStyle->Directive[24,Black],
	FrameLabel->{"Eperp [V/\[Angstrom]]","Eb [meV]"},
	PlotTheme->"Detailed",
	ImageSize->{1600,900},
	PlotLegends->SwatchLegend[Flatten[Table[{ToString@StringForm["`1` at N=`2`",type,nd]},{type,{"Min","Max"}},{nd,FPDdim[prc]}]],LegendMarkerSize->20]
]

mkf0plt[prc_]:=ListPlot[
	Flatten[Table[{FPeperp[prc, ne], FPf0[prc, type, nd, ne, 2]},{type,2},{nd,FPDdim[prc]},{ne,FPEdim[prc]}],1],
	PlotLabel->"f0 vs. Eperp",
	LabelStyle->Directive[24,Black],
	FrameLabel->{"Eperp [V/\[Angstrom]]","f0"},
	PlotTheme->"Detailed",
	ImageSize->{1600,900},
	PlotLegends->SwatchLegend[Flatten[Table[{ToString@StringForm["`1` at N=`2`",type,nd]},{type,{"Min","Max"}},{nd,FPDdim[prc]}]],LegendMarkerSize->20]
]

mkabsplt[prc_]:=ListPlot[
	Flatten[Table[{FPeperp[prc, ne], FPabs[prc, type, nd, ne, 2]},{type,2},{nd,FPDdim[prc]},{ne,FPEdim[prc]}],1],
	PlotLabel->"\[Alpha] vs. Eperp",
	LabelStyle->Directive[24,Black],
	FrameLabel->{"Eperp [V/\[Angstrom]]","\[Alpha] [m^(-1)]"},
	PlotTheme->"Detailed",
	ImageSize->{1600,900},
	PlotLegends->SwatchLegend[Flatten[Table[{ToString@StringForm["`1` at N=`2`",type,nd]},{type,{"Min","Max"}},{nd,FPDdim[prc]}]],LegendMarkerSize->20]
]

mkafacplt[prc_]:=ListPlot[
	Flatten[Table[{FPeperp[prc, ne], FPafac[prc, type, nd, ne, 2]},{type,2},{nd,FPDdim[prc]},{ne,FPEdim[prc]}],1],
	PlotLabel->"\[ScriptCapitalA] vs. Eperp",
	LabelStyle->Directive[24,Black],
	FrameLabel->{"Eperp [V/\[Angstrom]]","\[ScriptCapitalA]"},
	PlotTheme->"Detailed",
	ImageSize->{1600,900},
	PlotLegends->SwatchLegend[Flatten[Table[{ToString@StringForm["`1` at N=`2`",type,nd]},{type,{"Min","Max"}},{nd,FPDdim[prc]}]],LegendMarkerSize->20]
]

Nintfromfile[rn_,ef1_,ef2_]:=NIntegrate[(r^rn)*ef1*ef2, {r,0,10^6},MinRecursion->10,MaxRecursion->50]
params={33,0.0676,620000,0.45,16};
etab={0,2.75,0.25};
dtab={1,9,1};
Export["diag1.txt","Params, Dtab, and Etab initialized"];
Export["inp.m",{params,4.89,etab,dtab}];
Export["suite.m",SGSIndSuite[params,4.89,etab,dtab,IndCoul]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, processing..."];
Export["proc.m",ProcessSuite[]];
prc=Import["proc.m"];
Export["muplt.pdf",mkmuplt[prc]];
Export["egapplt.pdf",mkegplt[prc]];
Export["ebplt.pdf",mkebplt[prc]];
Export["f0plt.pdf",mkf0plt[prc]];
Export["absplt.pdf",mkabsplt[prc]];
Export["afacplt.pdf",mkafacplt[prc]];
Quit[]
