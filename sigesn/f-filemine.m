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
