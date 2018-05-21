ProcessSuite[]:=Module[
	{
		data=Import["suite.m"],
		params=Import["inp.m"][[1]],
		kappa=Import["inp.m"][[2]],
		etab=Import["inp.m"][[3]],
		dtab=Import["inp.m"][[4]]
	},
	Nd=Dimensions[ data ][[1]];
	Ne=Dimensions[ data ][[2]];
	{
	params,
	Table[(*Big outer table, 2 elements, one for each type (min/max) *)
		{
			labels[[type]],
			Table[(* This is the table of D's, with Nd elements *)
				{
					dfromsuite[ data[[Dn]][[1]][[type]] ],
					Table[
					{
						Quantity[eperpfromsuite[ data[[Dn]][[En]][[type]] ],"Volts"/"Angstroms"],
						Quantity[egapfromsuite[ data[[Dn]][[En]][[type]] ], "Electronvolts"],
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
								UnitConvert[Quantity[absfromsuite[ data[[Dn]][[En]][[type]], 1, nf, 0, 1],1/"BohrRadius"],1/"Meters"],
								afacfromsuite[ data[[Dn]][[En]][[type]], 1, nf, 0, 1 ]
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
f0fromsuite[onerun_,ni_,nf_,li_,lf_]:=2*mufromsuite[onerun]*Etrfromsuite[ onerun, ni, nf, li, lf ]*Nintfromfile[2, EFfromsuite[onerun, ni, li], EFfromsuite[onerun, nf, lf]]
absfromsuite[onerun_,params_,ni_,nf_,li_,lf_]:=2*((4*\[Pi])/(Sqrt[kappa]*(137)))*(na/((params[[4]]/B2nm)*mufromsuite[onerun]))*f0fromsuite[onerun,ni,nf,li,lf]*(2/damp)
afacfromsuite[onerun_,params_,ni_,nf_,li_,lf_]:=1-Exp[-absfromsuite[onerun,params,ni,nf,li,lf]*(params[[4]]/B2nm)]


Nintfromfile[rn_,ef1_,ef2_]:=NIntegrate[(r^rn)*ef1*ef2, {r,0,10^6},MinRecursion->10,MaxRecursion->50]
