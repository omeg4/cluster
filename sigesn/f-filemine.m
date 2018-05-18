BuildData[]:=Module[
		{
				resultfiles = FileNames["*calcs*"],
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
			thismax=Import[resultfiles[[i]]][[2]][[2]]
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
			AFACTABMIN=Append[AFACTABMIN, {ep, {Exp[-absfromfile[thismin,2,minmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm],Exp[-absfromfile[thismin,2,minmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm]}}];
			MUTABMAX=Append[MUTABMAX, {ep,maxmu}];
			EGTABMAX=Append[EGTABMAX, {ep,maxegap}];
			EVTABMAX=Append[EVTABMAX, {ep,maxevs}];
			F0TABMAX=Append[F0TABMAX, {ep, {f0fromfile[thismax,2],f0fromfile[thismax,3]}}];
			ABSTABMAX=Append[ABSTABMAX, {ep, {absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]/B2nm,absfromfile[thismax,3,maxmu,matpars[[4]],kappa,na,damp]/B2nm}}];
			AFACTABMAX=Append[AFACTABMAX, {ep, {Exp[-absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm],Exp[-absfromfile[thismax,2,maxmu,matpars[[4]],kappa,na,damp]*matpars[[4]]/B2nm]}}];
			]
		],
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

f0fromfile[filedata_,nf_]:=Module[{mu=filedata[[1]][[2]],gsev=filedata[[2]][[1]][[1]],exev=filedata[[2]][[nf]][[2]],gsef=filedata[[3]][[1]][[1]],exef=filedata[[3]][[nf]][[2]]},2*mu*(QuantityMagnitude@UnitConvert[exev-gsev,"Hartrees"])*Nintfromfile[gsef,exef]]

absfromfile[filedata_,nf_,mu_,hl_,kappa_,na_,damp_]:=2*((4*\[Pi])/(Sqrt[kappa]*(137)))*(na/((hl/B2nm)*mu))*f0fromfile[filedata,nf]*(2/damp)
