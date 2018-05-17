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
