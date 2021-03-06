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

IndHeteroSuite[pars1_,pars2_,kappa_,ezrange_,drange_,funcname_]:=Module[
	{
		ei = ezrange[[1]],
		ef = ezrange[[2]],
		estep = ezrange[[3]],
		en = 1 + (ezrange[[2]]-ezrange[[1]])/ezrange[[3]],
		di = drange[[1]],
		df = drange[[2]],
		dstep = drange[[3]],
		dn = 1 + (drange[[2]]-drange[[1]])/drange[[3]]
	},
	Table[
	{
		funcname[3, pars1, pars2, ep, kappa, (d*lBN)+((pars1[[4]]+pars2[[4]])/(2*B2nm)), 1],
		funcname[3, pars1, pars2, ep, kappa, (d*lBN)+((pars1[[4]]+pars2[[4]])/(2*B2nm)), -1]
	},
	{d,di,df,dstep},{ep,ei,ef,estep}
	]
]
