SetDirectory["$(pwd)"];
day:=DateList[][[3]];
month:=DateList[][[2]];
hour:=DateList[][[4]];
minute:=DateList[][[5]];
(* This chunk is just re-doing my long Keldysh calculations (need to check that the optical results are the same, otherwise I'm kinda boned *)
mindmax = 4;
dindmax = 8;
Table[
	Module[
		{
			rev,ref,ndestats,nef,norms,normstats,onc,oncstats,fox,foxstats,foy,foystats,
			nmax=12,
			mx=mus[[muind,1]],
			my=mus[[muind,2]],
			niopts={Method->"GlobalAdaptive",MinRecursion->1000,MaxRecursion->10^6},
			pot=VKeld,
			s=10^4,
			kappa=4.89,
			eps=10^-3,
			c=10,
			maxiter=10^15
			},
		ToString[StringForm["At `3`h`4`m on `5`/`6`, starting calculation for (RK, hbn-enc) mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Solve the schrodinger equation *)
		{{rev,ref},ndestats}=med[Unevaluated[ndesscale[nmax,mx,my,pot,kappa,chiphos,rightd[d],eps,c,maxiter,None]]];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, NDE calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* save the data *)
		(* adding a third list element to the "NDE results" where I can store computational parameters *)
		Export[ToString@StringForm["nderes_rk_he_mu`1`_d`2`.m",muind,(d+1)],{rev,ref,
				Association[
					"mu index" -> muind,
					"d" -> d,
					"dval" -> rightd[d],
					"mux" -> mx,
					"muy" -> my,
					"pot" -> pot,
					"kappa" -> kappa,
					"chi" -> chiphos,
					"nde eps" -> eps,
					"c rescale" -> c,
					"maxiter" -> maxiter,
					"int box s" -> s,
					"int opts" -> ToString[niopts]
				]
			}
		];
		Export[ToString@StringForm["ndestats_rk_he_mu`1`_d`2`.m",muind,(d+1)],ndestats];
		(* Calculate the normalization constants *)
		{norms,normstats}=Flatten[Table[med[Unevaluated[NIntegrate[
							(Conjugate[#]*#)&[ref[[n]][x,y]],
						{x,-s,s},{y,-s,s},
						Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax}],{{2}}];
		Export[ToString@StringForm["normconsts_rk_he_mu`1`_d`2`.m",muind,(d+1)],norms];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, norms calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Make the normalized EFs *)
		nef = Table[
			Function[{x,y},
				Evaluate[
					ref[[ n ]][ x,y ] / Sqrt[ norms[[ n ]] ]
				]
			],
			{n, nmax}
		];
		(* Save the normalized EFs because why not? *)
		Export[ToString@StringForm["normedefs_rk_he_mu`1`_d`2`.m",muind,(d+1)],nef];
		(* check orthonormality *)
		{onc,oncstats}=Flatten[Table[med[Unevaluated[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["onc_rk_he_mu`1`_d`2`.m",muind,(d+1)],onc];
		Export[ToString@StringForm["oncstats_rk_he_mu`1`_d`2`.m",muind,(d+1)],oncstats];
		(* calculate f0x *)
		{fox,foxstats}=Flatten[Table[med[Unevaluated[2 * mx * (rev[[m]] - rev[[n]]) * NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * x * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["fox_rk_he_mu`1`_d`2`.m",muind,(d+1)],fox];
		Export[ToString@StringForm["foxstats_rk_he_mu`1`_d`2`.m",muind,(d+1)],foxstats];
		(* Calculate f0y *)
		{foy,foystats}=Flatten[Table[med[Unevaluated[2 * my * (rev[[m]] - rev[[n]]) * NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * y * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["foy_rk_he_mu`1`_d`2`.m",muind,(d+1)],foy];
		Export[ToString@StringForm["foystats_rk_he_mu`1`_d`2`.m",muind,(d+1)],foystats];
		ToString[StringForm["`3`h`4`m on `5`/`6`: back to top of loop",muind,d,hour,minute,month,day]]>>>mylog.txt;
		ToString[StringForm["--------------------------------------------------------------------------------------------",muind,d,hour,minute,month,day]]>>>mylog.txt;
	],
	{muind,mindmax},{d,-1,dindmax}
];
ToString[StringForm["============================================================================================",muind,d,hour,minute,month,day]]>>>mylog.txt;
ToString[StringForm["`3`h`4`m on `5`/`6`: RK, h-BN encap done!!",muind,d,hour,minute,month,day]]>>>mylog.txt;
ToString[StringForm["============================================================================================",muind,d,hour,minute,month,day]]>>>mylog.txt;
(* Now do a quick version of the above for ONLY direct excitons in MLBP with kappa = 1, (3.8 + 1)/2, (4.89 + 1)/2 *)
Table[
	Module[
		{
			rev,ref,ndestats,nef,norms,normstats,onc,oncstats,fox,foxstats,foy,foystats,
			d = -1,
			nmax=12,
			mx=mus[[muind,1]],
			my=mus[[muind,2]],
			niopts={Method->"GlobalAdaptive",MinRecursion->1000,MaxRecursion->10^6},
			pot=VKeld,
			s=10^4,
			kappa=kapiter[[1]],
			envname=ToString[kapiter[[2]]],
			eps=9*10^-4,
			c=10,
			maxiter=10^15
			},
		ToString[StringForm["At `3`h`4`m on `5`/`6`, starting calculation for (RK, `7`) mu=`1`, d=`2`",muind,d,hour,minute,month,day,envname]]>>>mylog.txt;
		(* Solve the schrodinger equation *)
		{{rev,ref},ndestats}=med[Unevaluated[ndesscale[nmax,mx,my,pot,kappa,chiphos,rightd[d],eps,c,maxiter,None]]];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, NDE calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* save the data *)
		(* adding a third list element to the "NDE results" where I can store computational parameters *)
		Export[ToString@StringForm["nderes_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],{rev,ref,
				Association[
					"mu index" -> muind,
					"d" -> d,
					"dval" -> rightd[d],
					"mux" -> mx,
					"muy" -> my,
					"pot" -> pot,
					"kappa" -> kappa,
					"chi" -> chiphos,
					"nde eps" -> eps,
					"c rescale" -> c,
					"maxiter" -> maxiter,
					"int box s" -> s,
					"int opts" -> ToString[niopts]
				]
			}
		];
		Export[ToString@StringForm["ndestats_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],ndestats];
		(* Calculate the normalization constants *)
		{norms,normstats}=Flatten[Table[med[Unevaluated[NIntegrate[
							(Conjugate[#]*#)&[ref[[n]][x,y]],
						{x,-s,s},{y,-s,s},
						Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax}],{{2}}];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, norms calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Make the normalized EFs *)
		nef = Table[
			Function[{x,y},
				Evaluate[
					ref[[ n ]][ x,y ] / Sqrt[ norms[[ n ]] ]
				]
			],
			{n, nmax}
		];
		(* Save the normalized EFs because why not? *)
		Export[ToString@StringForm["normedefs_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],nef];
		(* check orthonormality *)
		{onc,oncstats}=Flatten[Table[med[Unevaluated[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["onc_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],onc];
		Export[ToString@StringForm["oncstats_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],oncstats];
		(* calculate f0x *)
		{fox,foxstats}=Flatten[Table[med[Unevaluated[2 * mx * (rev[[m]] - rev[[n]]) * NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * x * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["fox_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],fox];
		Export[ToString@StringForm["foxstats_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],foxstats];
		(* Calculate f0y *)
		{foy,foystats}=Flatten[Table[med[Unevaluated[2 * my * (rev[[m]] - rev[[n]]) * NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * y * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["foy_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],foy];
		Export[ToString@StringForm["foystats_rk_`3`_mu`1`_d`2`.m",muind,(d+1),envname],foystats];
		ToString[StringForm["`3`h`4`m on `5`/`6`: back to top of loop",muind,d,hour,minute,month,day]]>>>mylog.txt;
		ToString[StringForm["--------------------------------------------------------------------------------------------",muind,d,hour,minute,month,day]]>>>mylog.txt;
	],
	{muind,mindmax},{kapiter,{{1,ToString["fs"]},{((3.8+1)/2),ToString["ss"]},{((4.89+1)/2),ToString["hs"]}}}
];
ToString[StringForm["============================================================================================",muind,d,hour,minute,month,day]]>>>mylog.txt;
ToString[StringForm["`3`h`4`m on `5`/`6`: Direct X, RK, different kappa done!!",muind,d,hour,minute,month,day]]>>>mylog.txt;
ToString[StringForm["============================================================================================",muind,d,hour,minute,month,day]]>>>mylog.txt;
(* Do the same as the first big table here for the coulomb potential. Wasn't planning on doing direct excitons here but its a relatively minor inclusion considering how long this will run, so whatever. *)
Table[
	Module[
		{
			rev,ref,ndestats,nef,norms,normstats,onc,oncstats,fox,foxstats,foy,foystats,
			nmax=12,
			mx=mus[[muind,1]],
			my=mus[[muind,2]],
			niopts={Method->"GlobalAdaptive",MinRecursion->1000,MaxRecursion->10^6},
			pot=VCoul,
			s=10^4,
			kappa=4.89,
			eps=9*10^-4,
			c=10,
			maxiter=10^15
			},
		ToString[StringForm["At `3`h`4`m on `5`/`6`, starting calculation for (coul, hbn-enc) mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Solve the schrodinger equation *)
		{{rev,ref},ndestats}=med[Unevaluated[ndesscale[nmax,mx,my,pot,kappa,chiphos,rightd[d],eps,c,maxiter,None]]];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, NDE calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* save the data *)
		(* adding a third list element to the "NDE results" where I can store computational parameters *)
		Export[ToString@StringForm["nderes_c_he_mu`1`_d`2`.m",muind,(d+1)],{rev,ref,
				Association[
					"mu index" -> muind,
					"d" -> d,
					"dval" -> rightd[d],
					"mux" -> mx,
					"muy" -> my,
					"pot" -> pot,
					"kappa" -> kappa,
					"chi" -> chiphos,
					"nde eps" -> eps,
					"c rescale" -> c,
					"maxiter" -> maxiter,
					"int box s" -> s,
					"int opts" -> ToString[niopts]
				]
			}
		];
		Export[ToString@StringForm["ndestats_c_he_mu`1`_d`2`.m",muind,(d+1)],ndestats];
		(* Calculate the normalization constants *)
		{norms,normstats}=Flatten[Table[med[Unevaluated[NIntegrate[
							(Conjugate[#]*#)&[ref[[n]][x,y]],
						{x,-s,s},{y,-s,s},
						Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax}],{{2}}];
		ToString[StringForm["At `3`h`4`m on `5`/`6`, norms calc finished for mu=`1`, d=`2`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Make the normalized EFs *)
		nef = Table[
			Function[{x,y},
				Evaluate[
					ref[[ n ]][ x,y ] / Sqrt[ norms[[ n ]] ]
				]
			],
			{n, nmax}
		];
		(* Save the normalized EFs because why not? *)
		Export[ToString@StringForm["normedefs_c_he_mu`1`_d`2`.m",muind,(d+1)],nef];
		(* check orthonormality *)
		{onc,oncstats}=Flatten[Table[med[Unevaluated[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["onc_c_he_mu`1`_d`2`.m",muind,(d+1)],onc];
		Export[ToString@StringForm["oncstats_c_he_mu`1`_d`2`.m",muind,(d+1)],oncstats];
		(* calculate f0x *)
		{fox,foxstats}=Flatten[Table[med[Unevaluated[2 * mx * (rev[[m]] - rev[[n]]) * NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * x * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["fox_c_he_mu`1`_d`2`.m",muind,(d+1)],fox];
		Export[ToString@StringForm["foxstats_c_he_mu`1`_d`2`.m",muind,(d+1)],foxstats];
		(* Calculate f0y *)
		{foy,foystats}=Flatten[Table[med[Unevaluated[2 * my * (rev[[m]] - rev[[n]]) * NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * y * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["foy_c_he_mu`1`_d`2`.m",muind,(d+1)],foy];
		Export[ToString@StringForm["foystats_c_he_mu`1`_d`2`.m",muind,(d+1)],foystats];
		ToString[StringForm["`3`h`4`m on `5`/`6`: back to top of loop",muind,d,hour,minute,month,day]]>>>mylog.txt;
		ToString[StringForm["--------------------------------------------------------------------------------------------",muind,d,hour,minute,month,day]]>>>mylog.txt;
	],
	{muind,mindmax},{d,-1,dindmax}
];
