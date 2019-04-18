SetDirectory["$(pwd)"];
day:=DateList[][[3]];
month:=DateList[][[2]];
hour:=DateList[][[4]];
minute:=DateList[][[5]];
Table[
	Module[
		{
			rev,ref,ndestats,nef,norms,normstats,onc,oncstats,fox,foxstats,foy,foystats,
			nmax,
			mx=mus[[muind,1]],
			my=mus[[muind,2]],
			niopts={Method->"GlobalAdaptive",MinRecursion->1000,MaxRecursion->10^6},
			s=10^4
		},
		ToString[StringForm["Starting calculation for mu=`1`, d=`2`; at `3`h`4`m on `5`/`6`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Solve the schrodinger equation *)
		{{rev,ref},ndestats}=med[Unevaluated[ndes[nmax,mx,my,VKeld,4.89,chiphos,rightd[d],9*10^-6,10^12,None]]];
		ToString[StringForm["NDE calc finished for mu=`1`, d=`2`; at `3`h`4`m on `5`/`6`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* save the data *)
		Export[ToString@StringForm["nderes_mu`1`_d`2`.m",muind,(d+1)],{rev,ref}];
		Export[ToString@StringForm["ndestats_mu`1`_d`2`.m",muind,(d+1)],ndestats];
		(* normalize the EFs *)
		{norms,normstats}=Flatten[Table[med[Unevaluated[NIntegrate[
					Conjugate[ ref[[n]][ ArcTan[x], ArcTan[y] ] ] * ref[[n]][ ArcTan[x], ArcTan[y] ],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax}],{{2}}];
		ToString[StringForm["Norms calc finished for mu=`1`, d=`2`; at `3`h`4`m on `5`/`6`",muind,d,hour,minute,month,day]]>>>mylog.txt;
		(* Make the normalized EFs *)
		nef = Table[ Function[{x,y}, Evaluate[ref[[ n ]][ ArcTan[x], ArcTan[y] ] / Sqrt[ norms[[ n ]] ]] ], {n, nmax}];
		(* check orthonormality *)
		{onc,oncstats}=Flatten[Table[med[Unevaluated[NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["onc_mu`1`_d`2`.m",muind,(d+1)],onc];
		Export[ToString@StringForm["oncstats_mu`1`_d`2`.m",muind,(d+1)],oncstats];
		(* calculate f0x *)
		{fox,foxstats}=Flatten[Table[med[Unevaluated[2 * mx * (rev[[m]] - rev[[n]]) * NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * x * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["fox_mu`1`_d`2`.m",muind,(d+1)],fox];
		Export[ToString@StringForm["foxstats_mu`1`_d`2`.m",muind,(d+1)],foxstats];
		(* Calculate f0y *)
		{foy,foystats}=Flatten[Table[med[Unevaluated[2 * my * (rev[[m]] - rev[[n]]) * NIntegrate[
					Conjugate[ nef[[m]][x,y] ] * y * nef[[n]][x,y],
					{x,-s,s},{y,-s,s},
					Evaluate@FilterRules[{niopts},Options[NIntegrate]]
				]]],
		{n,nmax},{m,n,nmax}],{{3}}];
		Export[ToString@StringForm["foy_mu`1`_d`2`.m",muind,(d+1)],foy];
		Export[ToString@StringForm["foystats_mu`1`_d`2`.m",muind,(d+1)],foystats];
	],
	{muind,4},{d,-1,6}
];
