SetDirectory["$(pwd)"];


mx = mus[[1,1]];
my = mus[[1,2]];

ToString[(StringForm["Script starting on `1`/`2` at `3`h-`4`m", #[[2]], #[[3]], #[[4]], #[[5]]] &[DateList[]])] >>> mylog.txt;

Table[
	{ndeopts,{rev,ref},ndestats}=Import[ToString@StringForm["nde`1`.m",]];
	ToString[(StringForm["Import complete on `1`/`2` at `3`h-`4`m", #[[2]], #[[3]], #[[4]], #[[5]]] &[DateList[]])] >>> mylog.txt;
	Table[
		Module[
			{
				s=10^ess,
				norm,normst,normtime,nef,
				onc,oncst,onctime,
				fox,foxst,foxtime,
				foy,foyst,foytime,
				imax=10,
				mei=Automatic,
				maxr=10000
			},
			(* Normalize the EFs for this set of NInt options *)
			{normtime,{norm,normst}}=AbsoluteTiming[Flatten[Table[
							med[Unevaluated@NIntegrate[
								Conjugate[ref[[ i ]][ ArcTan[x], ArcTan[y] ] ] * ref[[ i ]][ ArcTan[x], ArcTan[y] ],
			{x,-s,s},{y,-s,s},
			Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
					]],
			{i,imax}],{{2}}]];
			ToString[(StringForm["Norm complete on `1`/`2` at `3`h-`4`m", #[[2]], #[[3]], #[[4]], #[[5]]] &[DateList[]])] >>> mylog.txt;
			(* Make normed EFs *)
			nef = Table[ Function[{x,y}, Evaluate[ref[[ i ]][ ArcTan[x], ArcTan[y] ] / Sqrt[ norm[[ i ]] ]] ], {i, imax}];
			(* Do a norm check *)
			{onctime,{onc,oncst}}=AbsoluteTiming[Flatten[Table[
						med[Unevaluated@NIntegrate[
							Conjugate[nef[[ j ]][ x, y ] ] * nef[[ i ]][ x, y ],
			{x,-s,s},{y,-s,s},
			Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
					]],
			{i,imax},{j,i,imax}],{{3}}]];
			ToString[(StringForm["Norm check complete on `1`/`2` at `3`h-`4`m", #[[2]], #[[3]], #[[4]], #[[5]]] &[DateList[]])] >>> mylog.txt;
			(* f0x table *)
			{foxtime,{fox, foxst}}=AbsoluteTiming[Flatten[Table[
						med[Unevaluated@(2 * (rev[[j]] - rev[[i]]) * mx * (NIntegrate[
							Conjugate[nef[[ j ]][ x, y ] ] * x * nef[[ i ]][ x, y ],
			{x,-s,s},{y,-s,s},
			Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
					]^2))],
			{i,imax},{j,i,imax}],{{3}}]];
			ToString[(StringForm["fox complete on `1`/`2` at `3`h-`4`m", #[[2]], #[[3]], #[[4]], #[[5]]] &[DateList[]])] >>> mylog.txt;
			(* f0y table *)
			{foytime,{foy, foyst}}=AbsoluteTiming[Flatten[Table[
						med[Unevaluated@(2 * (rev[[j]] - rev[[i]]) * my * (NIntegrate[
							Conjugate[nef[[ j ]][ x, y ] ] * y * nef[[ i ]][ x, y ],
			{x,-s,s},{y,-s,s},
			Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
					]^2))],
			{i,imax},{j,i,imax}],{{3}}]];
			ToString[(StringForm["foy complete on `1`/`2` at `3`h-`4`m", #[[2]], #[[3]], #[[4]], #[[5]]] &[DateList[]])] >>> mylog.txt;
			intout={
				{eps, maxi, vn, mei, minr, maxr, s},
				{nef},
				{norm, onc, fox, foy},
				{normtime, onctime, foxtime, foytime}
			};
			Export[ToString@StringForm["0403-nintdiag2-minr`1`-s`2`",minr,ess],intout];
			Export[ToString@StringForm["0403-nintstats2-minr`1`-s`2`",minr,ess],{normst, oncst, foxst, foyst}];
		],
		{ess,{7,6,8}},{minr,{20,50,100}}
	],
	{eeps,{9,10}}
];
