SetDirectory["$(pwd)"];

{ndeopts,{rev,ref},ndestats}=Import["../0323-short/0323-ndeout-2-eps1_maxi1000.m"];

mx = mus[[1,1]];
my = mus[[1,2]];

bigintout=Table[
	Module[
		{
			s=(10^ess),
			norm,normst,normtime,nef,
			onc,oncst,onctime,
			fox,foxst,foxtime,
			foy,foyst,foytime,
			imax=10
		},
		(* Normalize the EFs for this set of NInt options *)
		{normtime,{norm,normst}}=AbsoluteTiming[Flatten[ParallelTable[
						med[Unevaluated@NIntegrate[
							Conjugate[ref[[ i ]][ ArcTan[x], ArcTan[y] ] ] * ref[[ i ]][ ArcTan[x], ArcTan[y] ],
		{x,-s,s},{y,-s,s},
		Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
				]],
		{i,imax}],{{2}}]];
		(* Make normed EFs *)
		nef = Table[ Function[{x,y}, Evaluate[ref[[ i ]][ ArcTan[x], ArcTan[y] ] / Sqrt[ norm[[ i ]] ]] ], {i, imax}];
		(* Do a norm check *)
		{onctime,{onc,oncst}}=AbsoluteTiming[Flatten[ParallelTable[
					med[Unevaluated@NIntegrate[
						Conjugate[nef[[ j ]][ x, y ] ] * nef[[ i ]][ x, y ],
		{x,-s,s},{y,-s,s},
		Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
				]],
		{i,imax},{j,i,imax}],{{3}}]];
		(* f0x table *)
		{foxtime,{fox, foxst}}=AbsoluteTiming[Flatten[ParallelTable[
					med[Unevaluated@(2 * (rev[[j]] - rev[[i]]) * mx * NIntegrate[
						Conjugate[nef[[ j ]][ x, y ] ] * x * nef[[ i ]][ x, y ],
		{x,-s,s},{y,-s,s},
		Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
				])],
		{i,imax},{j,i,imax}],{{3}}]];
		(* f0y table *)
		{foytime,{foy, foyst}}=AbsoluteTiming[Flatten[ParallelTable[
					med[Unevaluated@(2 * (rev[[j]] - rev[[i]]) * my * NIntegrate[
						Conjugate[nef[[ j ]][ x, y ] ] * y * nef[[ i ]][ x, y ],
		{x,-s,s},{y,-s,s},
		Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
				])],
		{i,imax},{j,i,imax}],{{3}}]];
		intout={
			{eps, maxi, vn, mei, minr, maxr, s},
			{nef},
			{norm, onc, fox, foy},
			{normst, oncst, foxst, foyst},
			{normtime, onctime, foxtime, foytime}
		};
		Export[ToString@StringForm["0325-intout-nint_mei`1`_minr`2`_maxr`3`_s`4`.m",mei,minr,maxr,ess],intout];
		intout
	],
	{mei,{Automatic,5000}},{minr,{Automatic,10}},{maxr,{100,10000}},{ess,{4,5,6,7}}
];

Export["0325-bigintout.m",bigintout];
