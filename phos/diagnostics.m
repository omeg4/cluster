SetDirectory["$(pwd)"];

Table[
	Module[
		{
			imax=10, mx=mus[[1,1]], my=mus[[1,2]], pot=VKeld, chi=chiphos, d=0, eps=(eeps * 10^-5), maxi=(maxii*10^9), vn=None,
			rev,ref,ndestats,bigintout
		},
		(* Make NDEigensystem *)
		{{rev,ref},ndestats}=med[Unevaluated@ndes[imax, mx, my, pot, 4.89, chi, d, eps, maxi, vn]];
		Export[ToString@StringForm["0323-ndeout-2-eps`1`_maxi`2`.m",eeps,maxii],{{eps,maxi,vn},{rev, ref},{ndestats}}];
		(* Now for this NDE we want to test different NIntegrate approaches *)
		bigintout=Table[
			Module[
				{
					s=(10^ess),
					norm,normst,normtime,nef,
					onc,oncst,onctime,
					fox,foxst,foxtime,
					foy,foyst,foytime
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
				nef = Table[ Function[{x,y}, ref[[ i ]][ ArcTan[x], ArcTan[y] ] / Sqrt[ norm[[ i ]] ] ], {i, imax}];
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
					{norm, onc, fox, foy},
					{normst, oncst, foxst, foyst},
					{normtime, onctime, foxtime, foytime}
				};
				Export[ToString@StringForm["0323-intout-2-nde_eps`1`_maxi`2`-nint_mei`3`_minr`4`_maxr`5`_s`6`.m",eeps,maxii,mei,minr,maxr,ess],intout];
				intout
			],
			{mei,{Automatic}},{minr,{Automatic}},{maxr,{N[10000]}},{ess,{6}}
		];
		Export[ToString@StringForm["0323-bigintout-2-nde_eps`1`_maxi`2`.m",eeps,maxii],bigintout];
		bigintout
	],
	{eeps,{5,4,3,2,1,0.9,0.8,0.7,0.6,0.5}},{maxii,{1000}}
];
