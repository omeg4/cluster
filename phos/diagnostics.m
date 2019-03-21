SetDirectory["$(pwd)"];

Table[
	Module[
		{
			imax=10, mx=mus[[1,1]], my=mus[[1,2]], pot=VKeld, chi=chiphos, d=0, eps=(eeps * 10^-5), maxi=(maxii*10^9), vn=None,
			rev,ref,ndestats
		},
		(* Make NDEigensystem *)
		{{rev,ref},ndestats}=med[Unevaluated[ndes[nmax, mx, my, pot, chi, d, eps, maxi, vn]]];
		Export[ToString@StringForm["0321-ndeout-eps`1`_maxi`2`.m",eeps,maxii],{{eps,maxi,vn},{rev, ref},ndestats}];
		(* Now for this NDE we want to test different NIntegrate approaches *)
		Table[
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
							med[Unevaluated[NIntegrate[
								Conjugate[ref[[ i ]][ ArcTan[x], ArcTan[y] ] ] * ref[[ i ]][ ArcTan[x], ArcTan[y] ],
				{x,-s,s},{y,-s,s},
				Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
						]]],
				{i,imax}],{{2}}]];
				(* Make normed EFs *)
				nef = Table[ Function[{x,y}, ref[[ i ]][ ArcTan[x], ArcTan[y] ] / Sqrt[ norm[[ i ]] ] ], {i, imax}];
				(* Do a norm check *)
				{onctime,{onc,oncst}}=AbsoluteTiming[Flatten[ParallelTable[
							med[Unevaluated[NIntegrate[
								Conjugate[nef[[ j ]][ x, y ] ] * nef[[ i ]][ x, y ],
				{x,-s,s},{y,-s,s},
				Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
						]]],
				{i,imax},{j,i,imax}],{{3}}]];
				(* f0x table *)
				{foxtime,{fox, foxst}}=AbsoluteTiming[Flatten[ParallelTable[
							med[Unevaluated[2 * (rev[[j]] - rev[[i]]) * mx * NIntegrate[
								Conjugate[nef[[ j ]][ x, y ] ] * x * nef[[ i ]][ x, y ],
				{x,-s,s},{y,-s,s},
				Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
						]]],
				{i,imax},{j,i,imax}],{{3}}]];
				(* f0y table *)
				{foytime,{foy, foyst}}=AbsoluteTiming[Flatten[ParallelTable[
							med[Unevaluated[2 * (rev[[j]] - rev[[i]]) * my * NIntegrate[
								Conjugate[nef[[ j ]][ x, y ] ] * y * nef[[ i ]][ x, y ],
				{x,-s,s},{y,-s,s},
				Method->{"GlobalAdaptive","MaxErrorIncreases"->mei},MinRecursion->minr,MaxRecursion->maxr
						]]],
				{i,imax},{j,i,imax}],{{3}}]];
				intout={
					{eps, maxi, vn, mei, minr, maxr, s},
					{norm, onc, fox, foy},
					{normst, oncst, foxst, foyst},
					{normtime, onctime, foxtime, foytime}
				};
				Export[ToString@StringForm["0321-intout-nde_eps`1`_maxi`2`-nint_mei`3`_minr`4`_maxr`5`_s`6`.m",eeps,maxii,mei,minr,maxr,ess],intout];
				intout
			],
			{mei,{Automatic,N[10^5]}},{minr,{Automatic,N[5],N[10]}},{maxr,{N[200],N[500],N[1000]}},{ess,3,7}
		]
	],
	{eeps,5,1,-1},{maxii,{1,1000}}
];
