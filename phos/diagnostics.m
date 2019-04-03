SetDirectory["$(pwd)"];

Table[
	Module[
		{
			imax=10, mx=mus[[1,1]], my=mus[[1,2]], pot=VKeld, chi=chiphos, d=0, eps=(eeps * 10^-6), maxi=(10^15), vn=None,
			rev,ref,ndestats,bigintout
		},
		(* Make NDEigensystem *)
		{{rev,ref},ndestats}=med[Unevaluated@ndes[imax, mx, my, pot, 4.89, chi, d, eps, maxi, vn]];
		Export[ToString@StringForm["0403-ndeout-eps`1`.m",eeps],{{eps,maxi,vn},{rev, ref}}];
		Export[ToString@StringForm["0403-ndestats-eps`1`.m",eeps],ndestats];
	],
	{eeps,{10,9,8,7,6,5,4,3,2,1}}
];
