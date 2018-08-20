(* Set up phosphorene parameters *)
mey = { 1.2, 1.3, 0.7527, 1.12 };
mex = { 0.17, 0.1, 0.1990, 0.17 };
mhy = { 5.0, 2.8, 5.3525, 6.35 };
mhx = { 0.1, 0.2, 0.1678, 0.15 };
masses = { mey, mex, mhy, mhx };
mus = Table[
	{
		mex[[i]]*mhx[[i]]/(mex[[i]] + mhx[[i]]),
		mey[[i]]*mhy[[i]]/(mey[[i]] + mhy[[i]])
	},
	{i, Length[mey]}
];
chiphos = 0.41 / B2nm;
rho = 2 * Pi * chiphos / kappa;
