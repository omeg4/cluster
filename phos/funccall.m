s = 10000
ns = 10000
projdir = ToString[2018-04-26_sns10K_2]
filename = ToString[sns10K_2]
result=Table[
	{
		mux,
		muy,
		Table[{nhbn,CompPhosKeld[mux,muy,4.89,rho,nhbn*lBN,10,s,ns]},{nhbn,10}]
	},
	{mux,{0.062963,0.0910365}},{muy,{0.659901,0.967742}}
]//AbsoluteTiming;
Export[
	"/home/mbrunetti/phos/results/2018-04-26_sns10K_2/2018-04-26_sns10K_2.m",
	result
];
