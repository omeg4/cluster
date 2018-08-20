data = Import['/home/mbrunetti/phos/files/2-12_sns1k.m'];
Module[
	{
		EVs = data[[2]][[1]],
		EFs = data[[2]][[2]],
		mux = mus[[1]][[1]],
		muy = mus[[1]][[2]]
	},
	ornorm = Table[
			NIntegrate[EFs[[i]]*EFs[[j]],{x,-1000,1000},{y,-1000,1000},MaxRecursion->100,WorkingPrecision->100],
		{i,10},{j,10}
	];
	f0 = Table[
		{
			(EVs[[i]]-EVs[[1]])*H2eV,
			2*mux*(EVs[[i]]-EVs[[1]])*(NIntegrate[EFs[[i]]*x*EFs[[1]],{x,-1000,1000},{y,-1000,1000},MaxRecursion->100,WorkingPrecision->100])^2,
			2*muy*(EVs[[i]]-EVs[[1]])*(NIntegrate[EFs[[i]]*y*EFs[[1]],{x,-1000,1000},{y,-1000,1000},MaxRecursion->100,WorkingPrecision->100])^2
		},
		{i,10}
	];
	plotEn = ListPlot[
		Table[
			-1000*H2eV*EVs[[i]],
			{i,10}
		],
		PlotTheme->"Detailed"
	];
	plotEtr = ListPlot[
		Table[
			1000*H2eV*(EVs[[i]]-EVs[[1]]),
			{i,2,10}
		],
		PlotTheme->"Detailed"
	];
	Export['/home/mbrunetti/phos/files/2-12_sns1k_analysis.m',{ornorm,f0}];
	Export['/home/mbrunetti/phos/files/2-12_sns1k_Ebplot.pdf',plotEn];
	Export['/home/mbrunetti/phos/files/2-12_sns1k_Etrplot.pdf',plotEtr];
]
