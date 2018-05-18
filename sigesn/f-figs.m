FileFig2D[data_,savename_,title_,axes_]:=Module[
	{
		xax = data[[1]],
		yax = Flatten[{data[[2]]},1],
		xlab = axes[[1]],
		ylab = axes[[2]]
	},
	(* module *)
	Export["plt-"<>ToString@savename<>".pdf",
		ListPlot[
			Table[
				Table[
					{
						xax[[j]],
						yax[[i]][[j]]
					},
					{j,Length[xax]}
				],
				{i,Length[yax]}
			],
		PlotLabel->title,
		AxesLabel->{xlab,ylab},
		ImageSize->{1440,900},
		LabelStyle->Directive[24,Black],
		PlotTheme->"Detailed"
		]
	]
]
