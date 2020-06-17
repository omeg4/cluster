(* Not entirely sure what I wanna do here *)

(* Convert EVs to Etr's *)

(*
	"ndeevs" /* QuantityMagnitude /* Abs /* (#["i1"] - # &)
*)

dslp[frmlbl_,{lblsize_:24,mrkrsize_:24,size_:{800,450}}][opts:OptionsPattern[]]:=(ListPlot[
			#,
			PlotTheme->"Detailed",
			ImageSize->size,
			AspectRatio->(size[[2]]/size[[1]]),
			PlotStyle->randash,
			PlotMarkers->Transpose[{Riffle[emmarks4,fullmarks4],Table[mrkrsize,{8}]}],
			FrameLabel->frmlbl,
			LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
			Evaluate@FilterRules[{opts},Options[ListPlot]],
			PlotLegends->Placed[
				PointLegend[
					Automatic,
					LegendLayout->{"Row",2},
					LegendFunction->(Framed[#,RoundingRadius->5]&),
					LegendMarkerSize->mrkrsize
				],
				Above
			]
	]&)
