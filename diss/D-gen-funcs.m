(* /*{{{*/ Standard plot styles/labels/fonts/markers/colors/etc. *)
fullmarks4 = {"\[FilledSquare]", "\[FilledCircle]", 
   "\[FilledDiamond]", "\[FilledRectangle]"};
emmarks4 = {"\[EmptySquare]", "\[EmptyCircle]", "\[EmptyDiamond]", 
   "\[EmptyRectangle]"};
feordered = Flatten[{fullmarks4,emmarks4}];
feriff = Riffle[fullmarks4,emmarks4];
marksize[marks_,size_] := Module[{len=Length[marks]},Transpose[{marks,Table[size,{len}]}]]
labsty = {Black, 24, FontFamily -> "Arial"};
cols4 = {Red, Blue, Purple, Orange};
(* /*{{{*/ Special fold for 'randash' *)
randash = Transpose@{
	{
		RGBColor[0.65, 0., 0.],
		RGBColor[0.212151, 0.39271, 0.8],
		RGBColor[1, 0.212037, 0.300311],
		RGBColor[0.111025, 0.56, 0.418696],
		RGBColor[0.752461, 0.362306, 0.125339],
		RGBColor[0.7055434835026511, 0.5895048315002945, 0.],
		RGBColor[0.0504678, 0.526626, 0.627561],
		RGBColor[0.435888, 0.259065, 0.71028],
		RGBColor[0.4078757751993936, 0.27540579780035374,0.780310562001516],
		RGBColor[0.784922, 0.524612, 0.0407096],
		RGBColor[0.515278, 0.224, 0.530342],
		RGBColor[0.12461179749810647, 0.4603105620015149,0.7095341569977278]
	},
	{
		Dashing[{Large, Tiny, Small}],
		Dashing[{Large, Medium, Tiny}],
		Dashing[{Tiny}],
		Dashing[{Tiny, Small}],
		Dashing[{Large}],
		Dashing[{Large}],
		Dashing[{Large, Tiny, Medium}],
		Dashing[{Large, Small, Tiny, Medium}],
		Dashing[{Large}],
		Dashing[{Tiny, Small}],
		Dashing[{Small, Medium, Large}],
		Dashing[{Large, Small}]
	}
}
(* /*}}}*/*)
(* /*}}}*/*)

(* /*{{{*/ Output formatting helpers (background color, columns, framed grid, etc. *)
bg[expr_, clr_RGBColor] := Framed[expr, FrameStyle -> None, Background -> clr]
bg[clr_RGBColor, expr_] := Framed[expr, FrameStyle -> None, Background -> clr]
ccol[figs_] := Column[figs, Center]
nvv[assoc_] := Values[Values[Normal[assoc]]]
dfd = "~/datadrop/Dropbox/Berman-Research/dis/figs/"; (* Might as well save all figs into disstertation figures folder *)
combopltexp[filename_, plt_] := {
  Export[dfd <> "phos-comb-" <> filename <> ".pdf", plt, "PDF"(*,"PreviewFormat"\[Rule]"TIFF"*)],
  Export[dfd <> "phos-plt-" <> filename <> ".pdf", plt[[1]], "PDF"(*,"PreviewFormat"\[Rule]"TIFF"*)],
  Export[dfd <> "phos-leg-" <> filename <> ".pdf", plt[[2, 1]], "PDF"(*,"PreviewFormat"\[Rule]"TIFF"*)]
}
showandexp[filename_, plt_] := TableForm@{
   bg[{combopltexp[filename, plt]}, LightRed],
   plt
  }
framegrid[data_,opts:OptionsPattern[]] := Grid[data, Frame -> All, ItemSize->Full, Evaluate@FilterRules[{opts},Options[Grid]]]
(* /*}}}*/*)

tssf[x___]:=ToString@StringForm[x]
med[f_]:=({#["Result"],KeyTake[#,{"AbsoluteTiming","MessagesText","Timing"}]}&[KeyTake[EvaluationData[f],{"Result","AbsoluteTiming","MessagesText","Timing"}]])
med[f_,x___]:=({#["Result"],KeyTake[#,{"AbsoluteTiming","MessagesText","Timing"}]}&[KeyTake[EvaluationData[f[x]],{"Result","AbsoluteTiming","MessagesText","Timing"}]])
bg[expr_, clr_RGBColor] := Framed[expr, FrameStyle -> None, Background -> clr]
bg[clr_RGBColor, expr_] := Framed[expr, FrameStyle -> None, Background -> clr]
fulldate:=DateString[{"Month","/","Day"," @ ","Time"," | "}]
time:=DateString[{"Time"}]

nestedgridWheaders[data_,{innerrow_, innercolumn_,innergropts : OptionsPattern[]},{outerrow_, outercolumn_,outergropts : OptionsPattern[]}]:=Grid[
  Join[
   {PadLeft[outerrow, Length[outerrow] + 1, " "]} // Transpose,
   Join[
    {outercolumn},
    Map[
     Grid[
       Join[
        {PadLeft[innerrow, Length[innerrow] + 1, " "]} // Transpose,
        Join[
         {innercolumn},
         #
         ],
        2
        ],
       Dividers -> {{False, True}, {False, True}},
       ItemStyle -> 
        Directive[FontSize -> 14, Black, FontFamily -> "Arial"],
       Evaluate@FilterRules[{innergropts}, Options[Grid]]
       ] &,
     data,
     {2}
     ]
    ],
   2
   ],
  Dividers -> {{False, True}, {False, True}},
  ItemStyle -> Directive[FontSize -> 14, Black, FontFamily -> "Arial"],
  Evaluate@FilterRules[{outergropts}, Options[Grid]]
]

gridWheaders[data_, rowheads_, columnheads_, gropts : OptionsPattern[]] := Grid[
  Join[
   {PadLeft[rowheads, Length[rowheads] + 1, " "]} // Transpose,
   Join[
    {columnheads},
    data
    ],
   2
   ],
  Dividers -> {{False, True}, {False, True}},
  ItemStyle -> Directive[FontSize -> 20, Black, FontFamily -> "Arial"],
  Evaluate@FilterRules[{gropts}, Options[Grid]]
]
