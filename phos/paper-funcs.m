(* /*{{{*/ directory aliases *)
figdir = "../../figures/";
dfd = "~/datadrop/Dropbox/Berman-Research/dis/figs/"; (* Might as well save all figs into disstertation figures folder *)
(* /*}}}*/*)

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
combopltexp[filename_, plt_] := {
  Export[dfd <> "phos-comb-" <> filename <> ".pdf", plt, "PDF"(*,"PreviewFormat"\[Rule]"TIFF"*)],
  Export[dfd <> "phos-plt-" <> filename <> ".pdf", plt[[1]], "PDF"(*,"PreviewFormat"\[Rule]"TIFF"*)],
  Export[dfd <> "phos-leg-" <> filename <> ".pdf", plt[[2, 1]], "PDF"(*,"PreviewFormat"\[Rule]"TIFF"*)]
}
showandexp[filename_, plt_] := TableForm@{
   bg[{combopltexp[filename, plt]}, LightRed],
   plt
  }
framegrid[data_] := Grid[data, Frame -> All]
(* /*}}}*/*)

(* /*{{{*/ Key converters *)
keytox[key_] := <|
   Flatten[{
     Table[tssf["i`1`", i] -> i, {i, 12}],
     Table[tssf["j`1`", i] -> i, {i, 12}],
     (* Table[tssf["d`1`", d] -> rightd[d - 1], {d, 0, 9}]*)
		 Table[tssf["d`1`", d] -> d-1, {d,0,9}]
     }]|>[key]
xymaker[data_] := Map[
  KeyValueMap[{keytox[#1], Abs@QuantityMagnitude[#2]} &],
  data,
  {-3}
  ]
keytoleg[key_] := <|Flatten[{
     "fs" -> "Freestanding", 
     "ss" -> "\!\(\*SubscriptBox[\(SiO\), \(2\)]\) Subst.", 
     "hs" -> "\!\(\*
StyleBox[\"h\",\nFontSlant->\"Italic\"]\)-BN Subst.", "he" -> "\!\(\*
StyleBox[\"h\",\nFontSlant->\"Italic\"]\)-BN Encap.",
     "rk" -> "RK", "c" -> "Coulomb",
     Table[
      tssf["mu`1`", i] -> 
       tssf["\!\(\*SubscriptBox[\(\[Mu]\), \(`1`\)]\)", i], {i, 4}],
     "d0" -> "Direct", 
     Table[tssf["d`1`", d] -> 
       tssf["Indirect (N = `1`, D = `2` \[Angstrom])", d, 
        rightd[d - 1]], {d, 1, 9}]
     }]|>[key]
makeleg[strlist_] :=  TraditionalForm[StringRiffle[keytoleg /@ strlist, "; "]]
quantkeys = {"ndeparams", "ndeevs", "ndenefs", "onc", "fox", "foy", "alphax", "alphay", "afax", "afay"};
envkeys = {"fs", "ss", "hs", "he"};
potkeys = {"rk", "c"};
mukeys = Table[ToString@StringForm["mu`1`", i], {i, 4}];
dukeys = Table[ToString@StringForm["d`1`", i], {i, 0, 9}];
potleg = <|"rk" -> "RK, ", "c" -> "Coulomb, "|>;
envlegabbr = <|
	"fs" -> "FS",
	"ss" -> "SS",
	"hs" -> "HS",
	"he" -> "HE"
|>;
envleg = <|
	"fs" -> "Freestanding",
	"ss" -> TraditionalForm[
		"\!\(\*SubscriptBox[\(SiO\), \(2\)]\) Supported"],
	"hs" -> TraditionalForm["\!\(\*
		StyleBox[\"h\",\nFontSlant->\"Italic\"]\)-BN Supported"],
  "he" -> TraditionalForm["\!\(\*
		StyleBox[\"h\",\nFontSlant->\"Italic\"]\)-BN Encapsulated"]
|>;
dleg = <|Flatten[{"d0" -> "Direct ", Table[tssf["d`1`", i] -> tssf["N = `1`, D = `2`", i, rightd[i]], {i, 1, 9}]}]|>;
muleg = <|
	"mu1" ->
			TraditionalForm["\!\(\*SubscriptBox[\(\[Mu]\), \(1\)]\), "],
  "mu2" ->
    TraditionalForm["\!\(\*SubscriptBox[\(\[Mu]\), \(2\)]\), "],
  "mu3" ->
		TraditionalForm["\!\(\*SubscriptBox[\(\[Mu]\), \(3\)]\), "],
  "mu4" ->
    TraditionalForm["\!\(\*SubscriptBox[\(\[Mu]\), \(4\)]\), "]
|>;
quantleg = <|
	"ndeevs" -> "EVs",
	"ndenefs" -> "EFs",
  "onc" -> "Orthonormality",
	"fox" -> "Osc Str x",
  "foy" -> "Osc Str y",
	"alphax" -> "alpha x",
	"alphay" -> "alpha y",
  "afax" -> "abs fac x",
	"afay" -> "abs fac y"
|>;
xleg = <|
	"i1" -> "n",
  "d0" -> TraditionalForm[
     "\!\(\*SubscriptBox[\(N\), \(h - BN\)]\)"
	]
|>;
muinds = Table[i, {i, 4}];
dinds = Table[i, {i, 0, 9}];
ixax = <|Table[tssf["i`1`", i] -> i, {i, 12}]|>;
(* /*}}}*/*)

(* /*{{{*/ Dataset helpers (Select[Chop[...]], MergeAssoc, AssocDepth, etc...)*)

noemptyassoc = (Select[#, UnsameQ[#, <||>] &] &);

selchop[chp_: (10^-4)] := Select[(# > chp) &]
selchpy[chp_: (10^-4)] := Select[(#[[2]]>chp)&]

makeebars = Map[Around[Mean[#], {Mean[#] - Min[#], Max[#] - Mean[#]}] &];

(* flatrans does what I expect Transpose to do w/in a Dataset *)
flatrans = (Flatten[#, {{2},{1}}]&);
flatassmerge = Association[KeyValueMap[
		Module[{muky = #1, inass = #2},
    	AssociationThread[
      	Map[muky <> # &, Keys[inass]],
      	Values[inass]
      ]
]&]];

flasme = KeyValueMap[
	Module[{outkey = #1, outvals = #2},
		AssociationThread[
			Map[outkey <> # &, Keys[outvals]],
			Values[outvals]
		]
	]&]

assocdepth[ass_] := Module[
  {
   hed = Part[
       ass,
       ##
       ] & @@ {0},
   headpart = {0},
   assdepth = 0
   },
  While[hed == Association,
   assdepth += 1;
   headpart = Prepend[headpart, 1];
   hed = Part[
       ass,
       ##
       ] & @@ headpart;
   ];
  assdepth
]

ev2etr = (("ndeevs" /* QuantityMagnitude /* Abs) /* (#["i1"] - # &));
gsfox = "fox"/*"i1";
gsabx = "alphax"/*"i1";
gsafx = "afax"/*"i1";
gsfoy = "foy"/*"i1";
gsaby = "alphay"/*"i1";
gsafy = "afay"/*"i1";
(* /*}}}*/*)

(* /*{{{*/ Helper functions for ArrayPlots *)
logblendcols = {
  {0, Black},
  {10^-4, Blue},
  {10^-3, Green},
  {10^-2, Orange},
  {10^-1, Red},
  {1, Purple}
  }
barleglogcont = Reverse@Flatten[{1, Table[i*10^(-n), {n, 1, 4}, {i, {5, 1}}], 0}]
logblend[f_] := Blend[logblendcols,f]
(* /*}}}*/*)


pdiff[v1_,v2_]:= Abs[100*Abs[v1-v2]/Mean[{v1,v2}]]
pdiff[vals_List]:=With[{v1=Max[vals],v2=Min[vals]},Abs[100*Abs[v1-v2]/Mean[{v1,v2}]]]

dslpb[plopts:OptionsPattern[]]:=(Module[{
	imsi=First@Values@Merge[{Evaluate@FilterRules[{plopts},{"is"}],{"is"->{800,450}}},(#[[1]]&)], (* Custom option for ImageSize + AspectRatio *)
	masi=First@Values@Merge[{Evaluate@FilterRules[{plopts},{"ms"}],{"ms"->24}},(#[[1]]&)], (* Custom option for specifying marker size *)
	lasi=First@Values@Merge[{Evaluate@FilterRules[{plopts},{"ls"}],{"ls"->24}},(#[[1]]&)], (* Custom option for specifying label size *)
	mast=First@Values@Merge[{Evaluate@FilterRules[{plopts},{"mt"}],{"mt"->feordered}},(#[[1]]&)], (* Custom option for Marker types *)
	legopts=Values@Evaluate@FilterRules[{plopts},{lo}] (* Custom option for PlotLegends option specification. Actually can just put PlotLegends->{} as an option to ListPlot.. duh *)
	},
	ListPlot[
		#,
		PlotTheme->"Detailed",
		ImageSize->imsi,AspectRatio->(imsi[[2]]/imsi[[1]]),
		PlotStyle->randash,
		PlotMarkers->marksize[mast,masi],
		GridLinesStyle->{Thin,Gray},
		LabelStyle->Directive[lasi,Black,FontFamily->"Arial"],
		IntervalMarkers->"Bands",
		IntervalMarkersStyle->Directive[Dashing[None]],
		Evaluate@FilterRules[{plopts},Options[ListPlot]]
	]
]&)

dslpe[plopts:OptionsPattern[]]:=(Module[{
	imsi=First@Values@Merge[{Evaluate@FilterRules[{plopts},{"is"}],{"is"->{800,450}}},(#[[1]]&)], (* Custom option for ImageSize + AspectRatio *)
	masi=First@Values@Merge[{Evaluate@FilterRules[{plopts},{"ms"}],{"ms"->24}},(#[[1]]&)], (* Custom option for specifying marker size *)
	lasi=First@Values@Merge[{Evaluate@FilterRules[{plopts},{"ls"}],{"ls"->24}},(#[[1]]&)], (* Custom option for specifying label size *)
	mast=First@Values@Merge[{Evaluate@FilterRules[{plopts},{"mt"}],{"mt"->feordered}},(#[[1]]&)], (* Custom option for Marker types *)
	legopts=Values@Evaluate@FilterRules[{plopts},{lo}] (* Custom option for PlotLegends option specification. Actually can just put PlotLegends->{} as an option to ListPlot.. duh *)
	},
	ListPlot[
		#,
		PlotTheme->"Detailed",
		ImageSize->imsi,AspectRatio->(imsi[[2]]/imsi[[1]]),
		PlotStyle->randash,
		PlotMarkers->marksize[mast,masi],
		GridLinesStyle->{Thin,Gray},
		LabelStyle->Directive[lasi,Black,FontFamily->"Arial"],
		IntervalMarkers->"Ellipses",
		IntervalMarkersStyle->Directive[Dashing[None]],
		Evaluate@FilterRules[{plopts},Options[ListPlot]]
	]
]&)

(* /*{{{*/ Old attempt at a generalized plotting function *)
datncap[data_] := Module[
  {
   xydata = Normal[data[xymaker]],
   asdep, dat, leg
   },
  asdep = assocdepth[xydata];
  {dat, leg} = Flatten[
     Map[
      KeyValueMap[#2 &],
      MapIndexed[
       {#1, #2} &,
       xydata,
       {asdep}
       ],
      {0, asdep - 1}
      ],
     asdep - 1
     ] // Transpose;
  {dat, leg[[;; , ;; , 1]]}
  ]

makellp[dataset_, plopts : OptionsPattern[]] := Module[
  {
   data, legs
   },
  {data, legs} = datncap[dataset];
  ListLinePlot[
   data,
   ImageSize -> {1600, 900}, PlotTheme -> "Detailed", 
   PlotRange -> All, PlotMarkers -> {Automatic, 20},
   Evaluate@FilterRules[{plopts}, Options[ListLinePlot]],
   PlotLegends -> 
    Placed[LineLegend[Automatic, makeleg /@ legs, 
      LabelStyle -> Directive[16, Black, FontFamily -> "Arial"], 
      LegendFunction -> "Panel", LegendLayout -> {"Column", 4}], Above]
   ]
]
(* /*}}}*/*)
