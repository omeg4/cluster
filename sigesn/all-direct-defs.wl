(* ::Package:: *)

(* ::Input::Initialization:: *)
SetDirectory@NotebookDirectory[];
Import["f-filemine.m"];


(* ::Input::Initialization:: *)
bhprc=Import["results/big-sihbn-DIR-si-2018-05-27/proc.m"];
shprc=Import["results/small-sihbn-DIR-si-2018-05-27/proc.m"];
fssiprc=Import["results/fssi-DIR-si-2018-06-05/reproc.m"];
fsgeprc=Import["results/fsge-DIR-ge-2018-06-05/reproc.m"];
fssnprc=Import["results/fssn-DIR-sn-2018-06-05/reproc.m"];


(* ::Input::Initialization:: *)
HB=Quantity["ReducedPlanckConstant"];
HH=Quantity["PlanckConstant"];
CC=Quantity["SpeedOfLight"];
KK=1/(4*Pi*Quantity["ElectricConstant"]);
EE=Quantity["ElementaryCharge"];
KB=Quantity["BoltzmannConstant"];


configs={"mono","open"};
methods={"old","new"};


(* ::Input::Initialization:: *)
QConc[conc_]:=Quantity[conc,"Meters"^-2]
QFreq[freq_]:=Quantity[freq,"Seconds"^-1]


(* ::Input::Initialization:: *)
newvarplt[data_,bounds_,frmlbl_,lgndlbl_,colors_,size_:{1600,900},lblsize_:30,plopts:OptionsPattern[]]:=
{
Framed@LineLegend[colors,lgndlbl,LegendLayout->{"Row",1},LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],LegendFunction->"Frame"],
Plot[
data,
bounds,
PlotTheme->"Detailed",
FrameLabel->frmlbl,
GridLinesStyle->Directive[Gray,Thin],
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
PlotLegends->None,
Evaluate@FilterRules[{plopts},Options[Plot]]
]
}


(* ::Input::Initialization:: *)
newlineplt[data_,frmlbl_,lgndlbl_,colors_,size_:{1600,900},lblsize_:30,plopts:OptionsPattern[]]:=
{
Framed@LineLegend[colors,lgndlbl,LegendLayout->{"Row",1},LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],LegendFunction->"Frame",LegendMarkerSize->{75,25}],
ListLinePlot[
data,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
FrameLabel->frmlbl,
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
Evaluate@FilterRules[{plopts},Options[ListLinePlot]]
]
}


(* ::Input::Initialization:: *)
newplt[data_,frmlbl_,lgndlbl_,colors_,{markers_,msize_:24,lmsize_:30},size_:{1600,900},lblsize_:30,plopts:OptionsPattern[]]:=
{
Framed@PointLegend[colors,lgndlbl,LegendMarkers->Table[Style[ToString@markers[[i]],lmsize],{i,Length[markers]}],LegendMarkerSize->lmsize,LegendLayout->{"Row",1},LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],LegendFunction->"Frame"],
ListPlot[
data,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
FrameLabel->frmlbl,
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
PlotMarkers->Table[{markers[[i]],msize},{i,Length[markers]}],
Evaluate@FilterRules[{plopts},Options[ListPlot]]
]
}


(* ::Input::Initialization:: *)
insetvarplt[data_,bounds_,colors_,ipos_,size_:{800,450},lblsize_:24,plopts:OptionsPattern[]]:=
Inset[
Framed[
Plot[
data,
bounds,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
PlotLegends->None,
Evaluate@FilterRules[{plopts},Options[Plot]]
],
Background->White,FrameStyle->None
],
Scaled[{ipos[[1]],ipos[[2]]}],Scaled[{ipos[[3]],ipos[[4]]}]
]


(* ::Input::Initialization:: *)
insetlineplt[data_,colors_,ipos_,size_:{800,450},lblsize_:24,plopts:OptionsPattern[]]:=
Inset[
Framed[
ListLinePlot[
data,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
Evaluate@FilterRules[{plopts},Options[ListLinePlot]]
],
Background->White,FrameStyle->None
],
Scaled[{ipos[[1]],ipos[[2]]}],Scaled[{ipos[[3]],ipos[[4]]}]
]


(* ::Input::Initialization:: *)
insetlistplt[data_,colors_,markers_,ipos_,size_:{800,450},lblsize_:24,plopts:OptionsPattern[]]:=
Inset[
Framed[
ListPlot[
data,
PlotTheme->"Detailed",
GridLinesStyle->Directive[Gray,Thin],
LabelStyle->Directive[lblsize,Black,FontFamily->"Arial"],
ImageSize->size,
AspectRatio->size[[2]]/size[[1]],
PlotStyle->colors,
PlotMarkers->Table[{markers[[i]],lblsize},{i,Length[markers]}],
Evaluate@FilterRules[{plopts},Options[ListPlot]]
],
Background->White,FrameStyle->None
],
Scaled[{ipos[[1]],ipos[[2]]}],Scaled[{ipos[[3]],ipos[[4]]}]
]


(* ::Input::Initialization:: *)
NoNegEtr[prc_,type_]:=Module[
{
i=1
},
While[
QuantityMagnitude@Exres[prc,type,i]<0,
i++
];
i
]


(* ::Input::Initialization:: *)
nsoleffeps[prc_,t_,ep_]:=NSolve[-FPev[prc,t,1,ep,1,0]==(2*(Quantity["Coulomb Constant"]*Quantity["ElementaryCharge"]^2)^2*FPmu[prc,t,ep])/((Quantity["ReducedPlanckConstant"]^2)*(\[Epsilon]^2))&&\[Epsilon]>0,\[Epsilon]][[1]][[1]]//Last


(* ::Input::Initialization:: *)
newlines5={
Directive[Black,AbsoluteThickness[2.5]],
Directive[Red,AbsoluteThickness[2.5],AbsoluteDashing[{1,10}]],
Directive[Blue,AbsoluteThickness[2.5],AbsoluteDashing[{4,8}]],
Directive[Purple,AbsoluteThickness[2.5],AbsoluteDashing[{6,12}]],
Directive[RGBColor[0,0.6,0.2],AbsoluteThickness[2.5],AbsoluteDashing[{2,7,12,7}]]
};
newlines10={
Directive[Black,AbsoluteThickness[2.5]],Directive[Black,AbsoluteThickness[1.5]],
Directive[Red,AbsoluteThickness[2.5],AbsoluteDashing[{1,10}]],Directive[Red,AbsoluteThickness[2.5],AbsoluteDashing[{10,1}]],
Directive[Blue,AbsoluteThickness[2.5],AbsoluteDashing[{4,8}]],Directive[Blue,AbsoluteThickness[2.5],AbsoluteDashing[{8,4}]],
Directive[Purple,AbsoluteThickness[2.5],AbsoluteDashing[{6,12}]],Directive[Purple,AbsoluteThickness[2.5],AbsoluteDashing[{12,6}]],
Directive[RGBColor[0,0.6,0.2],AbsoluteThickness[2.5],AbsoluteDashing[{2,7,12,7}]],Directive[RGBColor[0,0.6,0.2],AbsoluteThickness[2.5],AbsoluteDashing[{7,12,7,2}]]
};
newlines4={
Directive[Red,AbsoluteThickness[1.5]],
Directive[Blue,AbsoluteThickness[1.5]],
Directive[Black,AbsoluteThickness[1.5],AbsoluteDashing[{8,4}]],
Directive[Purple,AbsoluteThickness[1.5],AbsoluteDashing[{8,4}]]
};


(* ::Input::Initialization:: *)
allmatsAB={"FS Si, B","FS Si, A","FS Ge, B","FS Ge, A","FS Sn, B","FS Sn, A","Type I, B","Type I, A","Type II, B","Type II, A"};
allmatsA={"FS Si","FS Ge","FS Sn","Type I","Type II"};
eplabel="\!\(\*SubscriptBox[\(E\), \(\[Perpendicular]\)]\) [V/\[Angstrom]]";


(* ::Input::Initialization:: *)
colors3={Red,Blue,Purple};
markers3={"\[FilledSquare]","\[FilledCircle]","\[FilledUpTriangle]"};
colors4={Red,Red,Blue,Blue};
markers4={"\[FilledSquare]","\[EmptySquare]","\[FilledCircle]","\[EmptyCircle]"};
colors6={Red,Red,Blue,Blue,Purple,Purple};
markers6={"\[FilledSquare]","\[EmptySquare]","\[FilledCircle]","\[EmptyCircle]","\[FilledUpTriangle]","\[EmptyUpTriangle]"};
colors8={Red,Red,Blue,Blue,Purple,Purple,Orange,Orange};
markers8={"\[FilledSquare]","\[EmptySquare]","\[FilledCircle]","\[EmptyCircle]","\[FilledUpTriangle]","\[EmptyUpTriangle]","\[FilledDownTriangle]","\[EmptyDownTriangle]"};
colors20={Red,Red,Blue,Blue,Purple,Purple,Orange,Orange,Green,Green,Red,Red,Blue,Blue,Purple,Purple,Orange,Orange,Green,Green};
markers20={"\[FilledSquare]","\[EmptySquare]","\[FilledCircle]","\[EmptyCircle]","\[FilledUpTriangle]","\[EmptyUpTriangle]","\[FilledDownTriangle]","\[EmptyDownTriangle]","\[FilledDiamond]","\[EmptyDiamond]","\[FilledDiamond]","\[EmptyDiamond]","\[FilledDownTriangle]","\[EmptyDownTriangle]","\[FilledUpTriangle]","\[EmptyUpTriangle]","\[FilledCircle]","\[EmptyCircle]","\[FilledSquare]","\[EmptySquare]"};


(* ::Input::Initialization:: *)
plotehmass[pm1_,pars_][ez_]:=Module[
{
dso=UnitConvert[Quantity[pars[[1]],"Millielectronvolts"],"Joules"],
d0=UnitConvert[Quantity[pars[[2]],"Nanometers"],"Meters"],
vf=Quantity[pars[[3]],"Meters"/"Seconds"]
},
UnitConvert[Abs[(pm1*dso-d0*Quantity[ez,"Electronvolts"/"Angstroms"])]/(vf^2),"ElectronMass"]
]


(* ::Input::Initialization:: *)
plotgap[pm1_,pars_][ez_]:=Module[
{
dso=UnitConvert[Quantity[pars[[1]],"Millielectronvolts"],"Joules"],
d0=UnitConvert[Quantity[pars[[2]],"Nanometers"],"Meters"]
},
2*UnitConvert[Abs[pm1*dso-d0*Quantity[ez,"Electronvolts"/"Angstroms"]],"Millielectronvolts"]
]


(* ::Input::Initialization:: *)
Lc[epscav_,N_][Ephot_]:=UnitConvert[(N*Pi*HB*CC)/(Ephot*Sqrt[epscav]),"Micrometers"]


(* ::Input::Initialization:: *)
Ldbr[epscav_,n1_,n2_][Ephotc_]:=UnitConvert[((HH*CC/Ephotc)*n1*n2)/(2*Sqrt[epscav]*Abs[n2-n1]),"Micrometers"]


(* ::Input::Initialization:: *)
Leff[epscav_,n1_,n2_,ndbr_,m_][Ephot_]:=Lc[epscav,m][Ephot]+ndbr*Ldbr[epscav,n1,n2][Ephot]


(* ::Input::Initialization:: *)
Mcav[epscav_][Ephot_]:=(Ephot*epscav)/CC^2


(* ::Input::Initialization:: *)
Mex[prc_,type_,ep_]:=4*FPmu[prc,type,ep]


(* ::Input::Initialization:: *)
HopXsq[Ex_,Ec_,V_]:=1/2 (1-(Ex-Ec)/Sqrt[(Ex-Ec)^2+(2*HB*V)^2])


(* ::Input::Initialization:: *)
HopCsq[Ex_,Ec_,V_]:=1/2 (1+(Ex-Ec)/Sqrt[(Ex-Ec)^2+(2*HB*V)^2])


(* ::Input::Initialization:: *)
gamcav[epscav_,R_,Leff_]:=UnitConvert[((1-Sqrt[R])/Sqrt[R])*(CC/(Sqrt[epscav]*Leff)),"Seconds"^-1]


gamcavreff[epscav_,r1_,r2_,leff_]:=UnitConvert[((1/2)*((1-Sqrt[r1])/Sqrt[r1]+(1-Sqrt[r2])/Sqrt[r2]))*(CC/(Sqrt[epscav]*leff)),"Seconds"^-1]


(* ::Input::Initialization:: *)
gampol[gamex_,gamcav_,HX_,HC_]:={HX*gamex +HC*gamcav,HC*gamex+HX*gamcav}


(* ::Input::Initialization:: *)
Exk[prc_,type_,ep_][k_]:=UnitConvert[Exres[prc,type,ep]+(HB*Quantity[k,"Micrometers"^(-1)])^2/(2*Mex[prc,type,ep]),"Millielectronvolts"]


(* ::Input::Initialization:: *)
Eck[epscav_][Ephot_][k_]:=UnitConvert[Ephot+(HB*Quantity[k,"Micrometers"^(-1)])^2/(2*Mcav[epscav][Ephot]),"Millielectronvolts"]


(* ::Input::Initialization:: *)
Ephotpct[prc_,type_,ep_][dE_]:=Exres[prc,type,ep]*(1+dE)


(* ::Input::Initialization:: *)
Elpup[Ex_,Ec_,gamex_,gamcav_,V_]:=UnitConvert[{(Ex+Ec+I*HB*(gamex+gamcav))/2-Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)],(Ex+Ec+I*HB*(gamex+gamcav))/2+Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)],
2*Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)]},"Millielectronvolts"]


Elp[Ex_,Ec_,gamex_,gamcav_,V_]:=UnitConvert[(Ex+Ec+I*HB*(gamex+gamcav))/2-Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)],"Millielectronvolts"]


Eup[Ex_,Ec_,gamex_,gamcav_,V_]:=UnitConvert[(Ex+Ec+I*HB*(gamex+gamcav))/2+Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)],"Millielectronvolts"]


(* ::Input::Initialization:: *)
Mlpup[Mex_,Mcav_,HX_,HC_]:={(HX/Mex+HC/Mcav)^-1,(HX/Mcav+HC/Mex)^-1}


(* ::Input::Initialization:: *)
Exres[prc_,type_,ep_]:=UnitConvert[2*FPegap[prc,type,ep]-Abs[FPev[prc,type,1,ep,1,0]],"Millielectronvolts"]


(* ::Input::Initialization:: *)
Vsav[epscav_,n1_,n2_,Nharm_,rmir_][prc_,type_,ep_][Ephot_]:=Sqrt[UnitConvert[(1+Sqrt[rmir])/Sqrt[rmir] ((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*Exres[prc,type,ep]*Leff[epscav,n1,n2,Nharm][Ephot]))*(FPr0[prc,type,1,ep,1,0]^2),"Seconds"^-2]]


Vsavmod[epscav_,n1_,n2_,m_,reff_][prc_,type_,ep_][Ephot_]:=Sqrt[UnitConvert[(1+Sqrt[rmir])/Sqrt[rmir] ((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*Exres[prc,type,ep]*Leff[epscav,n1,n2,m][Ephot]))*(FPr0[prc,type,1,ep,1,0]^2),"Seconds"^-2]]


(* ::Input::Initialization:: *)
Vflat[epscav_,n1_,n2_,Nharm_,Ndbr_,rmir_][prc_,type_,ep_][conc2d_,gamex_,epsbg_][Ephot_]:=Sqrt@UnitConvert[(1+Sqrt[rmir])/Sqrt[rmir]*((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2)*conc2d)/(epscav*Exres[prc,type,ep]*Leff[epscav,n1,n2,Nharm][Ephot]))+(((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2)*conc2d)/(CC*Sqrt[epscav]*Exres[prc,type,ep]))*(gamex+gamcav[epscav,rmir,Leff[epscav,n1,n2,Nharm][Ephot]]))+((Exres[prc,type,ep]*Quantity[prc[[1]][[4]],"Nanometers"]*epsbg)/(2*HB*Lc[epscav,Nharm][Ephot]*epscav*Sqrt[rmir]))^2,"Seconds"^-2]


(* ::Input::Initialization:: *)
EfromT[temp_]:=UnitConvert[Quantity["BoltzmannConstant"]*Quantity[temp,"Kelvins"],"Millielectronvolts"]


(* ::Input::Initialization:: *)
TfromE[E_]:=UnitConvert[E/Quantity["BoltzmannConstant"],"Kelvins"]


(* ::Input::Initialization:: *)
U0[prc_,type_,ep_,eps_]:=6*KK*(EE^2)*FPbohrad[prc,type,1,ep,1,0]/eps


(* ::Input::Initialization:: *)
Cs[Ueff_,npol_,Mlp_]:=UnitConvert[Sqrt[Ueff*npol/Mlp],"Meters"/"Seconds"]


(* ::Input::Initialization:: *)
Tcpol[npol_,cs_,Mp_,s_]:=Module[{Tc0=1/KB*((Pi*HB^2*npol*cs^4*Mp)/(6*s*Zeta[3]))^(1/3)},{((1+Sqrt[32/27 ((Mp*KB*Tc0)/(Pi*HB^2*npol))^3+1])^(1/3)-(Sqrt[32/27 ((Mp*KB*Tc0)/(Pi*HB^2*npol))^3+1]-1)^(1/3)) Tc0/2^(1/3),Tc0}//UnitSimplify]


Tcpol2[npol_,cs_,Mp_,s_]:=Module[{Tc0=1/KB*((2*Pi*HB^2*npol*cs^4*Mp)/(3*s*Zeta[3]))^(1/3)},{((1+Sqrt[32/27 ((Mp*KB*Tc0)/(Pi*HB^2*npol))^3+1])^(1/3)-(Sqrt[32/27 ((Mp*KB*Tc0)/(Pi*HB^2*npol))^3+1]-1)^(1/3)) Tc0/2^(1/3),Tc0}//UnitSimplify]


(* ::Input::Initialization:: *)
critdensfromT[tc_,cs_,mp_,s_]:=(Quantity[FindRoot[QuantityMagnitude@Tcpol2[na,cs,mp,s]==tc,{na,10^11}][[1]]//Last,"Meters"^-2])


(* ::Input::Initialization:: *)
normDens[cs_,mp_,s_,T_]:=UnitConvert[(3*Zeta[3])/(2*Pi*HB^2)*(s*(KB*Quantity[T,"Kelvins"])^3)/((cs^4)*mp),"Meters"^-2]


Tcdirect[mlp_,nlp_,ueff_]:=NSolve[tc==(((Pi*(HB^2))/(2*KB*mlp))*(nlp-((3*Zeta[3])/(2*Pi*(HB^2))*(((KB*tc)^3)*mlp)/(ueff*nlp)^2))),tc,Reals][[1]][[1]]//Last


Tcpostsolve[mlp_,nlp_,ueff_]:=Module[
	{
		a=(3*Zeta[3])/(4*(nlp^2)*(ueff^2)),
		c=(Pi*(HB^2)*nlp)/(2*mlp),
		A
	},
	A=(Sqrt[(3*(a^3))*(4+(27*a*(c^2)))]-(9*(a^2)*c));
	UnitConvert[(1/KB)*(((2/3)^(1/3)*(A^(-1/3)))-((1/(18*a^3))^(1/3)*(A^(1/3)))),"Kelvins"]
]


ns[mlp_,nlp_,ueff_,t_]:=nlp-(((3*Zeta[3])/(2*Pi*(HB^2)))*(((KB*t)^3*mlp)/(ueff^2*nlp^2)))


tc0[mlp_,nlp_,ueff_]:=nlp/KB ((2*Pi*(HB^2)*(ueff^2))/(3*Zeta[3]*mlp))^(1/3)


ncrit[mlp_,ueff_,t_]:=(KB*t)*((3*Zeta[3]*mlp)/(2*Pi*(HB^2)*(ueff^2)))^(1/3)


MakeMC[epscav_, n1_, n2_, rmir_, Ndbr_, prc_] := Module[
  {
   exrestab = Table[
     Exres[prc, type, ep],
     {type, 2}, {ep, FPEdim[prc]}
     ],
   lcfunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[dE,
       Lc[epscav, Nharm][Ephotpct[prc, type, ep][dE]]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   ldbrfunc = Table[
     With[
      {type = typee, ep = epp},
      Function[dE,
       Ldbr[epscav, n1, n2,Ndbr][Ephotpct[prc, type, ep][dE]]
       ]
      ],
     {typee, 2}, {epp, FPEdim[prc]}
     ],
   lefffunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[dE,
       Leff[epscav, n1, n2, Nharm,Ndbr][Ephotpct[prc, type, ep][dE]]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   gamcavfunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[dE,
       gamcav[epscav, rmir,
        Leff[epscav, n1, n2, Nharm,Ndbr][Ephotpct[prc, type, ep][dE]]]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   phi0sqtab = Table[
     FPr0[prc, type, 1, ep, 1, 0]^2,
     {type, 2}, {ep, FPEdim[prc]}
     ],
   Vsavfunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       dE,
       Vsav[epscav, n1, n2,Ndbr,Nharm, rmir][prc, type, ep][
        Ephotpct[prc, type, ep][dE]]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   Vflatfunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE},
       Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][conc2d,
          gamex, epsbg][Ephotpct[prc, type, ep][dE]]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   HopVsfunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {dE, k},
       {
        HopXsq[
         Exk[prc, type, ep][k],
         Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
         Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
          Ephotpct[prc, type, ep][dE]]
         ],
        HopCsq[
         Exk[prc, type, ep][k],
         Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
         Vsav[epscav, n1, n2, Ndbr,Nharm, rmir][prc, type, ep][
          Ephotpct[prc, type, ep][dE]]
         ]
        }
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   HopVffunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE, k},
       {
        HopXsq[
         Exk[prc, type, ep][k],
         Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
         Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
           conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
         ],
        HopCsq[
         Exk[prc, type, ep][k],
         Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
         Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
           conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
         ]
        }
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   ElpupVsfunc = Table[
     With[
      {type = ttype, Nharm = Nharmm, ep = epp},
      Function[{gamex, dE, k},
       Elpup[
        Exk[prc, type, ep][k],
        Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
        gamex,
        gamcav[epscav, rmir, 
         Leff[epscav, n1, n2, Nharm][Ephotpct[prc, type, ep][dE]]],
        Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
         Ephotpct[prc, type, ep][dE]]
        ]
       ]
      ],
     {ttype, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   ElpupVffunc = Table[
     With[
      {type = ttype, Nharm = Nharmm, ep = epp},
      Function[{conc2d, gamex, epsbg, dE, k},
       Elpup[
        Exk[prc, type, ep][k],
        Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
        gamex,
        gamcav[epscav, rmir,
         Leff[epscav, n1, n2, Nharm][Ephotpct[prc, type, ep][dE]]],
        Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
          conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
        ]
       ]
      ],
     {ttype, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   MlpupVsfunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {dE, k},
       Mlpup[
        Mex[prc, type, ep],
        Mcav[epscav][Eck[epscav][Ephotpct[prc, type, ep][dE]][0]],
        HopXsq[
         Exk[prc, type, ep][k],
         Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
         Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
          Ephotpct[prc, type, ep][dE]]
         ],
        HopCsq[
         Exk[prc, type, ep][k],
         Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
         Vsav[epscav, n1, n2, Ndbr,Nharm, rmir][prc, type, ep][
          Ephotpct[prc, type, ep][dE]]
         ]
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   MlpupVffunc = Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE, k},
       Mlpup[
        Mex[prc, type, ep],
        Mcav[epscav][Eck[epscav][Ephotpct[prc, type, ep][dE]][0]],
        HopXsq[
         Exk[prc, type, ep][k],
         Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
         Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
           conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
         ],
        HopCsq[
         Exk[prc, type, ep][k],
         Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
         Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
           conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
         ]
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   U0tab = Table[
     (6*KK*(EE^2))/ prc[[1]][[6]]*FPbohrad[prc, type, 1, ep, 1, 0],
     {type, 2}, {ep, FPEdim[prc]}
     ],
   U0kav = 
    Table[6*(Abs[
        FPev[prc, type, 1, ep, 1, 0]])*(FPbohrad[prc, type, 1, ep, 1, 
         0]^2), {type, 2}, {ep, FPEdim[prc]}],
   CsOlegVsfunc = Table[
     With[
      {type = ttype, Nharm = Nharmm, ep = epp},
      Function[{npol, dE},
       Cs[
        ((6*KK*(EE^2))/ prc[[1]][[6]]*
           FPbohrad[prc, type, 1, ep, 1, 0])*
         (HopXsq[
            Exk[prc, type, ep][0],
            Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
            
            Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
             Ephotpct[prc, type, ep][dE]]
            ])^2,
        npol,
        Mlpup[
          Mex[prc, type, ep],
          Mcav[epscav][Eck[epscav][Ephotpct[prc, type, ep][dE]][0]],
          HopXsq[
           Exk[prc, type, ep][0],
           Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
           
           Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
            Ephotpct[prc, type, ep][dE]]
           ],
          HopCsq[
           Exk[prc, type, ep][0],
           Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
           
           Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
            Ephotpct[prc, type, ep][dE]]
           ]
          ][[1]]
        ]
       ]
      ],
     {ttype, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   CsOlegVffunc = Table[
     With[
      {type = ttype, Nharm = Nharmm, ep = epp},
      Function[{conc2d, gamex, epsbg, npol, dE},
       Cs[
        ((6*KK*(EE^2))/ prc[[1]][[6]]*
           FPbohrad[prc, type, 1, ep, 1, 0])*
         (HopXsq[
            Exk[prc, type, ep][0],
            Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
            
            Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
              conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
            ])^2,
        npol,
        Mlpup[
          Mex[prc, type, ep],
          Mcav[epscav][Eck[epscav][Ephotpct[prc, type, ep][dE]][0]],
          HopXsq[
           Exk[prc, type, ep][0],
           Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
           
           Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
             conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
           ],
          HopCsq[
           Exk[prc, type, ep][0],
           Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
           
           Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
             conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
           ]
          ][[1]]
        ]
       ]
      ],
     {ttype, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   CsKavVsfunc = Table[
     With[
      {type = ttype, Nharm = Nharmm, ep = epp},
      Function[{npol, dE},
       Cs[
        (6*(Abs[
             FPev[prc, type, 1, ep, 1, 0]])*(FPbohrad[prc, type, 1, 
              ep, 1, 0]^2))*
         (HopXsq[
            Exk[prc, type, ep][0],
            Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
            
            Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
             Ephotpct[prc, type, ep][dE]]
            ])^2,
        npol,
        Mlpup[
          Mex[prc, type, ep],
          Mcav[epscav][Eck[epscav][Ephotpct[prc, type, ep][dE]][0]],
          HopXsq[
           Exk[prc, type, ep][0],
           Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
           
           Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
            Ephotpct[prc, type, ep][dE]]
           ],
          HopCsq[
           Exk[prc, type, ep][0],
           Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
           
           Vsav[epscav, n1, n2,Ndbr, Nharm, rmir][prc, type, ep][
            Ephotpct[prc, type, ep][dE]]
           ]
          ][[1]]
        ]
       ]
      ],
     {ttype, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   CsKavVffunc = Table[
     With[
      {type = ttype, Nharm = Nharmm, ep = epp},
      Function[{conc2d, gamex, epsbg, npol, dE},
       Cs[
        (6*(Abs[
             FPev[prc, type, 1, ep, 1, 0]])*(FPbohrad[prc, type, 1, 
              ep, 1, 0]^2))*
         (HopXsq[
            Exk[prc, type, ep][0],
            Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
            
            Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
              conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
            ])^2,
        npol,
        Mlpup[
          Mex[prc, type, ep],
          Mcav[epscav][Eck[epscav][Ephotpct[prc, type, ep][dE]][0]],
          HopXsq[
           Exk[prc, type, ep][0],
           Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
           
           Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
             conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
           ],
          HopCsq[
           Exk[prc, type, ep][0],
           Eck[epscav][Ephotpct[prc, type, ep][dE]][0],
           
           Vflat[epscav, n1, n2, Nharm, Ndbr, rmir][prc, type, ep][
             conc2d, gamex, epsbg][Ephotpct[prc, type, ep][dE]]
           ]
          ][[1]]
        ]
       ]
      ],
     {ttype, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ]
   },
  Association[
   "Inputs" -> 
    Association[
     "\!\(\*SubscriptBox[\(\[Epsilon]\), \(cav\)]\)" -> epscav, 
     "n1" -> n1, "n2" -> n2, "Rmir" -> rmir, "Ndbr" -> Ndbr, 
     "prc" -> ToString[prc[[1]]]],
   "prc" -> prc,
   "NNEtr" -> Table[NoNegEtr[prc, t], {t, 2}],
   "Emax" -> FPEdim[prc],
   "Lcav" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[dE,
       lcfunc[[type]][[Nharm]][[ep]][dE]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "Ldbr" -> Table[
     With[
      {type = typee, ep = epp},
      Function[dE,
       ldbrfunc[[type]][[ep]][dE]
       ]
      ],
     {typee, 2}, {epp, FPEdim[prc]}
     ],
   "Leff" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[dE,
       lefffunc[[type]][[Nharm]][[ep]][dE]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "DBRband" -> Table[
     Module[
      {type = typee, ep = epp, bw, ec, emax, emin, lc, lmin, lmax, 
       lbw},
      ec = exrestab[[type]][[ep]];
      bw = 
       UnitConvert[ec*(4/Pi)*ArcSin[(n2 - n1)/(n1 + n2)], 
        "Millielectronvolts"];
      emax = ec + (bw/2);
      emin = ec - (bw/2);
      lc = 
       UnitConvert[Quantity["PlanckConstant"]*CC/(ec*Sqrt[epscav]), 
        "Micrometers"];
      lmin = 
       UnitConvert[Quantity["PlanckConstant"]*CC/(emin*Sqrt[epscav]), 
        "Micrometers"];
      lmax = 
       UnitConvert[Quantity["PlanckConstant"]*CC/(emax*Sqrt[epscav]), 
        "Micrometers"];
      lbw = lmax - lmin;
      {
       {emin, ec, emax, bw},
       {lmin, lc, lmax, lbw}
       }
      ],
     {typee, 2}, {epp, FPEdim[prc]}
     ],
   "gamcav" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[dE,
       gamcavfunc[[type]][[Nharm]][[ep]][dE]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "exres" -> exrestab,
   "Exk" -> Table[
     With[
      {type = typee, ep = epp},
      Function[k,
       Exk[prc, type, ep][k]
       ]
      ],
     {typee, 2}, {epp, FPEdim[prc]}
     ],
   "Eck" -> Table[
     With[
      {type = typee, ep = epp},
      Function[{dE, k}, Eck[epscav][Ephotpct[prc, type, ep][dE]][k]]
      ],
     {typee, 2}, {epp, FPEdim[prc]}
     ],
   "Vsav" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[dE,
       UnitConvert[HB*Vsavfunc[[type]][[Nharm]][[ep]][dE], 
        "Millielectronvolts"]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "Vflat" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE},
       UnitConvert[
        HB*Vflatfunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE],
         "Millielectronvolts"]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "HopVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {dE, k},
       HopVsfunc[[type]][[Nharm]][[ep]][dE, k]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "HopVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE, k},
       HopVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, k]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "ElpupVs" -> Table[
     With[
      {type = ttype, Nharm = Nharmm, ep = epp},
      Function[{gamex, dE, k},
       Elpup[
        Exk[prc, type, ep][k],
        Eck[epscav][Ephotpct[prc, type, ep][dE]][k],
        gamex,
        gamcavfunc[[type]][[Nharm]][[ep]][dE],
        Vsavfunc[[type]][[Nharm]][[ep]][dE, k]
        ]
       ]
      ],
     {ttype, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "ElpupVf" -> Table[
     With[
      {type = ttype, Nharm = Nharmm, ep = epp},
      Function[{conc2d, gamex, epsbg, dE, k},
       ElpupVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, 
        k]
       ]
      ],
     {ttype, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "MpolVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {dE, k},
       MlpupVsfunc[[type]][[Nharm]][[ep]][dE, k]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "MpolVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE, k},
       MlpupVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, 
        k]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "gampolVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {gamex, dE, k},
       gampol[
        gamex,
        gamcavfunc[[type]][[Nharm]][[ep]][dE, k],
        HopVsfunc[[type]][[Nharm]][[ep]][dE, k][[1]],
        HopVsfunc[[type]][[Nharm]][[ep]][dE, k][[2]]
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "gampolVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE, k},
       gampol[
        gamex,
        gamcavfunc[[type]][[Nharm]][[ep]][dE, k],
        HopVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, 
          k][[1]],
        HopVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, 
          k][[2]]
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "UeffOlegVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       dE,
       UnitConvert[
        U0tab[[type]][[
          ep]]*(HopVsfunc[[type]][[Nharm]][[ep]][dE, 0][[1]])^2, 
        "Millielectronvolts"*"BohrRadius"^2]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "UeffOlegVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE},
       UnitConvert[
        U0tab[[type]][[
          ep]]*(HopVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg,
              dE, 0][[1]])^2, "Millielectronvolts"*"BohrRadius"^2]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "UeffKavVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       dE,
       UnitConvert[
        U0kav[[type]][[
          ep]]*(HopVsfunc[[type]][[Nharm]][[ep]][dE, 0][[1]])^2, 
        "Millielectronvolts"*"BohrRadius"^2]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "UeffKavVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, dE},
       UnitConvert[
        U0kav[[type]][[
          ep]]*(HopVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg,
              dE, 0][[1]])^2, "Millielectronvolts"*"BohrRadius"^2]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "CsOlegVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {npol, dE},
       CsOlegVsfunc[[type]][[Nharm]][[ep]][npol, dE]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "CsOlegVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, npol, dE},
       CsOlegVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, npol,
         dE]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "CsKavVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {npol, dE},
       CsKavVsfunc[[type]][[Nharm]][[ep]][npol, dE]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "CsKavVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, npol, dE},
       CsKavVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, npol, 
        dE]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "TcOlegVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {npol, dE},
       Tcpol2[
        npol,
        CsOlegVsfunc[[type]][[Nharm]][[ep]][npol, dE],
        MlpupVsfunc[[type]][[Nharm]][[ep]][dE, 0][[1]],
        1
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "TcOlegVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, npol, dE},
       Tcpol2[
        npol,
        CsOlegVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, 
         npol, dE],
        MlpupVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, 
          0][[1]],
        1
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "TcKavVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {npol, dE},
       Tcpol2[
        npol,
        CsKavVsfunc[[type]][[Nharm]][[ep]][npol, dE],
        MlpupVsfunc[[type]][[Nharm]][[ep]][dE, 0][[1]],
        1
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "TcKVmax"->Function[
     {type,m,ep,npol,de},
     With[
     {
     mph=Mcav[epscav][Ephotpct[prc,type,ep][de]],
     mex=4*FPmu[prc,type,ep],
     xsq=HopXsq[Exres[prc,type,ep],Ephotpct[prc,type,ep][de],Vsav[epscav,n1,n2,m,rmir][prc,type,ep][Ephotpct[prc,type,ep][de]]]
     },
     tc0[Mlpup[mex,mph,xsq,1-xsq][[1]], npol, 6*(Exres[prc,type,ep])*(FPbohrad[prc,type,1,ep,1,0]^2)*(xsq^2)]
     ]
     ],
   "TcKVmax2"->Function[
     {type,m,ep,npol,de},
     With[
     {
     mph=Mcav[epscav][Ephotpct[prc,type,ep][de]],
     mex=4*FPmu[prc,type,ep],
     xsq=HopXsq[
         Exk[prc, type, ep][0],
         Eck[epscav][Ephotpct[prc, type, ep][de]][0],
         Vsav[epscav, n1, n2, m, rmir][prc, type, ep][
          Ephotpct[prc, type, ep][de]]
         ]
     },
       Tcpol2[
          npol,
          Sqrt[(6*(Exres[prc,type,ep])*(FPbohrad[prc,type,1,ep,1,0]^2)*(xsq^2))*npol/(Mlpup[mex,mph,xsq,1-xsq][[1]])],
          Mlpup[mex,mph,xsq,1-xsq][[1]],
          1
       ]
      ]
     ],
   "TcKavVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, npol, dE},
       Tcpol2[
        npol,
        CsKavVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, npol,
          dE],
        MlpupVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, 
          0][[1]],
        1
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "CritnOlegVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {Tc, dE},
       (Quantity[FindRoot[QuantityMagnitude@Tcpol2[
               Quantity[na, "Meters"^-2],
               
               CsOlegVsfunc[[type]][[Nharm]][[ep]][
                Quantity[na, "Meters"^-2], dE],
               MlpupVsfunc[[type]][[Nharm]][[ep]][dE, 0][[1]],
               1][[1]] == Tc,
            {na, 10^13}][[1]] // Last, "Meters"^-2])
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "CritnOlegVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, Tc, dE},
       (Quantity[FindRoot[QuantityMagnitude@Tcpol2[
               Quantity[na, "Meters"^-2],
               
               CsOlegVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, 
                epsbg, Quantity[na, "Meters"^-2], dE],
               
               MlpupVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, 
                 epsbg, dE, 0][[1]],
               1][[1]] == Tc,
            {na, 10^13}][[1]] // Last, "Meters"^-2])
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "CritnKavVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {Tc, dE},
       (Quantity[FindRoot[QuantityMagnitude@Tcpol2[
               Quantity[na, "Meters"^-2],
               
               CsKavVsfunc[[type]][[Nharm]][[ep]][
                Quantity[na, "Meters"^-2], dE],
               MlpupVsfunc[[type]][[Nharm]][[ep]][dE, 0][[1]],
               1][[1]] == Tc,
            {na, 10^13}][[1]] // Last, "Meters"^-2])
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "CritnKavVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, Tc, dE},
       (Quantity[FindRoot[QuantityMagnitude@Tcpol2[
               Quantity[na, "Meters"^-2],
               
               CsKavVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, 
                epsbg, Quantity[na, "Meters"^-2], dE],
               
               MlpupVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, 
                 epsbg, dE, 0][[1]],
               1][[1]] == Tc,
            {na, 10^13}][[1]] // Last, "Meters"^-2])
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "normdensOlegVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {na, dE, T},
       normDens[
        CsOlegVsfunc[[type]][[Nharm]][[ep]][na, dE],
        MlpupVsfunc[[type]][[Nharm]][[ep]][dE, 0][[1]],
        1,
        T
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "normdensOlegVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, na, dE, T},
       normDens[
        CsOlegVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, na, 
         dE],
        MlpupVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, 
          0][[1]],
        1,
        T
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "normdensKavVs" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {na, dE, T},
       normDens[
        CsKavVsfunc[[type]][[Nharm]][[ep]][na, dE],
        MlpupVsfunc[[type]][[Nharm]][[ep]][dE, 0][[1]],
        1,
        T
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ],
   "normdensKavVf" -> Table[
     With[
      {type = typee, Nharm = Nharmm, ep = epp},
      Function[
       {conc2d, gamex, epsbg, na, dE, T},
       normDens[
        CsKavVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, na, 
         dE],
        MlpupVffunc[[type]][[Nharm]][[ep]][conc2d, gamex, epsbg, dE, 
          0][[1]],
        1,
        T
        ]
       ]
      ],
     {typee, 2}, {Nharmm, 10}, {epp, FPEdim[prc]}
     ]
   ]
  ]


MCall[epscav_, n1_, n2_, r1_, r2_, prc_, genname_] := Module[
  {
   exrestab = Table[
     Exres[prc, type, ep],
     {type, 2}, {ep, FPEdim[prc]}
     ],
   Mex = Function[{type, ep}, 4*FPmu[prc, type, ep]],
   Mph = Function[{type, ep, de}, 
     Ephotpct[prc, type, ep][de]*epscav/(CC^2)],
   Reff = (4*r1*r2)/((Sqrt[r1] + Sqrt[r2])^2),
   Leff1 = 
    Function[{type, ep, m, de}, 
     UnitConvert[Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
      Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]],"Micrometers"],
   Leff2 = 
    Function[{type, ep, m, de}, 
     UnitConvert[Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
      2*Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]],"Micrometers"],
   MOgamcav = 
    Function[{type, ep, m, de},
    UnitConvert[(1 - Sqrt[r1])/
      Sqrt[r1]*(CC/(Sqrt[
           epscav]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
            Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))),"Seconds"^-1]],
   MNgamcav = 
    Function[{type, ep, m, de},
    UnitConvert[(1 - Sqrt[r1])/
      Sqrt[r1]*(CC/(Sqrt[
           epscav]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
            2*Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))),"Seconds"^-1]],
   OOgamcav = 
    Function[{type, ep, m, de}, 
    UnitConvert[
    (1 - Sqrt[r2])/Sqrt[r2]*(CC/(2*Sqrt[epscav]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
            Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))),"Seconds"^-1]],
   ONgamcav = 
    Function[{type, ep, m, 
      de}, 
      UnitConvert[(1/2*((1 - Sqrt[r1])/Sqrt[r1] + (1 - Sqrt[r2])/
          Sqrt[r2]))*(CC/(Sqrt[
           epscav]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
            Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))),"Seconds"^-1]],
   MOV = Function[{type, ep, m, de}, 
     Sqrt[UnitConvert[((1 + Sqrt[r1])/
          Sqrt[r1])*((4*Pi*
            KK*(EE^2)*(Quantity[prc[[1]][[3]], 
               "Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*
            Exres[prc, type, 
             ep]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
              Ldbr[epscav, n1, n2][
               Ephotpct[prc, type, ep][0]])))*(FPr0[prc, type, 1, ep, 
           1, 0]^2), "Seconds"^-2]]],
   MNV = Function[{type, ep, m, de}, 
     Sqrt[UnitConvert[((1 + Sqrt[r1])/
          Sqrt[r1])* ((4*Pi*
            KK*(EE^2)*(Quantity[prc[[1]][[3]], 
               "Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*
            Exres[prc, type, 
             ep]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
              2*Ldbr[epscav, n1, n2][
                Ephotpct[prc, type, ep][0]])))*(FPr0[prc, type, 1, ep,
            1, 0]^2), "Seconds"^-2]]],
   OOV = Function[{type, ep, m, de}, 
     Sqrt[UnitConvert[((1 + Sqrt[r2])/
          Sqrt[r2])* ((4*Pi*
            KK*(EE^2)*(Quantity[prc[[1]][[3]], 
               "Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*
            Exres[prc, type, 
             ep]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
              Ldbr[epscav, n1, n2][
               Ephotpct[prc, type, ep][0]])))*(FPr0[prc, type, 1, ep, 
           1, 0]^2), "Seconds"^-2]]],
   ONV = Function[{type, ep, m, de}, 
     Sqrt[UnitConvert[((1 + 
            Sqrt[(4*r1*r2)/((Sqrt[r1] + Sqrt[r2])^2)])/
          Sqrt[(4*r1*r2)/((Sqrt[r1] + Sqrt[r2])^2)])* ((4*Pi*
            KK*(EE^2)*(Quantity[prc[[1]][[3]], 
               "Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*
            Exres[prc, type, 
             ep]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 
              Ldbr[epscav, n1, n2][
               Ephotpct[prc, type, ep][0]])))*(FPr0[prc, type, 1, ep, 
           1, 0]^2), "Seconds"^-2]]],
   U0 = Function[{type, ep}, 
     6*(Abs@FPev[prc, type, 1, ep, 1, 0])*(FPbohrad[prc, type, 1, ep, 
         1, 0]^2)],
   HXsq = Function[{V, type, ep, m, de, k},
     HopXsq[
      Exres[prc, type, ep],
      Ephotpct[prc, type, ep][de],
      V[type, ep, m, de]
      ]
     ],
   HCsq = Function[{V, type, ep, m, de, k},
     HopCsq[
      Exres[prc, type, ep],
      Ephotpct[prc, type, ep][de],
      V[type, ep, m, de]
      ]
     ],
   MLP = Function[{V, type, ep, m, de, k},
     ((HopXsq[Exres[prc, type, ep], Ephotpct[prc, type, ep][de], 
          V[type, ep, m, de]]/(4*FPmu[prc, type, ep]) + 
         HopCsq[Exres[prc, type, ep], Ephotpct[prc, type, ep][de], 
          V[type, ep, m, de]]/(
         Ephotpct[prc, type, ep][de]*epscav/(CC^2)))^(-1))
     ]
   },
  Association[
   "inps" -> <|"epscav" -> epscav, "n1" -> n1, "n2" -> n2|>,
   "genname" -> ToString@genname,
   "prc" -> prc,
   "NNEtr" -> Function[{t}, NoNegEtr[prc, t]],
   "Emax" -> FPEdim[prc],
   "exres" -> Function[{type, ep}, Exres[prc, type, ep]],
   "lc" -> 
    Function[{type, ep, m, de}, 
     Lc[epscav, m][Ephotpct[prc, type, ep][de]]],
   "ldbr" -> 
    Function[{type, ep}, 
     Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]],
   "dbrband" -> Function[{type, ep}, With[
      {ec = Exres[prc, type, ep], 
       bw = Exres[prc, type, ep]*(4/Pi)*ArcSin[(n2 - n1)/(n1 + n2)]},
      UnitConvert[{ec, ec + (bw/2), ec - (bw/2)}, "Millielectronvolts"]
      ]
     ],
   "reff" -> Reff,
   "mph" -> Mph,
   "mex" -> Mex,
   "Exk" -> 
    Function[{type, ep, k}, 
     Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(
      2*Mex[type, ep])
     ],
   "Eck" -> Function[{type, ep, de, k},
     Ephotpct[prc, type, ep][de] + ((HB^2)*
       Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de])
     ],
   "mono" -> Association[
     "confname" -> genname + ", 2 DBR",
     "old" -> Association[
       "methname" -> genname + ", 2 DBR, Standard approach",
       "leff" -> Function[{type, ep, m, de}, Leff1[type, ep, m, de]],
       "gamcav" -> 
        Function[{type, ep, m, de}, MOgamcav[type, ep, m, de]],
       "V" -> Function[{type, ep, m, de}, MOV[type, ep, m, de]],
       "HX" -> 
        Function[{type, ep, m, de, k}, 
         HXsq[MOV, type, ep, m, de, k]],
       "HC" -> 
        Function[{type, ep, m, de, k}, 
         HCsq[MOV, type, ep, m, de, k]],
       "mlp" -> 
        Function[{type, ep, m, de, 
          k}, (HXsq[MOV, type, ep, m, de, k]/Mex[type, ep] + 
            HCsq[MOV, type, ep, m, de, k]/Mph[type, ep, de])^(-1)
         ],
       "mup" -> 
        Function[{type, ep, m, de, 
          k}, (HCsq[MOV, type, ep, m, de, k]/Mex[type, ep] + 
            HXsq[MOV, type, ep, m, de, k]/Mph[type, ep, de])^(-1) 
         ],
       "Elp" -> Function[{type, ep, m, gamex, de, k},
         Elp[
          
          Exres[prc, type, ep] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          
          Ephotpct[prc, type, ep][de] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          MOgamcav[type, ep, m, de],
          MOV[type, ep, m, de]
          ]
         ],
       "Eup" -> Function[{type, ep, m, gamex, de, k},
         Eup[
          
          Exres[prc, type, ep] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          
          Ephotpct[prc, type, ep][de] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          MOgamcav[type, ep, m, de],
          MOV[type, ep, m, de]
          ]
         ],
       "Elpup" -> Function[{type, ep, m, gamex, de, k},
         Eup[
           
           Exres[prc, type, ep] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           
           Ephotpct[prc, type, ep][de] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           MOgamcav[type, ep, m, de],
           MOV[type, ep, m, de]
           ] - Elp[
           
           Exres[prc, type, ep] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           
           Ephotpct[prc, type, ep][de] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           MOgamcav[type, ep, m, de],
           MOV[type, ep, m, de]
           ]
         ],
       "gamlp" -> Function[{type, ep, m, gamex, de, k},
         (HXsq[MOV, type, ep, m, de, k]*
            gamex) + (HCsq[MOV, type, ep, m, de, k]*
            MOgamcav[type, ep, m, de])
         ],
       "gamup" -> Function[{type, ep, m, gamex, de, k},
         (HCsq[MOV, type, ep, m, de, k]*
            gamex) + (HXsq[MOV, type, ep, m, de, k]*
            MOgamcav[type, ep, m, de])
         ],
       "Ueff" -> 
        Function[{type, ep, m, de}, 
         U0[type, ep]*(HXsq[MOV, type, ep, m, de, 0]^2)],
       "cs" -> 
        Function[{type, ep, m, npol, de}, 
         Sqrt[
          U0[type, ep]*(HXsq[MOV, type, ep, m, de, 0]^2)*
           npol/(MLP[MOV, type, ep, m, de, 0])]],
       "tc0" -> Function[{type, ep, m, npol, de},
         tc0[
          MLP[MOV, type, ep, m, de, 0],
          npol,
          U0[type, ep]*(HXsq[MOV, type, ep, m, de, 0]^2)
          ]
         ],
       "ncrit" -> Function[{type, ep, m, t, de},
         ncrit[
          MLP[MOV, type, ep, m, de, 0],
          U0[type, ep]*(HXsq[MOV, type, ep, m, de, 0]^2),
          t
          ]
         ],
       "tcbkt" -> Function[{type, ep, m, npol, de},
         Tcpostsolve[
          MLP[MOV, type, ep, m, de, 0],
          npol,
          U0[type, ep]*(HXsq[MOV, type, ep, m, de, 0]^2)
          ]
         ],
       "tcmax" -> Function[{t, ep, m, conc},
         Module[{max},
          max = FindMaximum[
            QuantityMagnitude[
             (Tcpol2[
                conc,
                
                Sqrt[U0[t, ep]*(HXsq[MOV, t, ep, m, dE, 0]^2)*
                  conc/(MLP[MOV, t, ep, m, dE, 0])], 
                MLP[MOV, t, ep, m, dE, 0],
                1
                ])[[1]]],
            {dE, 0.01}
            ];
          {dE /. Last@max, max[[1]]}
          ]
         ]
       ],
     "new" -> Association[
       "methname" -> genname + ", 2 DBR, Modified",
       "leff" -> Function[{type, ep, m, de}, Leff2[type, ep, m, de]],
       "gamcav" -> 
        Function[{type, ep, m, de}, MNgamcav[type, ep, m, de]],
       "V" -> Function[{type, ep, m, de}, MNV[type, ep, m, de]],
       "HX" -> 
        Function[{type, ep, m, de, k}, 
         HXsq[MNV, type, ep, m, de, k]],
       "HC" -> 
        Function[{type, ep, m, de, k}, 
         HCsq[MNV, type, ep, m, de, k]],
       "mlp" -> 
        Function[{type, ep, m, de, 
          k}, (HXsq[MNV, type, ep, m, de, k]/Mex[type, ep] + 
            HCsq[MNV, type, ep, m, de, k]/Mph[type, ep, de])^(-1)
         ],
       "mup" -> 
        Function[{type, ep, m, de, 
          k}, (HCsq[MNV, type, ep, m, de, k]/Mex[type, ep] + 
            HXsq[MNV, type, ep, m, de, k]/Mph[type, ep, de])^(-1) 
         ],
       "Elp" -> Function[{type, ep, m, gamex, de, k},
         Elp[
          
          Exres[prc, type, ep] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          
          Ephotpct[prc, type, ep][de] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          MNgamcav[type, ep, m, de],
          MNV[type, ep, m, de]
          ]
         ],
       "Eup" -> Function[{type, ep, m, gamex, de, k},
         Eup[
          
          Exres[prc, type, ep] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          
          Ephotpct[prc, type, ep][de] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          MNgamcav[type, ep, m, de],
          MNV[type, ep, m, de]
          ]
         ],
       "Elpup" -> Function[{type, ep, m, gamex, de, k},
         Eup[
           
           Exres[prc, type, ep] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           
           Ephotpct[prc, type, ep][de] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           MNgamcav[type, ep, m, de],
           MNV[type, ep, m, de]
           ] - Elp[
           
           Exres[prc, type, ep] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           
           Ephotpct[prc, type, ep][de] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           MNgamcav[type, ep, m, de],
           MNV[type, ep, m, de]
           ]
         ],
       "gamlp" -> Function[{type, ep, m, gamex, de, k},
         (HXsq[MNV, type, ep, m, de, k]*
            gamex) + (HCsq[MNV, type, ep, m, de, k]*
            MNgamcav[type, ep, m, de])
         ],
       "gamup" -> Function[{type, ep, m, gamex, de, k},
         (HCsq[MNV, type, ep, m, de, k]*
            gamex) + (HXsq[MNV, type, ep, m, de, k]*
            MNgamcav[type, ep, m, de])
         ],
       "Ueff" -> 
        Function[{type, ep, m, de}, 
         U0[type, ep]*(HXsq[MNV, type, ep, m, de, 0]^2)],
       "cs" -> 
        Function[{type, ep, m, npol, de}, 
         Sqrt[U0[type, ep]*(HXsq[MNV, type, ep, m, de, 0]^2)*
           npol/(MLP[MNV, type, ep, m, de, 0])]],
       "tc0" -> Function[{type, ep, m, npol, de},
         tc0[
          MLP[MNV, type, ep, m, de, 0],
          npol,
          U0[type, ep]*(HXsq[MNV, type, ep, m, de, 0]^2)
          ]
         ],
       "ncrit" -> Function[{type, ep, m, t, de},
         ncrit[
          MLP[MNV, type, ep, m, de, 0],
          U0[type, ep]*(HXsq[MNV, type, ep, m, de, 0]^2),
          t
          ]
         ],
       "tcbkt" -> Function[{type, ep, m, npol, de},
         Tcpostsolve[
          MLP[MNV, type, ep, m, de, 0],
          npol,
          U0[type, ep]*(HXsq[MNV, type, ep, m, de, 0]^2)
          ]
         ],
       "tcmax" -> Function[{t, ep, m, conc},
         Module[{max},
          max = FindMaximum[
            QuantityMagnitude[
             (Tcpol2[
                conc,
                
                Sqrt[U0[t, ep]*(HXsq[MNV, t, ep, m, dE, 0]^2)*
                  conc/(MLP[MNV, t, ep, m, dE, 0])], 
                MLP[MNV, t, ep, m, dE, 0],
                1
                ])[[1]]],
            {dE, 0.01}
            ];
          {dE /. Last@max, max[[1]]}
          ]
         ]
       ]
     ],
   "open" -> Association[
     "confname" -> genname + ", 1 DBR",
     "old" -> Association[
       "methname" -> genname + ", 1 DBR, Standard",
       "leff" -> Function[{type, ep, m, de}, Leff1[type, ep, m, de]],
       "gamcav" -> 
        Function[{type, ep, m, de}, OOgamcav[type, ep, m, de]],
       "V" -> Function[{type, ep, m, de}, OOV[type, ep, m, de]],
       "HX" -> 
        Function[{type, ep, m, de, k}, 
         HXsq[OOV, type, ep, m, de, k]],
       "HC" -> 
        Function[{type, ep, m, de, k}, 
         HCsq[OOV, type, ep, m, de, k]],
       "mlp" -> 
        Function[{type, ep, m, de, 
          k}, (HXsq[OOV, type, ep, m, de, k]/Mex[type, ep] + 
            HCsq[OOV, type, ep, m, de, k]/Mph[type, ep, de])^(-1)
         ],
       "mup" -> 
        Function[{type, ep, m, de, 
          k}, (HCsq[OOV, type, ep, m, de, k]/Mex[type, ep] + 
            HXsq[OOV, type, ep, m, de, k]/Mph[type, ep, de])^(-1) 
         ],
       "Elp" -> Function[{type, ep, m, gamex, de, k},
         Elp[
          
          Exres[prc, type, ep] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          
          Ephotpct[prc, type, ep][de] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          OOgamcav[type, ep, m, de],
          OOV[type, ep, m, de]
          ]
         ],
       "Eup" -> Function[{type, ep, m, gamex, de, k},
         Eup[
          
          Exres[prc, type, ep] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          
          Ephotpct[prc, type, ep][de] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          OOgamcav[type, ep, m, de],
          OOV[type, ep, m, de]
          ]
         ],
       "Elpup" -> Function[{type, ep, m, gamex, de, k},
         Eup[
           
           Exres[prc, type, ep] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           
           Ephotpct[prc, type, ep][de] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           OOgamcav[type, ep, m, de],
           OOV[type, ep, m, de]
           ] - Elp[
           
           Exres[prc, type, ep] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           
           Ephotpct[prc, type, ep][de] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           OOgamcav[type, ep, m, de],
           OOV[type, ep, m, de]
           ]
         ],
       "gamlp" -> Function[{type, ep, m, gamex, de, k},
         (HXsq[OOV, type, ep, m, de, k]*
            gamex) + (HCsq[OOV, type, ep, m, de, k]*
            OOgamcav[type, ep, m, de])
         ],
       "gamup" -> Function[{type, ep, m, gamex, de, k},
         (HCsq[OOV, type, ep, m, de, k]*
            gamex) + (HXsq[OOV, type, ep, m, de, k]*
            OOgamcav[type, ep, m, de])
         ],
       "Ueff" -> 
        Function[{type, ep, m, de}, 
         U0[type, ep]*(HXsq[OOV, type, ep, m, de, 0]^2)],
       "cs" -> 
        Function[{type, ep, m, npol, de}, 
         Sqrt[U0[type, ep]*(HXsq[OOV, type, ep, m, de, 0]^2)*
           npol/(MLP[OOV, type, ep, m, de, 0])]],
       "tc0" -> Function[{type, ep, m, npol, de},
         tc0[
          MLP[OOV, type, ep, m, de, 0],
          npol,
          U0[type, ep]*(HXsq[OOV, type, ep, m, de, 0]^2)
          ]
         ],
       "ncrit" -> Function[{type, ep, m, t, de},
         ncrit[
          MLP[OOV, type, ep, m, de, 0],
          U0[type, ep]*(HXsq[OOV, type, ep, m, de, 0]^2),
          t
          ]
         ],
       "tcbkt" -> Function[{type, ep, m, npol, de},
         Tcpostsolve[
          MLP[OOV, type, ep, m, de, 0],
          npol,
          U0[type, ep]*(HXsq[OOV, type, ep, m, de, 0]^2)
          ]
         ],
       "tcmax" -> Function[{t, ep, m, conc},
         Module[{max},
          max = FindMaximum[
            QuantityMagnitude[
             (Tcpol2[
                conc,
                
                Sqrt[U0[t, ep]*(HXsq[OOV, t, ep, m, dE, 0]^2)*
                  conc/(MLP[OOV, t, ep, m, dE, 0])], 
                MLP[OOV, t, ep, m, dE, 0],
                1
                ])[[1]]],
            {dE, 0.01}
            ];
          {dE /. Last@max, max[[1]]}
          ]
         ]
       ],
     "new" -> Association[
       "methname" -> genname + ", 1 DBR, Modified",
       "leff" -> Function[{type, ep, m, de}, Leff1[type, ep, m, de]],
       "gamcav" -> 
        Function[{type, ep, m, de}, ONgamcav[type, ep, m, de]],
       "V" -> Function[{type, ep, m, de}, ONV[type, ep, m, de]],
       "HX" -> 
        Function[{type, ep, m, de, k}, 
         HXsq[ONV, type, ep, m, de, k]],
       "HC" -> 
        Function[{type, ep, m, de, k}, 
         HCsq[ONV, type, ep, m, de, k]],
       "mlp" -> 
        Function[{type, ep, m, de, 
          k}, (HXsq[ONV, type, ep, m, de, k]/Mex[type, ep] + 
            HCsq[ONV, type, ep, m, de, k]/Mph[type, ep, de])^(-1)
         ],
       "mup" -> 
        Function[{type, ep, m, de, 
          k}, (HCsq[ONV, type, ep, m, de, k]/Mex[type, ep] + 
            HXsq[ONV, type, ep, m, de, k]/Mph[type, ep, de])^(-1) 
         ],
       "Elp" -> Function[{type, ep, m, gamex, de, k},
         Elp[
          
          Exres[prc, type, ep] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          
          Ephotpct[prc, type, ep][de] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          ONgamcav[type, ep, m, de],
          ONV[type, ep, m, de]
          ]
         ],
       "Eup" -> Function[{type, ep, m, gamex, de, k},
         Eup[
          
          Exres[prc, type, ep] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          
          Ephotpct[prc, type, ep][de] + ((HB^2)*
            Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          ONgamcav[type, ep, m, de],
          ONV[type, ep, m, de]
          ]
         ],
       "Elpup" -> Function[{type, ep, m, gamex, de, k},
         Eup[
           
           Exres[prc, type, ep] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           
           Ephotpct[prc, type, ep][de] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           ONgamcav[type, ep, m, de],
           ONV[type, ep, m, de]
           ] - Elp[
           
           Exres[prc, type, ep] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           
           Ephotpct[prc, type, ep][de] + ((HB^2)*
             Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           ONgamcav[type, ep, m, de],
           ONV[type, ep, m, de]
           ]
         ],
       "gamlp" -> Function[{type, ep, m, gamex, de, k},
         (HXsq[ONV, type, ep, m, de, k]*
            gamex) + (HCsq[ONV, type, ep, m, de, k]*
            ONgamcav[type, ep, m, de])
         ],
       "gamup" -> Function[{type, ep, m, gamex, de, k},
         (HCsq[ONV, type, ep, m, de, k]*
            gamex) + (HXsq[ONV, type, ep, m, de, k]*
            ONgamcav[type, ep, m, de])
         ],
       "Ueff" -> 
        Function[{type, ep, m, de}, 
         U0[type, ep]*(HXsq[ONV, type, ep, m, de, 0]^2)],
       "cs" -> 
        Function[{type, ep, m, npol, de}, 
         Sqrt[U0[type, ep]*(HXsq[ONV, type, ep, m, de, 0]^2)*
           npol/(MLP[ONV, type, ep, m, de, 0])]],
       "tc0" -> Function[{type, ep, m, npol, de},
         tc0[
          MLP[ONV, type, ep, m, de, 0],
          npol,
          U0[type, ep]*(HXsq[ONV, type, ep, m, de, 0]^2)
          ]
         ],
       "ncrit" -> Function[{type, ep, m, t, de},
         ncrit[
          MLP[ONV, type, ep, m, de, 0],
          U0[type, ep]*(HXsq[ONV, type, ep, m, de, 0]^2),
          t
          ]
         ],
       "tcbkt" -> Function[{type, ep, m, npol, de},
         Tcpostsolve[
          MLP[ONV, type, ep, m, de, 0],
          npol,
          U0[type, ep]*(HXsq[ONV, type, ep, m, de, 0]^2)
          ]
         ],
       "tcmax" -> Function[{t, ep, m, conc},
         Module[{max},
          max = FindMaximum[
            QuantityMagnitude[
             (Tcpol2[
                conc,
                
                Sqrt[U0[t, ep]*(HXsq[ONV, t, ep, m, dE, 0]^2)*
                  conc/(MLP[ONV, t, ep, m, dE, 0])], 
                MLP[ONV, t, ep, m, dE, 0],
                1
                ])[[1]]],
            {dE, 0.01}
            ];
          {dE /. Last@max, max[[1]]}
          ]
         ]
       ]
     ]
   ]
  ]
