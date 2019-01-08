SetDirectory@NotebookDirectory[];
Import["f-filemine.m"];

bhprc=Import["results/big-sihbn-DIR-si-2018-05-27/proc.m"];
shprc=Import["results/small-sihbn-DIR-si-2018-05-27/proc.m"];
fssiprc=Import["results/fssi-DIR-si-2018-06-05/reproc.m"];
fsgeprc=Import["results/fsge-DIR-ge-2018-06-05/reproc.m"];
fssnprc=Import["results/fssn-DIR-sn-2018-06-05/reproc.m"];

HB=Quantity["ReducedPlanckConstant"];
HH=Quantity["PlanckConstant"];
CC=Quantity["SpeedOfLight"];
KK=1/(4*Pi*Quantity["ElectricConstant"]);
EE=Quantity["ElementaryCharge"];
KB=Quantity["BoltzmannConstant"];

configs={"mono","open"};
methods={"old","new"};

QConc[conc_]:=Quantity[conc,"Meters"^-2]
QFreq[freq_]:=Quantity[freq,"Seconds"^-1]
CTmev[quant_]:=UnitConvert[quant,"Millielectronvolts"]

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

nsoleffeps[prc_,t_,ep_]:=NSolve[-FPev[prc,t,1,ep,1,0]==(2*(Quantity["Coulomb Constant"]*Quantity["ElementaryCharge"]^2)^2*FPmu[prc,t,ep])/((Quantity["ReducedPlanckConstant"]^2)*(\[Epsilon]^2))&&\[Epsilon]>0,\[Epsilon]][[1]][[1]]//Last

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

allmatsAB={"FS Si, B","FS Si, A","FS Ge, B","FS Ge, A","FS Sn, B","FS Sn, A","Type I, B","Type I, A","Type II, B","Type II, A"};
allmatsA={"FS Si","FS Ge","FS Sn","Type I","Type II"};
eplabel="\!\(\*SubscriptBox[\(E\), \(\[Perpendicular]\)]\) [V/\[Angstrom]]";
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

plotehmass[pm1_,pars_][ez_]:=Module[
{
dso=UnitConvert[Quantity[pars[[1]],"Millielectronvolts"],"Joules"],
d0=UnitConvert[Quantity[pars[[2]],"Nanometers"],"Meters"],
vf=Quantity[pars[[3]],"Meters"/"Seconds"]
},
UnitConvert[Abs[(pm1*dso-d0*Quantity[ez,"Electronvolts"/"Angstroms"])]/(vf^2),"ElectronMass"]
]

plotgap[pm1_,pars_][ez_]:=Module[
{
dso=UnitConvert[Quantity[pars[[1]],"Millielectronvolts"],"Joules"],
d0=UnitConvert[Quantity[pars[[2]],"Nanometers"],"Meters"]
},
2*UnitConvert[Abs[pm1*dso-d0*Quantity[ez,"Electronvolts"/"Angstroms"]],"Millielectronvolts"]
]

Lc[epscav_,N_][Ephot_]:=UnitConvert[(N*Pi*HB*CC)/(Ephot*Sqrt[epscav]),"Micrometers"]

Ldbr[epscav_,n1_,n2_][Ephotc_]:=UnitConvert[((HH*CC/Ephotc)*n1*n2)/(2*Sqrt[epscav]*Abs[n2-n1]),"Micrometers"]

Leff[epscav_,n1_,n2_,ndbr_,m_][Ephot_]:=Lc[epscav,m][Ephot]+ndbr*Ldbr[epscav,n1,n2][Ephot]

Mcav[epscav_][Ephot_]:=(Ephot*epscav)/CC^2

Mex[prc_,type_,ep_]:=4*FPmu[prc,type,ep]

HopXsq[Ex_,Ec_,V_]:=1/2 (1-(Ex-Ec)/Sqrt[(Ex-Ec)^2+(2*HB*V)^2])

HopCsq[Ex_,Ec_,V_]:=1/2 (1+(Ex-Ec)/Sqrt[(Ex-Ec)^2+(2*HB*V)^2])

gamcav[epscav_,R_,Leff_]:=UnitConvert[((1-Sqrt[R])/Sqrt[R])*(CC/(Sqrt[epscav]*Leff)),"Seconds"^-1]

gamcavreff[epscav_,r1_,r2_,leff_]:=UnitConvert[((1/2)*((1-Sqrt[r1])/Sqrt[r1]+(1-Sqrt[r2])/Sqrt[r2]))*(CC/(Sqrt[epscav]*leff)),"Seconds"^-1]

gampol[gamex_,gamcav_,HX_,HC_]:={HX*gamex +HC*gamcav,HC*gamex+HX*gamcav}

Exk[prc_,type_,ep_][k_]:=UnitConvert[Exres[prc,type,ep]+(HB*Quantity[k,"Micrometers"^(-1)])^2/(2*Mex[prc,type,ep]),"Millielectronvolts"]

Eck[epscav_][Ephot_][k_]:=UnitConvert[Ephot+(HB*Quantity[k,"Micrometers"^(-1)])^2/(2*Mcav[epscav][Ephot]),"Millielectronvolts"]

Ephotpct[prc_,type_,ep_][dE_]:=Exres[prc,type,ep]*(1+dE)

Elpup[Ex_,Ec_,gamex_,gamcav_,V_]:=UnitConvert[{(Ex+Ec+I*HB*(gamex+gamcav))/2-Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)],(Ex+Ec+I*HB*(gamex+gamcav))/2+Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)],2*Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)]},"Millielectronvolts"]

Elp[Ex_,Ec_,gamex_,gamcav_,V_]:=UnitConvert[(Ex+Ec+I*HB*(gamex+gamcav))/2-Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)],"Millielectronvolts"]

Eup[Ex_,Ec_,gamex_,gamcav_,V_]:=UnitConvert[(Ex+Ec+I*HB*(gamex+gamcav))/2+Sqrt[(HB*V)^2+ (1/4)*((Ex-Ec-I*HB*(gamcav-gamex))^2)],"Millielectronvolts"]

Mlpup[Mex_,Mcav_,HX_,HC_]:={(HX/Mex+HC/Mcav)^-1,(HX/Mcav+HC/Mex)^-1}

Exres[prc_,type_,ep_]:=UnitConvert[2*FPegap[prc,type,ep]-Abs[FPev[prc,type,1,ep,1,0]],"Millielectronvolts"]

Vsav[epscav_,n1_,n2_,Nharm_,rmir_][prc_,type_,ep_][Ephot_]:=Sqrt[UnitConvert[(1+Sqrt[rmir])/Sqrt[rmir] ((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*Exres[prc,type,ep]*Leff[epscav,n1,n2,Nharm][Ephot]))*(FPr0[prc,type,1,ep,1,0]^2),"Seconds"^-2]]

Vsavmod[epscav_,n1_,n2_,m_,reff_][prc_,type_,ep_][Ephot_]:=Sqrt[UnitConvert[(1+Sqrt[rmir])/Sqrt[rmir] ((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*Exres[prc,type,ep]*Leff[epscav,n1,n2,m][Ephot]))*(FPr0[prc,type,1,ep,1,0]^2),"Seconds"^-2]]

Vflat[epscav_,n1_,n2_,Nharm_,Ndbr_,rmir_][prc_,type_,ep_][conc2d_,gamex_,epsbg_][Ephot_]:=Sqrt@UnitConvert[(1+Sqrt[rmir])/Sqrt[rmir]*((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2)*conc2d)/(epscav*Exres[prc,type,ep]*Leff[epscav,n1,n2,Nharm][Ephot]))+(((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2)*conc2d)/(CC*Sqrt[epscav]*Exres[prc,type,ep]))*(gamex+gamcav[epscav,rmir,Leff[epscav,n1,n2,Nharm][Ephot]]))+((Exres[prc,type,ep]*Quantity[prc[[1]][[4]],"Nanometers"]*epsbg)/(2*HB*Lc[epscav,Nharm][Ephot]*epscav*Sqrt[rmir]))^2,"Seconds"^-2]

EfromT[temp_]:=UnitConvert[Quantity["BoltzmannConstant"]*Quantity[temp,"Kelvins"],"Millielectronvolts"]

TfromE[E_]:=UnitConvert[E/Quantity["BoltzmannConstant"],"Kelvins"]

U0[prc_,type_,ep_,eps_]:=6*KK*(EE^2)*FPbohrad[prc,type,1,ep,1,0]/eps

Cs[Ueff_,npol_,Mlp_]:=UnitConvert[Sqrt[Ueff*npol/Mlp],"Meters"/"Seconds"]

Tcpol[npol_,cs_,Mp_,s_]:=Module[{Tc0=1/KB*((Pi*HB^2*npol*cs^4*Mp)/(6*s*Zeta[3]))^(1/3)},{((1+Sqrt[32/27 ((Mp*KB*Tc0)/(Pi*HB^2*npol))^3+1])^(1/3)-(Sqrt[32/27 ((Mp*KB*Tc0)/(Pi*HB^2*npol))^3+1]-1)^(1/3)) Tc0/2^(1/3),Tc0}//UnitSimplify]

Tcpol2[npol_,cs_,Mp_,s_]:=Module[{Tc0=1/KB*((2*Pi*HB^2*npol*cs^4*Mp)/(3*s*Zeta[3]))^(1/3)},{((1+Sqrt[32/27 ((Mp*KB*Tc0)/(Pi*HB^2*npol))^3+1])^(1/3)-(Sqrt[32/27 ((Mp*KB*Tc0)/(Pi*HB^2*npol))^3+1]-1)^(1/3)) Tc0/2^(1/3),Tc0}//UnitSimplify]

critdensfromT[tc_,cs_,mp_,s_]:=(Quantity[FindRoot[QuantityMagnitude@Tcpol2[na,cs,mp,s]==tc,{na,10^11}][[1]]//Last,"Meters"^-2])

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
   Leff1 = Function[{type, ep, m, de}, UnitConvert[Lc[epscav, m][Ephotpct[prc, type, ep][de]] + Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]],"Micrometers"]],
   Leff2 = Function[{type, ep, m, de}, UnitConvert[Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 2*Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]],"Micrometers"]],
   MOgamcav = Function[{type, ep, m, de},
			UnitConvert[((1 - Sqrt[r1])/Sqrt[r1])*(CC/(Sqrt[epscav]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))),"Seconds"^-1]
		],
   MNgamcav = Function[{type, ep, m, de},
			UnitConvert[((1 - Sqrt[r1])/Sqrt[r1])*(CC/(Sqrt[epscav]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 2*Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))),"Seconds"^-1]
		],
   OOgamcav = Function[{type, ep, m, de},
	   		UnitConvert[((1 - Sqrt[r2])/Sqrt[r2])*(CC/(2*Sqrt[epscav]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))),"Seconds"^-1]
		],
   ONgamcav = Function[{type, ep, m, de},
			UnitConvert[(1/2*(((1 - Sqrt[r1])/Sqrt[r1]) + ((1 - Sqrt[r2])/Sqrt[r2])))*(CC/(Sqrt[epscav]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + (1/2)*Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))),"Seconds"^-1]
		],
   MOV = Function[{type, ep, m, nx, de},
			Sqrt[nx*UnitConvert[((1 + Sqrt[r1])/Sqrt[r1])*((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*Exres[prc, type, ep]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]])))*(FPr0[prc, type, 1, ep, 1, 0]^2), "Seconds"^-2]]
		],
   MNV = Function[{type, ep, m, nx, de},
			Sqrt[nx*UnitConvert[((1 + Sqrt[r1])/Sqrt[r1])* ((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]],"Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*Exres[prc, type, ep]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + 2*Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]])))*(FPr0[prc, type, 1, ep, 1, 0]^2), "Seconds"^-2]]
		],
   OOV = Function[{type, ep, m, nx, de},
			Sqrt[nx*UnitConvert[((1 + Sqrt[r2])/Sqrt[r2])* ((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]], "Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*Exres[prc, type, ep]*(2*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]))))*(FPr0[prc, type, 1, ep, 1, 0]^2), "Seconds"^-2]]
		],
   ONV = Function[{type, ep, m, nx, de},
			Sqrt[nx*UnitConvert[((1 + Sqrt[(4*r1*r2)/((Sqrt[r1] + Sqrt[r2])^2)])/Sqrt[(4*r1*r2)/((Sqrt[r1] + Sqrt[r2])^2)])* ((4*Pi*KK*(EE^2)*(Quantity[prc[[1]][[3]], "Meters"/"Seconds"]^2))/(Sqrt[epscav*prc[[1]][[6]]]*Exres[prc, type, ep]*(Lc[epscav, m][Ephotpct[prc, type, ep][de]] + (1/2)*Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]])))*(FPr0[prc, type, 1, ep, 1, 0]^2), "Seconds"^-2]]
		],
   U0 = Function[{type, ep}, 6*(Abs@FPev[prc, type, 1, ep, 1, 0])*(FPbohrad[prc, type, 1, ep, 1, 0]^2)],
   HXsq = Function[{V, type, ep, m, nx, de, k},
		HopXsq[
			Exres[prc, type, ep],
			Ephotpct[prc, type, ep][de],
			V[type, ep, m, nx, de]
		]
	],
   HCsq = Function[{V, type, ep, m, nx, de, k},
     HopCsq[
      Exres[prc, type, ep],
      Ephotpct[prc, type, ep][de],
      V[type, ep, m, nx, de]
      ]
     ],
   MLP = Function[{V, type, ep, m, nx, de, k},
		(
			(HopXsq[Exres[prc, type, ep], Ephotpct[prc, type, ep][de], V[type, ep, m, nx, de]]/(4*FPmu[prc, type, ep]) +
			(HopCsq[Exres[prc, type, ep], Ephotpct[prc, type, ep][de], V[type, ep, m, nx, de]]/(Ephotpct[prc, type, ep][de]*epscav/(CC^2)))
		)^(-1))
	]
   },
  Association[
   "inps" -> Association["epscav" -> epscav, "n1" -> n1, "n2" -> n2, "r1"-> r1, "r2"-> r2, "reff"-> Reff, "genname" -> genname],
   "genname" -> ToString@genname,
   "prc" -> prc,
   "NNEtr" -> Function[{t}, NoNegEtr[prc, t]],
   "Emax" -> FPEdim[prc],
   "exres" -> Function[{type, ep}, Exres[prc, type, ep]],
   "lc" -> Function[{type, ep, m, de}, Lc[epscav, m][Ephotpct[prc, type, ep][de]]],
   "ldbr" -> Function[{type, ep}, Ldbr[epscav, n1, n2][Ephotpct[prc, type, ep][0]]],
   "dbrband" -> Function[{type, ep},
	   With[
		{
			ec = Exres[prc, type, ep],
   			bw = Exres[prc, type, ep]*(4/Pi)*ArcSin[(n2 - n1)/(n1 + n2)]
		},
			UnitConvert[{ec, ec + (bw/2), ec - (bw/2)}, "Millielectronvolts"]
      ]
     ],
   "reff" -> Reff,
   "mph" -> Mph,
   "mex" -> Mex,
   "Exk" -> Function[{type, ep, k}, UnitConvert[Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),"Millielectronvolts"]],
   "Eck" -> Function[{type, ep, de, k}, UnitConvert[Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),"Millielectronvolts"]],
   "mono" -> Association[
     "confname" -> StringJoin[genname,", 2 DBR"],
     "old" -> Association[
       "methname" -> "2 DBR, Old",
	   "fullname" -> genname <> ": 2 DBR: Old",
       "leff" -> Function[{type, ep, m, de}, Leff1[type, ep, m, de]],
       "gamcav" -> Function[{type, ep, m, de}, MOgamcav[type, ep, m, de]],
       "V" -> Function[{type, ep, m, nx, de}, MOV[type, ep, m, nx, de]],
       "HX" -> Function[{type, ep, m, nx, de, k}, HXsq[MOV, type, ep, m, nx, de, k]],
       "HC" -> Function[{type, ep, m, nx, de, k}, HCsq[MOV, type, ep, m, nx, de, k]],
       "mlp" -> Function[{type, ep, m, nx, de, k}, (HXsq[MOV, type, ep, m, nx, de, k]/Mex[type, ep] + HCsq[MOV, type, ep, m, nx, de, k]/Mph[type, ep, de])^(-1)],
	   "mlpfcheck" -> Function[{type, ep, m, nx, de, k}, MLP[MOV, type, ep, m, nx, de, k]],
       "mup" -> Function[{type, ep, m, nx, de, k}, (HCsq[MOV, type, ep, m, nx, de, k]/Mex[type, ep] + HXsq[MOV, type, ep, m, nx, de, k]/Mph[type, ep, de])^(-1)],
       "Elp" -> Function[{type, ep, m, nx, gamex, de, k},
         Elp[
          Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          MOgamcav[type, ep, m, de],
          MOV[type, ep, m, nx, de]
          ]
         ],
       "Eup" -> Function[{type, ep, m, nx, gamex, de, k},
         Eup[
          Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          MOgamcav[type, ep, m, de],
          MOV[type, ep, m, nx, de]
          ]
         ],
       "Elpup" -> Function[{type, ep, m, nx, gamex, de, k},
         Eup[
           Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           MOgamcav[type, ep, m, de],
           MOV[type, ep, m, nx, de]
        ] - Elp[
           Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           MOgamcav[type, ep, m, de],
           MOV[type, ep, m, nx, de]
           ]
         ],
       "gamlp" -> Function[{type, ep, m, nx, gamex, de, k},
         (HXsq[MOV, type, ep, m, nx, de, k]*
            gamex) + (HCsq[MOV, type, ep, m, nx, de, k]*
            MOgamcav[type, ep, m, de])
         ],
       "gamup" -> Function[{type, ep, m, nx, gamex, de, k},
         (HCsq[MOV, type, ep, m, nx, de, k]*gamex) + (HXsq[MOV, type, ep, m, nx, de, k]*MOgamcav[type, ep, m, de])
         ],
       "Ueff" -> Function[{type, ep, m, nx, de}, U0[type, ep]*(HXsq[MOV, type, ep, m, nx, de, 0]^2)],
       "cs" -> Function[{type, ep, m, nx, npol, de},
         Sqrt[
          U0[type, ep]*(HXsq[MOV, type, ep, m, nx, de, 0]^2)*
           npol/(MLP[MOV, type, ep, m, nx, de, 0])]],
       "tc0" -> Function[{type, ep, m, nx, npol, de}
         tc0[
          MLP[MOV, type, ep, m, nx, de, 0],
          npol,
          U0[type, ep]*(HXsq[MOV, type, ep, m, nx, de, 0]^2)
          ]
         ],
       "ncrit" -> Function[{type, ep, m, nx, t, de},
         ncrit[
          MLP[MOV, type, ep, m, nx, de, 0],
          U0[type, ep]*(HXsq[MOV, type, ep, m, nx, de, 0]^2),
          t
          ]
         ],
       "tcbkt" -> Function[{type, ep, m, nx, npol, de},
         Tcpostsolve[
          MLP[MOV, type, ep, m, nx, de, 0],
          npol,
          U0[type, ep]*(HXsq[MOV, type, ep, m, nx, de, 0]^2)
          ]
         ],
       "tcmax" -> Function[{t, ep, m, nx, conc},
         Module[{max},
          max = FindMaximum[
            QuantityMagnitude[
             (Tcpol2[
                conc,
                Sqrt[U0[t, ep]*(HXsq[MOV, t, ep, m, nx, dE, 0]^2)*
                  conc/(MLP[MOV, t, ep, m, nx, dE, 0])],
				  MLP[MOV, t, ep, m, nx, dE, 0],
                1
                ])[[1]]],
            {dE, 0.001},
			AccuracyGoal->4
            ];
          {dE /. Last@max, max[[1]]}
          ]
         ],
		 "SCstart"->Function[{t,m,nx,gamex},
			Module[{j=1},
				While[
					QuantityMagnitude[MOV[t,j,m,nx,0,0]] < QuantityMagnitude[Abs[(MOgamcav[t,j,m,0]-gamex)/2]],
					j++;
					If[j==29,Break[]]
				];
				j
			 ]
       ]
	],
     "new" -> Association[
       "methname" -> "2 DBR, New",
	   "fullname" -> genname <> ": 2 DBR: New",
       "leff" -> Function[{type, ep, m, de}, Leff2[type, ep, m, de]],
       "gamcav" -> Function[{type, ep, m, de}, MNgamcav[type, ep, m, de]],
       "V" -> Function[{type, ep, m, nx, de}, MNV[type, ep, m, nx, de]],
       "HX" -> Function[{type, ep, m, nx, de, k}, HXsq[MNV, type, ep, m, nx, de, k]],
       "HC" -> Function[{type, ep, m, nx, de, k}, HCsq[MNV, type, ep, m, nx, de, k]],
       "mlp" -> Function[{type, ep, m, nx, de, k}, (HXsq[MNV, type, ep, m, nx, de, k]/Mex[type, ep] + HCsq[MNV, type, ep, m, nx, de, k]/Mph[type, ep, de])^(-1)],
	   "mlpfcheck" -> Function[{type, ep, m, nx, de, k}, MLP[MNV, type, ep, m, nx, de, k]],
       "mup" -> Function[{type, ep, m, nx, de, k}, (HCsq[MNV, type, ep, m, nx, de, k]/Mex[type, ep] + HXsq[MNV, type, ep, m, nx, de, k]/Mph[type, ep, de])^(-1)],
       "Elp" -> Function[{type, ep, m, nx, gamex, de, k},
         Elp[
          Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          MNgamcav[type, ep, m, de],
          MNV[type, ep, m, nx, de]
          ]
         ],
       "Eup" -> Function[{type, ep, m, nx, gamex, de, k},
         Eup[
          Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          MNgamcav[type, ep, m, de],
          MNV[type, ep, m, nx, de]
          ]
         ],
       "Elpup" -> Function[{type, ep, m, nx, gamex, de, k},
         Eup[
           Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           MNgamcav[type, ep, m, de],
           MNV[type, ep, m, nx, de]
	   ] - Elp[
           Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           MNgamcav[type, ep, m, de],
           MNV[type, ep, m, nx, de]
           ]
         ],
       "gamlp" -> Function[{type, ep, m, nx, gamex, de, k},
         (HXsq[MNV, type, ep, m, nx, de, k]*
            gamex) + (HCsq[MNV, type, ep, m, nx, de, k]*
            MNgamcav[type, ep, m, de])
         ],
       "gamup" -> Function[{type, ep, m, nx, gamex, de, k},
         (HCsq[MNV, type, ep, m, nx, de, k]*
            gamex) + (HXsq[MNV, type, ep, m, nx, de, k]*
            MNgamcav[type, ep, m, de])
         ],
       "Ueff" -> Function[{type, ep, m, nx, de}, U0[type, ep]*(HXsq[MNV, type, ep, m, nx, de, 0]^2)],
       "cs" -> Function[{type, ep, m, nx, npol, de}, Sqrt[U0[type, ep]*(HXsq[MNV, type, ep, m, nx, de, 0]^2)*npol/(MLP[MNV, type, ep, m, nx, de, 0])]],
       "tc0" -> Function[{type, ep, m, nx, npol, de},
         tc0[
          MLP[MNV, type, ep, m, nx, de, 0],
          npol,
          U0[type, ep]*(HXsq[MNV, type, ep, m, nx, de, 0]^2)
          ]
         ],
       "ncrit" -> Function[{type, ep, m, nx, t, de},
         ncrit[
          MLP[MNV, type, ep, m, nx, de, 0],
          U0[type, ep]*(HXsq[MNV, type, ep, m, nx, de, 0]^2),
          t
          ]
         ],
       "tcbkt" -> Function[{type, ep, m, nx, npol, de},
         Tcpostsolve[
          MLP[MNV, type, ep, m, nx, de, 0],
          npol,
          U0[type, ep]*(HXsq[MNV, type, ep, m, nx, de, 0]^2)
          ]
         ],
       "tcmax" -> Function[{t, ep, m, nx, conc},
         Module[{max},
          max = FindMaximum[
            QuantityMagnitude[
             (Tcpol2[
                conc,
                Sqrt[U0[t, ep]*(HXsq[MNV, t, ep, m, nx, dE, 0]^2)*conc/(MLP[MNV, t, ep, m, nx, dE, 0])],
                MLP[MNV, t, ep, m, nx, dE, 0],
                1
                ])[[1]]],
            {dE, 0.01},
			AccuracyGoal->4
            ];
          {dE /. Last@max, max[[1]]}
          ]
         ],
		 "SCstart"->Function[{t,m,nx,gamex},
			Module[{j=1},
				While[
					QuantityMagnitude[MNV[t,j,m,nx,0,0]] < QuantityMagnitude[Abs[(MNgamcav[t,j,m,0]-gamex)/2]],
					j++;
					If[j==29,Break[]]
					];
				j
				]
       ]
     ]
	 ],
   "open" -> Association[
     "confname" -> StringJoin[genname, ", 1 DBR"],
     "old" -> Association[
       "methname" -> "1 DBR, Old",
	   "fullname" -> genname <> ": 1 DBR: Old",
       "leff" -> Function[{type, ep, m, de}, Leff1[type, ep, m, de]],
       "gamcav" -> Function[{type, ep, m, de}, OOgamcav[type, ep, m, de]],
       "V" -> Function[{type, ep, m, nx, de}, OOV[type, ep, m, nx, de]],
       "HX" -> Function[{type, ep, m, nx, de, k}, HXsq[OOV, type, ep, m, nx, de, k]],
       "HC" -> Function[{type, ep, m, nx, de, k}, HCsq[OOV, type, ep, m, nx, de, k]],
       "mlp" -> Function[{type, ep, m, nx, de, k}, (HXsq[OOV, type, ep, m, nx, de, k]/Mex[type, ep] + HCsq[OOV, type, ep, m, nx, de, k]/Mph[type, ep, de])^(-1)],
	   "mlpfcheck" -> Function[{type, ep, m, nx, de, k}, MLP[OOV, type, ep, m, nx, de, k]],
       "mup" -> Function[{type, ep, m, nx, de, k}, (HCsq[OOV, type, ep, m, nx, de, k]/Mex[type, ep] + HXsq[OOV, type, ep, m, nx, de, k]/Mph[type, ep, de])^(-1)],
       "Elp" -> Function[{type, ep, m, nx, gamex, de, k},
         Elp[
          Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          OOgamcav[type, ep, m, de],
          OOV[type, ep, m, nx, de]
          ]
         ],
       "Eup" -> Function[{type, ep, m, nx, gamex, de, k},
         Eup[
          Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          OOgamcav[type, ep, m, de],
          OOV[type, ep, m, nx, de]
          ]
         ],
       "Elpup" -> Function[{type, ep, m, nx, gamex, de, k},
         Eup[
           Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           OOgamcav[type, ep, m, de],
           OOV[type, ep, m, nx, de]
           ] - Elp[
           Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           OOgamcav[type, ep, m, de],
           OOV[type, ep, m, nx, de]
           ]
         ],
       "gamlp" -> Function[{type, ep, m, nx, gamex, de, k},
         (HXsq[OOV, type, ep, m, nx, de, k]*
            gamex) + (HCsq[OOV, type, ep, m, nx, de, k]*
            OOgamcav[type, ep, m, de])
         ],
       "gamup" -> Function[{type, ep, m, nx, gamex, de, k},
         (HCsq[OOV, type, ep, m, nx, de, k]*
            gamex) + (HXsq[OOV, type, ep, m, nx, de, k]*
            OOgamcav[type, ep, m, de])
         ],
       "Ueff" -> Function[{type, ep, m, nx, de}, U0[type, ep]*(HXsq[OOV, type, ep, m, nx, de, 0]^2)],
       "cs" -> Function[{type, ep, m, nx, npol, de}, Sqrt[U0[type, ep]*(HXsq[OOV, type, ep, m, nx, de, 0]^2)*npol/(MLP[OOV, type, ep, m, nx, de, 0])]],
       "tc0" -> Function[{type, ep, m, nx, npol, de},
         tc0[
          MLP[OOV, type, ep, m, nx, de, 0],
          npol,
          U0[type, ep]*(HXsq[OOV, type, ep, m, nx, de, 0]^2)
          ]
         ],
       "ncrit" -> Function[{type, ep, m, nx, t, de},
         ncrit[
          MLP[OOV, type, ep, m, nx, de, 0],
          U0[type, ep]*(HXsq[OOV, type, ep, m, nx, de, 0]^2),
          t
          ]
         ],
       "tcbkt" -> Function[{type, ep, m, nx, npol, de},
         Tcpostsolve[
          MLP[OOV, type, ep, m, nx, de, 0],
          npol,
          U0[type, ep]*(HXsq[OOV, type, ep, m, nx, de, 0]^2)
          ]
         ],
       "tcmax" -> Function[{t, ep, m, nx, conc},
         Module[{max},
          max = FindMaximum[
            QuantityMagnitude[
             (Tcpol2[
                conc,
                Sqrt[U0[t, ep]*(HXsq[OOV, t, ep, m, nx, dE, 0]^2)*conc/(MLP[OOV, t, ep, m, nx, dE, 0])],
                MLP[OOV, t, ep, m, nx, dE, 0],
                1
                ])[[1]]],
            {dE, 0.01},
			AccuracyGoal->4
            ];
          {dE /. Last@max, max[[1]]}
          ]
         ],
		 "SCstart"->Function[{t,m,nx,gamex},
			Module[{j=1},
				While[
					QuantityMagnitude[OOV[t,j,m,nx,0,0]] < QuantityMagnitude[Abs[(OOgamcav[t,j,m,0]-gamex)/2]],
					j++;
					If[j==29,Break[]]
				];
				j
			]
		]
       ],
     "new" -> Association[
       "methname" -> "1 DBR, New",
	   "fullname" -> genname <> ": 1 DBR: New",
       "leff" -> Function[{type, ep, m, de}, Leff1[type, ep, m, de]],
       "gamcav" -> Function[{type, ep, m, de}, ONgamcav[type, ep, m, de]],
       "V" -> Function[{type, ep, m, nx, de}, ONV[type, ep, m, nx, de]],
       "HX" -> Function[{type, ep, m, nx, de, k}, HXsq[ONV, type, ep, m, nx, de, k]],
       "HC" -> Function[{type, ep, m, nx, de, k}, HCsq[ONV, type, ep, m, nx, de, k]],
       "mlp" -> Function[{type, ep, m, nx, de, k}, (HXsq[ONV, type, ep, m, nx, de, k]/Mex[type, ep] + HCsq[ONV, type, ep, m, nx, de, k]/Mph[type, ep, de])^(-1)],
	   "mlpfcheck" -> Function[{type, ep, m, nx, de, k}, MLP[ONV, type, ep, m, nx, de, k]],
       "mup" -> Function[{type, ep, m, nx, de, k}, (HCsq[ONV, type, ep, m, nx, de, k]/Mex[type, ep] + HXsq[ONV, type, ep, m, nx, de, k]/Mph[type, ep, de])^(-1)],
       "Elp" -> Function[{type, ep, m, nx, gamex, de, k},
         Elp[
          Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          ONgamcav[type, ep, m, de],
          ONV[type, ep, m, nx, de]
          ]
         ],
       "Eup" -> Function[{type, ep, m, nx, gamex, de, k},
         Eup[
          Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
          Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
          gamex,
          ONgamcav[type, ep, m, de],
          ONV[type, ep, m, nx, de]
          ]
         ],
       "Elpup" -> Function[{type, ep, m, nx, gamex, de, k},
         Eup[
           Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           ONgamcav[type, ep, m, nx, de],
           ONV[type, ep, m, nx, de]
           ] - Elp[
           Exres[prc, type, ep] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mex[type, ep]),
           Ephotpct[prc, type, ep][de] + ((HB^2)*Quantity[k, "Micrometers"^-1]^2)/(2*Mph[type, ep, de]),
           gamex,
           ONgamcav[type, ep, m, de],
           ONV[type, ep, m, nx, de]
           ]
         ],
       "gamlp" -> Function[{type, ep, m, nx, gamex, de, k},
         (HXsq[ONV, type, ep, m, nx, de, k]*
            gamex) + (HCsq[ONV, type, ep, m, nx, de, k]*
            ONgamcav[type, ep, m, de])
         ],
       "gamup" -> Function[{type, ep, m, nx, gamex, de, k},
         (HCsq[ONV, type, ep, m, nx, de, k]*
            gamex) + (HXsq[ONV, type, ep, m, nx, de, k]*
            ONgamcav[type, ep, m, de])
         ],
       "Ueff" -> Function[{type, ep, m, nx, de}, U0[type, ep]*(HXsq[ONV, type, ep, m, nx, de, 0]^2)],
       "cs" -> Function[{type, ep, m, nx, npol, de}, Sqrt[U0[type, ep]*(HXsq[ONV, type, ep, m, nx, de, 0]^2)*npol/(MLP[ONV, type, ep, m, nx, de, 0])]],
       "tc0" -> Function[{type, ep, m, nx, npol, de},
         tc0[
          MLP[ONV, type, ep, m, nx, de, 0],
          npol,
          U0[type, ep]*(HXsq[ONV, type, ep, m, nx, de, 0]^2)
          ]
         ],
       "ncrit" -> Function[{type, ep, m, nx, t, de},
         ncrit[
          MLP[ONV, type, ep, m, nx, de, 0],
          U0[type, ep]*(HXsq[ONV, type, ep, m, nx, de, 0]^2),
          t
          ]
         ],
       "tcbkt" -> Function[{type, ep, m, nx, npol, de},
         Tcpostsolve[
          MLP[ONV, type, ep, m, nx, de, 0],
          npol,
          U0[type, ep]*(HXsq[ONV, type, ep, m, nx, de, 0]^2)
          ]
         ],
       "tcmax" -> Function[{t, ep, m, nx, conc},
         Module[{max},
          max = FindMaximum[
            QuantityMagnitude[
             (Tcpol2[
                conc,
                Sqrt[U0[t, ep]*(HXsq[ONV, t, ep, m, nx, dE, 0]^2)*conc/(MLP[ONV, t, ep, m, nx, dE, 0])],
                MLP[ONV, t, ep, m, nx, dE, 0],
                1
                ])[[1]]],
            {dE, 0.01},
			AccuracyGoal->4
            ];
          {dE /. Last@max, max[[1]]}
          ]
         ],
		 "SCstart"->Function[{t,m,nx,gamex},
			Module[{j=1},
				While[
					QuantityMagnitude[ONV[t,j,m,nx,0,0]] < QuantityMagnitude[Abs[(ONgamcav[t,j,m,0]-gamex)/2]],
					j++;
					If[j==29,Break[]]
				];
				j
				]
			]
		]
   ]
   ]
]
