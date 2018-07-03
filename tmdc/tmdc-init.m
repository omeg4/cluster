#!/usr/bin/env wolframscript
(* ::Package:: *)

(* analytical solutions: *)
\[Beta][\[Mu]_,eps_]:= (2 * \[Mu]/eps);
\[Gamma][D_,\[Mu]_,eps_] := Sqrt[(\[Mu])/(eps*(D^3))]
gamkeld[mu_,kappa_,rho_,D_]:=-(Pi/(4*kappa*(rho^2)*D))*(StruveH[-1,D/rho]-BesselY[-1,D/rho])
a[mu_,kappa_,rho_,D_]:=Sqrt[1/(2*Sqrt[2*mu*gamkeld[mu,kappa,rho,D]])]
v0keld[mu_,kappa_,rho_,D_]:=(Pi/(2*kappa*rho))*(StruveH[0,D/rho]-BesselY[0,D/rho])
PsiDir[r_,\[Phi]_,n_,l_,\[Mu]_,eps_,d_] := ((\[Beta][\[Mu],eps]/(n-0.5))^(Abs[l]+1))*Sqrt[Factorial[n - Abs[l] - 1]/(Factorial[n+Abs[l]-1]*(2*n-1))]* (r^Abs[l]) *Exp[-( \[Beta][\[Mu],eps]/(n-0.5))*r/2]*LaguerreL[n-Abs[l]-1,2*Abs[l],\[Beta][\[Mu],eps]*r/(n-0.5)]*Exp[-I*l*\[Phi]]*Sqrt[1/(2*Pi)]
PsiInd[r_,\[Phi]_,n_,l_,\[Mu]_,eps_,d_]:= Sqrt[2*Factorial[(n-Abs[l])/2]/(Factorial[(n+Abs[l])/2])]*Sqrt[\[Gamma][d,\[Mu],eps]/2]*Exp[-\[Gamma][d,\[Mu],eps]*(r^2)/2]*((\[Gamma][d,\[Mu],eps]*(r^2))^(Abs[l]/2))*LaguerreL[(n-Abs[l])/2,Abs[l],\[Gamma][d,\[Mu],eps]*r^2]*Exp[I*l*\[Phi]]*((-1)^((n-Abs[l])/2))*(1/Sqrt[Pi])
PsiKeld[r_,phi_,n_,l_,mu_,kappa_,rho_,d_]:=(r^Abs[l])*Exp[-(r^2)/(4*a[mu,kappa,rho,d]^2)]*LaguerreL[(n-Abs[l])/2,Abs[l],(r^2)/(2*a[mu,kappa,rho,d]^2)]*Exp[I*l*phi]
PsiKeldNorm[r_,phi_,n_,l_,mu_,kappa_,rho_,d_]:=Module[
{norm},
norm=Integrate[rr*Conjugate[PsiKeld[rr,phii,n,l,mu,kappa,rho,d]]*PsiKeld[rr,phii,n,l,mu,kappa,rho,d],{rr,0,Infinity},{phii,0,2*Pi}];
(1/Sqrt[norm])*PsiKeld[r,phi,n,l,mu,kappa,rho,d]
]
enDir[n_,\[Mu]_,eps_,d_]:=-(\[Mu]/(2*(eps^2)*(n-1/2)^2)) 
enInd[n_,\[Mu]_,eps_,D_]:= SetPrecision[Sqrt[1/(\[Mu]*eps*(D^3))]*(n + 1) - (1/(eps*D)),30]
enKeld[n_,mu_,kappa_,rho_,D_]:=-v0keld[mu,kappa,rho,D]+ Sqrt[2*gamkeld[mu,kappa,rho,D]/mu]*(n+1)


(* ::Subsubsection:: *)
(*Conversion Factors:*)


B2nm=0.052917721092; (*nm/bohr*)
H2J = 4.35974417*(10^-18); (*J/Subscript[E, h]*)
H2eV = 27.21138602; (*eV/Subscript[E, h]*)
T2S = 2.41884326505*10^-17; (*\[HBar]/Subscript[E, h] ;; seconds/[atomic time]*)
lBN=(1/B2nm)*.333; (* Thickness of hBN in units of Bohr *)
(* hBN is 3.33 ANGS *)


(* Generic values for absorption: na, \[CapitalGamma], etc. *)
na = (5*10^15)*(10^-18)*(B2nm^2) ;
damp = (10^13)*(T2S);
damp2=(10^12)*T2S;
damp3 = (10^11)*(T2S);
l = (3.1 * 10^-10)*(10^9)/(B2nm) ;


compsysKeld[maxm_,mew_,kappa_,rhorho_]:=Module[
{
(*pltarray,evTab,efTab,efTabRenorm,ef,ev,evNew,
Kevdir,Kefdir,Cevdir,Cefdir,Cevind,Cefind,bigarray,evout,efout,*)
maxcell=10^-4,
maxiter=4*10^5,
(*e1=eee1,
e2=eee2,*)
\[Kappa]=kappa,
\[Rho]=rhorho,
e=1,
d=10,
\[Mu]=mew,
shift=10,
mmax=maxm,
radialeqs,
solnmat
(*radialEqKeld,radialEqDir,radialEqInd,radial\[Xi]Keld,radial\[Xi]Dir,radial\[Xi]Ind*)
},
radialEqKeld=-(1/(\[Mu]*2))f''[r]-(1/(2*\[Mu]*r))* f'[r]-(((((Pi*(e^2))/(2*\[Kappa]*\[Rho]))*(StruveH[0,r/\[Rho]]-BesselY[0,r/\[Rho]])))-(m^2/(2*\[Mu]* r^2)))* f[r];
radial\[Xi]Keld[m_]=Simplify[radialEqKeld/.f->(\[Psi][ArcTan[#]]&)/.r->(Tan[\[Xi]]),Pi/2>\[Xi]>0];
solnmat={};
evTab={};
efTab={};
bigarray={};
evout={};
efout={};
Do[
ev={};
{ev,ef}=NDEigensystem[{radial\[Xi]Keld[mind]+shift \[Psi][\[Xi]],DirichletCondition[\[Psi][\[Xi]]==0,\[Xi]==Pi/2]},\[Psi][\[Xi]],{\[Xi],0,Pi/2},mmax-mind+1,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->maxcell}}},"Eigensystem"->{"Arnoldi",MaxIterations->maxiter}}];
evTab=Append[evTab,ev-shift];
efTab=Append[efTab,ef],{mind,0,mmax}
];
{Table[evTab[[i-j+1]][[j]],{i,Dimensions[evTab][[1]]},{j,i,1,-1}],Table[efTab[[i-j+1]][[j]],{i,Dimensions[efTab][[1]]},{j,i,1,-1}]}
]


compsysDir[maxm_,mew_,eeps_]:=Module[
{
(*pltarray,evTab,efTab,efTabRenorm,ef,ev,evNew,
Kevdir,Kefdir,Cevdir,Cefdir,Cevind,Cefind,bigarray,evout,efout,*)
maxcell=10^-4,
maxiter=4*10^5,
eps=eeps,
\[Mu]=mew,
shift=10,
mmax=maxm,
radialeqs,
solnmat
(*radialEqKeld,radialEqDir,radialEqInd,radial\[Xi]Keld,radial\[Xi]Dir,radial\[Xi]Ind*)
},
radialEqDir=-(1/(\[Mu]*2))f''[r]-(1/(2\[Mu]*r))* f'[r]-((1/(eps*r))-(m^2/(2*\[Mu]* r^2)))* f[r];
radial\[Xi]Dir[m_]=Simplify[radialEqDir/.f->(\[Psi][ArcTan[#]]&)/.r->(Tan[\[Xi]]),Pi/2>\[Xi]>0];
solnmat={};
evTab={};
efTab={};
bigarray={};
evout={};
efout={};
Do[
ev={};{ev,ef}=NDEigensystem[{radial\[Xi]Dir[mind]+shift \[Psi][\[Xi]],DirichletCondition[\[Psi][\[Xi]]==0,\[Xi]==Pi/2]},\[Psi][\[Xi]],{\[Xi],0,Pi/2},mmax-mind+1,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->maxcell}}},"Eigensystem"->{"Arnoldi",MaxIterations->maxiter}}];
evTab=Append[evTab,ev-shift];
efTab=Append[efTab,ef],{mind,0,mmax}
];
{Table[evTab[[i-j+1]][[j]],{i,Dimensions[evTab][[1]]},{j,i,1,-1}],Table[efTab[[i-j+1]][[j]],{i,Dimensions[efTab][[1]]},{j,i,1,-1}]}
]


compsysInd[maxm_,mew_,eeps_,dee_]:=Module[
{
(*pltarray,evTab,efTab,efTabRenorm,ef,ev,evNew,
Kevdir,Kefdir,Cevdir,Cefdir,Cevind,Cefind,bigarray,evout,efout,*)
maxcell=10^-4,
maxiter=4*10^5,
eps=eeps,
d=dee,
\[Mu]=mew,
shift=10,
mmax=maxm,
radialeqs,
solnmat
(*radialEqKeld,radialEqDir,radialEqInd,radial\[Xi]Keld,radial\[Xi]Dir,radial\[Xi]Ind*)
},
radialEqInd=-(1/(\[Mu]*2))f''[r]-(1/(2\[Mu]*r))* f'[r]-((1/Sqrt[r^2 + d^2])-(m^2/(2*\[Mu]* r^2)))* f[r];
radial\[Xi]Ind[m_]=Simplify[radialEqInd/.f->(\[Psi][ArcTan[#]]&)/.r->(Tan[\[Xi]]),Pi/2>\[Xi]>0];
solnmat={};
evTab={};
efTab={};
bigarray={};
Do[
ev={};{ev,ef}=NDEigensystem[{radial\[Xi]Ind[mind]+shift \[Psi][\[Xi]],DirichletCondition[\[Psi][\[Xi]]==0,\[Xi]==Pi/2]},\[Psi][\[Xi]],{\[Xi],0,Pi/2},mmax-mind+1,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->maxcell}}},"Eigensystem"->{"Arnoldi",MaxIterations->maxiter}}];
evTab=Append[evTab,ev-shift];
efTab=Append[efTab,ef],
{mind,0,mmax}
];
{Table[evTab[[i-j+1]][[j]],{i,Dimensions[evTab][[1]]},{j,i,1,-1}],Table[efTab[[i-j+1]][[j]],{i,Dimensions[efTab][[1]]},{j,i,1,-1}]}
]


compsysKeldInd[maxm_,mew_,kappa_,rhorho_,dee_]:=Module[
{
\[Kappa]=kappa,
\[Rho]=rhorho,
d=dee,
\[Mu]=mew,
mmax=maxm,
maxcell=10^-4,
maxiter=4*10^5,
shift=10,
radialeqs,
solnmat
},
radialEqKInd=-(1/(\[Mu]*2))f''[r]-(1/(2\[Mu]*r))* f'[r]-(((Pi/(2*\[Kappa]*\[Rho])*(StruveH[0,Sqrt[r^2+d^2]/\[Rho]]-BesselY[0,Sqrt[r^2+d^2]/\[Rho]])))-(m^2/(2*\[Mu]* r^2)))* f[r];
radial\[Xi]KInd[m_]=Simplify[radialEqKInd/.f->(\[Psi][ArcTan[#]]&)/.r->(Tan[\[Xi]]),Pi/2>\[Xi]>0];
solnmat={};
evTab={};
efTab={};
bigarray={};
Do[
ev={};
{ev,ef}=NDEigensystem[{radial\[Xi]KInd[mind]+shift \[Psi][\[Xi]],DirichletCondition[\[Psi][\[Xi]]==0,\[Xi]==Pi/2]},\[Psi][\[Xi]],{\[Xi],0,Pi/2},mmax-mind+1,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->maxcell}}},"Eigensystem"->{"Arnoldi",MaxIterations->maxiter}}];
evTab=Append[evTab,ev-shift];
efTab=Append[efTab,ef],
{mind,0,mmax}
];
{Table[evTab[[i-j+1]][[j]],{i,Dimensions[evTab][[1]]},{j,i,1,-1}],Table[efTab[[i-j+1]][[j]],{i,Dimensions[efTab][[1]]},{j,i,1,-1}]}
]


compsysHO[maxm_,mew_,eeps_,dee_]:=Module[
{
(*pltarray,evTab,efTab,efTabRenorm,ef,ev,evNew,
Kevdir,Kefdir,Cevdir,Cefdir,Cevind,Cefind,bigarray,evout,efout,*)
maxcell=10^-4,
maxiter=4*10^5,
eps=eeps,
e=1,
d=dee,
\[Mu]=mew,
shift=10,
mmax=maxm,
radialeqs,
solnmat
(*radialEqKeld,radialEqDir,radialEqInd,radial\[Xi]Keld,radial\[Xi]Dir,radial\[Xi]Ind*)
},
radialEqHO=-(1/(\[Mu]*2))f''[r]-(1/(2\[Mu]*r))* f'[r]-((1/d)-(r^2/(2*d^3))-(m^2/(2*\[Mu]* r^2)))* f[r];
radial\[Xi]HO[m_]=Simplify[radialEqHO/.f->(\[Psi][ArcTan[#]]&)/.r->(Tan[\[Xi]]),Pi/2>\[Xi]>0];
solnmat={};
evTab={};
efTab={};
bigarray={};
Do[
ev={};{ev,ef}=NDEigensystem[{radial\[Xi]HO[mind]+shift \[Psi][\[Xi]],DirichletCondition[\[Psi][\[Xi]]==0,\[Xi]==Pi/2]},\[Psi][\[Xi]],{\[Xi],0,Pi/2},mmax-mind+1,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->maxcell}}},"Eigensystem"->{"Arnoldi",MaxIterations->maxiter}}];
evTab=Append[evTab,ev-shift];
efTab=Append[efTab,ef],
{mind,0,mmax}
];
{ReorgI2[evTab],ReorgI2[efTab]}
]


compsysKeldHO[maxm_,mew_,kappa_,rhorho_,dee_]:=Module[
{
(*pltarray,evTab,efTab,efTabRenorm,ef,ev,evNew,
Kevdir,Kefdir,Cevdir,Cefdir,Cevind,Cefind,bigarray,evout,efout,*)
maxcell=10^-4,
maxiter=4*10^5,
\[Kappa]=kappa,
\[Rho]=rhorho,
d=dee,
\[Mu]=mew,
shift=10,
mmax=maxm,
radialeqs,
solnmat,
zerosmat
(*radialEqKeld,radialEqDir,radialEqInd,radial\[Xi]Keld,radial\[Xi]Dir,radial\[Xi]Ind*)
},
v0 =(Pi/(2*\[Kappa]*\[Rho]))*(StruveH[0,d/\[Rho]]-BesselY[0,d/\[Rho]]);
gamma=(Pi/(4*\[Kappa]*(\[Rho]^2)*d))*(StruveH[-1,d/\[Rho]]-BesselY[-1,d/\[Rho]]);
radialEqKHO=-(1/(\[Mu]*2))f''[r]-(1/(2\[Mu]*r))* f'[r]-(v0+gamma*r^2-(m^2/(2*\[Mu]* r^2)))* f[r];
radial\[Xi]KHO[m_]=Simplify[radialEqKHO/.f->(\[Psi][ArcTan[#]]&)/.r->(Tan[\[Xi]]),Pi/2>\[Xi]>0];
solnmat={};
evTab={};
efTab={};
bigarray={};
zerosmat=
Do[
ev={};{ev,ef}=NDEigensystem[{radial\[Xi]KHO[mind]+shift \[Psi][\[Xi]],DirichletCondition[\[Psi][\[Xi]]==0,\[Xi]==Pi/2]},\[Psi][\[Xi]],{\[Xi],0,Pi/2},mmax-mind+1,Method->{"SpatialDiscretization"->{"FiniteElement",{"MeshOptions"->{"MaxCellMeasure"->maxcell}}},"Eigensystem"->{"Arnoldi",MaxIterations->maxiter}}];
evTab=Append[evTab,ev-shift];
efTab=Append[efTab,ef],
{mind,0,mmax}
];
{evTab,efTab}
]


exactsuite[inmat_,kappa_,drange_,mmax_,conc_,damping_,runprefix_,filext_]:=
Module[
{
leadprefix,
alldata={}
},
leadprefix=StringForm["`1`-`2`-`3`_`4`",DateList[][[1]],DateList[][[2]],DateList[][[3]],runprefix];
alldata=Table[
Module[
{
mu = inmat[[i]][[1]],
rho = 2*Pi*(inmat[[i]][[2]]/(B2nm*10*kappa)), (* matrix element is \[Chi] in Angstroms, this converts to \[Rho] = (2\[Pi]\[Chi])/\[Kappa] in units of bohr *)
l = inmat[[i]][[3]]/(B2nm*10),
matprefix = inmat[[i]][[4]],
Di = drange[[1]],
Df = drange[[2]],
Ds = drange[[3]],
na = conc, (* (bohr^-2)*)
damp = damping,
fileprefix,
dlist,
evXeV,
efXnorm,
exactsys
      },
fileprefix=StringForm["`1`_`2`",leadprefix,matprefix];
dlist=Table[d,{d,l,l+(Df*lBN),(Ds*lBN)}];
exactsys=Table[
Module[
{tempev,tempef,tempefnorm},
{tempev,tempef}=compsysKeldInd[mmax,mu,kappa,rho,d];
tempefnorm=norm[tempef,10^4,20][[1]];
{
tempev,
tempefnorm
}
],
{d,l,l+(Df*lBN),(Ds*lBN)}
];
{
Prepend[Table[{
dlist[[i]],
exactsys[[i]][[1]]
},
{i,Df}],
ToString[matprefix]
],
Prepend[Table[{
dlist[[i]],
exactsys[[i]][[2]]
},
{i,Df}],
ToString[matprefix]
]
}
],
{i,Dimensions[inmat][[1]]}
];
{
Table[alldata[[i]][[1]],{i,Dimensions[inmat][[1]]}],
Table[alldata[[i]][[2]],{i,Dimensions[inmat][[1]]}]
}
]


exactCOULsuite[inmat_,kappa_,drange_,mmax_,conc_,damping_,runprefix_,filext_]:=
Module[
{
leadprefix,
alldata={}
},
leadprefix=StringForm["`1`-`2`-`3`_`4`",DateList[][[1]],DateList[][[2]],DateList[][[3]],runprefix];
alldata=Table[
Module[
{
mu = inmat[[i]][[1]],
rho = 2*Pi*(inmat[[i]][[2]]/(B2nm*10*kappa)), (* matrix element is \[Chi] in Angstroms, this converts to \[Rho] = (2\[Pi]\[Chi])/\[Kappa] in units of bohr *)
l = inmat[[i]][[3]]/(B2nm*10),
matprefix = inmat[[i]][[4]],
Di = drange[[1]],
Df = drange[[2]],
Ds = drange[[3]],
na = conc, (* (bohr^-2)*)
damp = damping,
fileprefix,
dlist,
evXeV,
efXnorm,
exactsys
      },
fileprefix=StringForm["`1`_`2`",leadprefix,matprefix];
dlist=Table[d,{d,l,l+(Df*lBN),(Ds*lBN)}];
exactsys=Table[
Module[
{tempev,tempef,tempefnorm},
{tempev,tempef}=compsysInd[mmax,mu,kappa,d];
tempefnorm=norm[tempef,10^4,20][[1]];
{
tempev,
tempefnorm
}
],
{d,l,l+(Df*lBN),(Ds*lBN)}
];
{
Prepend[Table[{
dlist[[i]],
exactsys[[i]][[1]]
},
{i,Df}],
ToString[matprefix]
],
Prepend[Table[{
dlist[[i]],
exactsys[[i]][[2]]
},
{i,Df}],
ToString[matprefix]
]
}
],
{i,Dimensions[inmat][[1]]}
];
{
Table[alldata[[i]][[1]],{i,Dimensions[inmat][[1]]}],
Table[alldata[[i]][[2]],{i,Dimensions[inmat][[1]]}]
}
]


(* ::Input::Initialization:: *)
EXmatrange={
{
0.28,
6.6,
6.18,
"mos2-hi"
},
{
0.16,
7.112,
6.18,
"mos2-lo"
},
{
0.31,
8.23,
6.527,
"mose2-hi"
},
{
0.27,
8.461,
6.527,
"mose2-lo"
},
{
0.23,
6.03,
6.219,
"ws2-hi"
},
{
0.15,
6.393,
6.219,
"ws2-lo"
},
{
0.27,
7.18,
6.575,
"wse2-hi"
},
{
0.15,
7.571,
6.575,
"wse2-lo"
}
};


(* ::Input::Initialization:: *)
extracthiloEB[matrix_,materialind_,dind_]:=matrix[[1]][[materialind]][[dind]][[2]][[1]][[1]]
extracthiloETR[matrix_,materialind_,dind_,final_,initial_]:=matrix[[1]][[materialind]][[dind]][[2]][[final]][[2]]-matrix[[1]][[materialind]][[dind]][[2]][[initial]][[1]]
extracthiloF0[calcmatrix_,inmatrix_,materialind_,dind_,final_,initial_]:=2*inmatrix[[materialind]][[1]]*extracthiloETR[calcmatrix,materialind,dind,final,initial]*((1/2)*NIntegrate[(r^2)*(calcmatrix[[2]][[materialind]][[dind]][[2]][[final]][[2]]/.\[Xi]->ArcTan[r])*(calcmatrix[[2]][[materialind]][[dind]][[2]][[initial]][[1]]/.\[Xi]->ArcTan[r]),{r,0,10^4},MaxRecursion->100])^2
extractmaxconc[calcmatrix_,inmatrix_,materialind_,dind_,nind_]:=0.02/((((NIntegrate[(r^2)*(calcmatrix[[2]][[materialind]][[dind]][[2]][[nind]][[1]]/.\[Xi]->ArcTan[r])*(calcmatrix[[2]][[materialind]][[dind]][[2]][[nind]][[1]]/.\[Xi]->ArcTan[r]),{r,0,10^4},MaxRecursion->100]*(2/3))^2)))
extractradius[calcmatrix_,inmatrix_,materialind_,dind_,n_,l_]:=
Sqrt[
NIntegrate[
(r^3)*
(calcmatrix[[2]][[materialind]][[dind]][[2]][[n]][[l]]/.\[Xi]->ArcTan[r])*
(calcmatrix[[2]][[materialind]][[dind]][[2]][[n]][[l]]/.\[Xi]->ArcTan[r]),
{r,0,10^4},MaxRecursion->100
]
]
extractalpha[calcmatrix_,inmatrix_,materialind_,dind_,final_,initial_]:=2*((2*Pi^2)/(Sqrt[4.9]*(137)))*(na/(inmatrix[[materialind]][[1]]*(2*inmatrix[[materialind]][[3]])))*extracthiloF0[calcmatrix,inmatrix,materialind,dind,final,initial]*(2/damp)
extractAFAC[calcmatrix_,inmatrix_,materialind_,dind_,final_,initial_]:= 1 -Exp[-extractalpha[calcmatrix,inmatrix,materialind,dind,final,initial]*2*inmatrix[[materialind]][[3]]]


(* ::Input:: *)
(*alphascale[inmatrix_,materialind_]:=2*((2*Pi^2)/(Sqrt[4.9]*(137)))*(na/(inmatrix[[materialind]][[1]]*(2*inmatrix[[materialind]][[3]])))*(2/damp)*)
(*afacscale[inmatrix_,materialind_]:=1 -Exp[-extractalpha[calcmatrix,inmatrix,materialind,dind,final,initial]*2*inmatrix[[materialind]][[3]]]*)



