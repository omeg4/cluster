Nintfromfile[ef1_,ef2_]:=NIntegrate[r*ef1*ef2, {r,0,10^6},MinRecursion->10,MaxRecursion->50]

f0fromfile[filedata_,nf_]:=Module[{mu=filedata[[1]][[2]],gsev=filedata[[2]][[1]][[1]],exev=filedata[[2]][[nf]][[2]],gsef=filedata[[3]][[1]][[1]],exef=filedata[[3]][[nf]][[2]]},2*mu*(QuantityMagnitude@UnitConvert[exev-gsev,"Hartrees"])*Nintfromfile[gsef,exef]]

absfromfile[filedata_,nf_,mu_,hl_,kappa_,na_,damp_]:=2*((4*\[Pi])/(Sqrt[kappa]*(137)))*(na/((hl/B2nm)*mu))*f0fromfile[filedata,nf]*(2/damp)
