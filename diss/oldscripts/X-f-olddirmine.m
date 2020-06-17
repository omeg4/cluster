getdireb[matrix_, type_, eperp_] := 
 matrix[[2]][[eperp]][[type]][[4]][[1]][[1]]
getdirmu[matrix_, type_, eperp_] := matrix[[2]][[eperp]][[type]][[3]]
getdirgap[matrix_, type_, eperp_] := 
 matrix[[2]][[eperp]][[type]][[2]]
getdireperp[matrix_, eperp_] := matrix[[2]][[eperp]][[1]]
getdirEVs[matrix_, type_, eperp_] := 
 matrix[[2]][[eperp]][[type]][[4]]
getdirEFs[matrix_, type_, eperp_] := 
 matrix[[2]][[eperp]][[type]][[5]]
calcdirf0[matrix_, type_, eperp_] := 
 2*getdirmu[matrix, type, 
   eperp]*((getdirEVs[matrix, type, eperp][[2]][[2]] - 
      getdirEVs[matrix, type, eperp][[1]][[1]])/
    H2eV)*((1/2)*
     NIntegrate[(r^2)*getdirEFs[matrix, type, eperp][[2]][[2]]*
       getdirEFs[matrix, type, eperp][[1]][[1]], {r, 0, 10^6}, 
      MaxRecursion -> 20])^2
(* I think I need to add a factor of 2 here because direct excitons \
live in one monolayer *)

calcdiralpha[matrix_, pars_, type_, eperp_] := 
 2*((4*Pi)/(Sqrt[4.89]*(137)))*(na/((pars[[4]]/B2nm)*
      getdirmu[matrix, type, eperp]))*
  calcdirf0[matrix, type, eperp]*(2/damp)
calcdirafac[matrix_, pars_, type_, eperp_] := 
 1 - Exp[-calcdiralpha[matrix, pars, type, eperp]*(pars[[4]]/B2nm)]

