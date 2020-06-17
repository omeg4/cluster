getindeb[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[2]][[1]][[1]]
getindmu[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[1]][[2]]
getindeperp[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[1]]
getindgap[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[1]][[1]]
getindd[matrix_, type_, dind_] := matrix[[2]][[type]][[dind]][[1]]
getindEVs[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[2]]
getindEFs[matrix_, type_, dind_, eperp_] := 
 matrix[[2]][[type]][[dind]][[2]][[eperp]][[2]][[3]]
calcindf0[matrix_, type_, dind_, eperp_] := 
 2*getindmu[matrix, type, dind, 
   eperp]*((getindEVs[matrix, type, dind, eperp][[2]][[2]] - 
      getindEVs[matrix, type, dind, eperp][[1]][[1]])/
    H2eV)*((1/2)*
     NIntegrate[(r^2)*(getindEFs[matrix, type, dind, eperp][[2]][[
          2]] /. \[Xi] -> 
          ArcTan[r])*(getindEFs[matrix, type, dind, eperp][[1]][[
          1]] /. \[Xi] -> ArcTan[r]), {r, 0, 10^6}, 
      MaxRecursion -> 20])^2
calcindalpha[matrix_, pars_, type_, dind_, eperp_] := 
 2*((2*Pi)/(Sqrt[4.89]*(137)))*(na/(pars[[4]]*
      getindmu[matrix, type, dind, eperp]))*
  calcindf0[matrix, type, dind, eperp]*(2/damp)
(* Why wasn't I multiplying by 2l in the argument of the exponential???? \
that worries me. Need to dig though all my old calculations for \
TMDC's now.... *)
calcindafac[matrix_, pars_, type_, dind_, eperp_] := 
 1 - Exp[-calcindalpha[matrix, pars, type, dind, eperp]*2*pars[[4]]]

getTindeb[matrix_, type_, dind_, eperp_] := 
 matrix[[type]][[2]][[dind]][[2]][[eperp]][[4]][[1]][[1]][[1]] // 
  QuantityMagnitude
getTindmu[matrix_, type_, dind_, eperp_] := 
 matrix[[type]][[2]][[dind]][[2]][[eperp]][[3]] // QuantityMagnitude
getTindeperp[matrix_, type_, dind_, eperp_] := 
 matrix[[type]][[2]][[dind]][[2]][[eperp]][[1]] // QuantityMagnitude
getTindgap[matrix_, type_, dind_, eperp_] := 
 matrix[[type]][[2]][[dind]][[2]][[eperp]][[2]] // QuantityMagnitude
getTindd[matrix_, type_, dind_] := 
 matrix[[type]][[2]][[dind]][[1]] // QuantityMagnitude
getTindEVs[matrix_, type_, dind_, 
  eperp_] := (1/1000)*
   Table[matrix[[type]][[2]][[dind]][[2]][[eperp]][[4]][[i]][[j]][[
     1]], {i, 3}, {j, i}] // QuantityMagnitude
getTindEFs[matrix_, type_, dind_, eperp_] := 
 Table[matrix[[type]][[2]][[dind]][[2]][[eperp]][[4]][[i]][[j]][[
   2]], {i, 3}, {j, i}]
calcTindf0[matrix_, type_, dind_, eperp_, f_] := 
 2*getTindmu[matrix, type, dind, 
   eperp]*((getTindEVs[matrix, type, dind, eperp][[f]][[2]] - 
      getTindEVs[matrix, type, dind, eperp][[1]][[1]])/
    H2eV)*((1/2)*
     NIntegrate[(r^2)*getTindEFs[matrix, type, dind, eperp][[f]][[2]]*
       getTindEFs[matrix, type, dind, eperp][[1]][[1]], {r, 0, 10^6}, 
      MaxRecursion -> 20])^2
calcTindalpha[matrix_, pars_, type_, dind_, eperp_, f_] := 
 2*((2*Pi)/(Sqrt[4.89]*(137)))*(na/((pars[[4]]/B2nm)*
      getTindmu[matrix, type, dind, eperp]))*
  calcTindf0[matrix, type, dind, eperp, f]*(2/damp)
(* Why wasn't I multiplying by 2l in the argument of the exponential???? \
that worries me. Need to dig though all my old calculations for \
TMDC's now.... *)

calcTindafac[matrix_, pars_, type_, dind_, eperp_, f_] := 
 1 - Exp[-calcTindalpha[matrix, pars, type, dind, eperp, 
      f]*2*(pars[[4]]/B2nm)]
