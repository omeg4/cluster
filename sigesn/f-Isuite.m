SGSIndSuite[pars_, ezrange_, drange_, name_, projdir_] := Module[
  {make,
   normed,
   processed,
   ei = ezrange[[1]],
   ef = ezrange[[2]],
   estep = ezrange[[3]],
   di = drange[[1]],
   df = drange[[2]],
   dstep = drange[[3]]
   },
  make = {
    Table[{dd, 
      Table[{Eperp, 
        IndKeld[2, pars[[2 ;;]], pars[[1]], Eperp, 4.89, 
         dd*lBN + (pars[[4]]/B2nm), 1]}, {Eperp, ei, ef, 
        estep}]}, {dd, di, df, dstep}],
    Table[{dd, 
      Table[{Eperp, 
        IndKeld[2, pars[[2 ;;]], pars[[1]], Eperp, 4.89, 
         dd*lBN + (pars[[4]]/B2nm), -1]}, {Eperp, ei, ef, 
        estep}]}, {dd, di, df, dstep}]
    };
  normed = normindEFs[make, 10^5];
  processed = {
    pars,
    Table[
     {
      labels[[type]],
      Table[
       {
        Nindeperp[normed, type, 1, eperp],
        Nindgap[normed, type, 1, eperp],
        Nindmu[normed, type, 1, eperp],
        Table[
         {
          NindD[normed, type, d],
          Table[
           {
            NindEV[normed, type, d, eperp, n, l],
            Sqrt[NindNINT[normed, type, d, eperp, n, l, n, l, 3]]
            },
           {n, 3}, {l, 0, n - 1}
           ],
          Table[
           {
            (NindEV[normed, type, d, eperp, nf, 1] - 
              NindEV[normed, type, d, eperp, 1, 0]),
            Nindf0[normed, type, d, eperp, 1, 0, nf, 1],
            Nindalpha[normed, type, d, eperp, 1, 0, nf, 1, pars],
            Nindafac[normed, type, d, eperp, 1, 0, nf, 1]
            },
           {nf, {2, 3}}
           ]
          },
         {d, di, df, dstep}
         ]
        },
       {eperp, 2, ((ef - ei)/estep) + 1}
       ]
      },
     {type, 2}
     ]
    };
  Export[projdir <> 
    ToString[
     StringForm["`1`-`2`-`3`_`4`_normed.m", DateObject[][[1]][[1]], 
      DateObject[][[1]][[2]], DateObject[][[1]][[3]], name]], 
   normed];
  Export[projdir <> 
    ToString[
     StringForm["`1`-`2`-`3`_`4`_processed.m", DateObject[][[1]][[1]],
       DateObject[][[1]][[2]], DateObject[][[1]][[3]], name]], 
   processed];
  Export[projdir <> 
    ToString[
     StringForm["`1`-`2`-`3`_`4`_info.txt", DateObject[][[1]][[1]], 
      DateObject[][[1]][[2]], DateObject[][[1]][[3]], name]],
   ToString[
    StringForm[
      "input parameters: ``", pars];
    ]
	]
]
