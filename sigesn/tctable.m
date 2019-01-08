(megatclistplot =
    Module[
     {
      maxde = ParallelTable[
         Quiet[
          mc[conf, meth, "tcmax"][2, ep, 1, nx, 
            QConc[5*10^15]][[1]]],
         {mc, mcs}, {conf, configs}, {meth, methods}, {nx, 4}, {ep, 
          mc["Emax"]}
         ] // AbsoluteTiming,
      datatime,
      datadata
      },
     {datatime, datadata} = ParallelTable[
       {
        FPeperp[mc[[2]]["prc"], ep],
        mc[[2]][conf[[2]], meth[[2]], "tcbkt"][2, ep, 1, nx, 
         QConc[conc], 
         maxde[[2]][[mc[[1]], conf[[1]], meth[[1]], nx, ep]]]
        },
       {mc, MapIndexed[{#2[[1]], #1} &, mcs]},
       {conf, MapIndexed[{#2[[1]], #1} &, configs]},
       {meth, MapIndexed[{#2[[1]], #1} &, methods]},
       {freq, {10^13, 2.5*10^13, 5*10^13}},
       {nx, 4},
       {conc, Table[i*10^15, {i, {1, 3, 5}}]},
       {ep, 
        mc[[2]][conf[[2]], meth[[2]], "SCstart"][2, 1, nx, 
         QFreq[freq]], mc[[2]]["Emax"]}
       ];
     {{maxde[[1]], datatime}, {maxde[[2]], datadata}}
     ];) // AbsoluteTiming
megatclistplot[[1]]
megatclistplot[[2, 2]] // Dimensions
Export["0107-megatclistplot-mapindexed.m", megatclistplot[[2, 2]]]
