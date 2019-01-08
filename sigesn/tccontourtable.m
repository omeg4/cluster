(cplt4darray = Flatten[
     ParallelTable[
      With[
       {scut = mc[conf, meth, "SCstart"][2, 1, nx, QFreq[gamex]]},
       ListContourPlot[
        Flatten[
         Table[
          {
           FPeperp[mc["prc"], ep],
           de,
           Quiet@mc[conf, meth, "tcbkt"][2, ep, 1, nx, QConc[10^15],
             de]
           },
          {de, -0.1, 0.1, 0.01}, {ep, scut, mc["Emax"]}
          ],
         1],
        FrameLabel -> {{"\[CapitalDelta]E", 
           None}, {"\!\(\*SubscriptBox[\(E\), \(\[Perpendicular]\)]\) \
[V/\[Angstrom]]", None}},
        PlotLabel -> bsmc["mono", "old", "fullname"],
        LabelStyle -> Directive[24, Black, FontFamily -> "Arial"],
        ImageSize -> {500, 500},
        ColorFunction -> Function[{z}, Hue[Mod[z/300, 1]]],
        ColorFunctionScaling -> False,
        PlotLegends -> Automatic,
        Contours -> Table[i*10, {i, 100}],
        ContourLabels -> Automatic,
        PlotRange -> {{-0.05, 2.75}, {-0.105, 0.105}},
        InterpolationOrder -> 0,
        RegionFunction -> 
         Function[{x, y, z}, 
          QuantityMagnitude[FPeperp[bsmc["prc"], scut]] < x]
        ]
       ],
      {mc, mcs}, {conf, configs}, {meth, methods}, {nx, 
       4}, {gamex, {10^13, 2.5*10^13, 5*10^13}}
      ],
     {{1}, {2, 3}, {4}, {5}}
     ];) // AbsoluteTiming
Export["0107-bigcontourplotarray.m", cplt4darray]
