types={"MoS2-hi","MoS2-lo","MoSe2-hi","MoSe2-lo","WS2-hi","WS2-lo","WSe2-hi","WSe2-lo"};
Export["/home/mbrunetti/tmdc/results/2018-03-04_tmdc-coul-2/2018-03-04_tmdc-coul-2-eigenenergies.m",Table[{types[[t]],Table[{nhbn-2,thedata[[1]][[t]][[nhbn]][[2]]},{nhbn,3,11}]},{t,8}]]
(* Extract r^2 *)
Export["/home/mbrunetti/tmdc/results/2018-03-04_tmdc-coul-2/2018-03-04_tmdc-coul-2-radius.m",Table[{types[[t]],Table[{nhbn-2,extractradius[thedata,EXmatrange,t,nhbn,1,1]},{nhbn,3,11}]},{t,8}]]
(* Extract f0, 1\[Rule]2 *)
Export["/home/mbrunetti/tmdc/results/2018-03-04_tmdc-coul-2/2018-03-04_tmdc-coul-2-f012.m",Table[{types[[t]],Table[{nhbn-2,extracthiloF0[thedata,EXmatrange,t,nhbn,2,1]},{nhbn,3,11}]},{t,8}]]
(* Extract f0, 1\[Rule]3 *)
Export["/home/mbrunetti/tmdc/results/2018-03-04_tmdc-coul-2/2018-03-04_tmdc-coul-2-f013.m",Table[{types[[t]],Table[{nhbn-2,extracthiloF0[thedata,EXmatrange,t,nhbn,3,1]},{nhbn,3,11}]},{t,8}]]
(* Extract \[Alpha] 1\[Rule]2 *)
Export["/home/mbrunetti/tmdc/results/2018-03-04_tmdc-coul-2/2018-03-04_tmdc-coul-2-alpha12.m",Table[{types[[t]],Table[{nhbn-2,extractalpha[thedata,EXmatrange,t,nhbn,2,1]},{nhbn,3,11}]},{t,8}]]
(* Extract \[Alpha] 1\[Rule]3 *)
Export["/home/mbrunetti/tmdc/results/2018-03-04_tmdc-coul-2/2018-03-04_tmdc-coul-2-alpha13.m",Table[{types[[t]],Table[{nhbn-2,extractalpha[thedata,EXmatrange,t,nhbn,3,1]},{nhbn,3,11}]},{t,8}]]
(* Extract \[ScriptCapitalA] 1\[Rule]2 *)
Export["/home/mbrunetti/tmdc/results/2018-03-04_tmdc-coul-2/2018-03-04_tmdc-coul-2-afac12.m",Table[{types[[t]],Table[{nhbn-2,extractAFAC[thedata,EXmatrange,t,nhbn,2,1]},{nhbn,3,11}]},{t,8}]]
(* Extract \[ScriptCapitalA] 1\[Rule]3 *)
Export["/home/mbrunetti/tmdc/results/2018-03-04_tmdc-coul-2/2018-03-04_tmdc-coul-2-afac13.m",Table[{types[[t]],Table[{nhbn-2,extractAFAC[thedata,EXmatrange,t,nhbn,3,1]},{nhbn,3,11}]},{t,8}]]
