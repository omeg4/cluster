SetDirectory["/home/mbrunetti/cluster/sigesn/results/fssi-DIR-si-2018-06-05"]
params={0.95,0.046,650000,0.4,11.9};
etab={0,2.7,0.1};
Export["inp.m",{params,1,etab,{0,0,1}}]
Export["diag1.txt","Params and Etab initialized"]
Export["suite.m",SiGeSuite[3,params,etab,1]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, beginning analysis of data files."]
Export["proc.m",ProcessSuite[]];
Quit[]
