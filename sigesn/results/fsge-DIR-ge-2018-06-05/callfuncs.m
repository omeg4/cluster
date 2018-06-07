SetDirectory["/home/mbrunetti/cluster/sigesn/results/fsge-DIR-ge-2018-06-05"]
params={16.5,0.0676,620000,0.45,16};
etab={0,2.7,0.1};
Export["inp.m",{params,1,etab,{0,0,1}}]
Export["diag1.txt","Params and Etab initialized"]
Export["suite.m",SiGeSuite[3,params,etab,1]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, beginning analysis of data files."]
Export["proc.m",ProcessSuite[]];
Quit[]
