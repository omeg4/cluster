SetDirectory["/home/mbrunetti/cluster/sigesn/results/drf-DIR-si-2018-05-21"]
params={1.9,0.046,650000,0.4,11.9};
etab={0,1,0.5};
Export["inp.m",{params,4.89,etab,{0,0,1}}]
Export["diag1.txt","Params and Etab initialized"]
Export["suite.m",SiGeSuite[3,params,etab,4.89]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, beginning analysis of data files."]
Export["proc.m",ProcessSuite[]];
Quit[]
