SetDirectory["/home/mbrunetti/cluster/sigesn/results/small-sihbn-DIR-si-2018-05-27"]
params={13.5,0.046,433000,0.333,11.9};
etab={0,2.7,0.1};
Export["inp.m",{params,4.89,etab,{0,0,1}}]
Export["diag1.txt","Params and Etab initialized"]
Export["suite.m",SiGeSuite[3,params,etab,4.89]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, beginning analysis of data files."]
Export["proc.m",ProcessSuite[]];
Quit[]
