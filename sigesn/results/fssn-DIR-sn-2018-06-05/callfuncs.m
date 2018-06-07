SetDirectory["/home/mbrunetti/cluster/sigesn/results/fssn-DIR-sn-2018-06-05"]
params={50.5,0.085,550000,0.5,24};
etab={0,2.7,0.1};
Export["inp.m",{params,1,etab,{0,0,1}}]
Export["diag1.txt","Params and Etab initialized"]
Export["suite.m",SiGeSuite[3,params,etab,1]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, beginning analysis of data files."]
Export["proc.m",ProcessSuite[]];
Quit[]
