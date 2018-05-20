SetDirectory["/home/mbrunetti/cluster/sigesn/results/DIR-ge-2018-05-15-imitate"]
params={33,0.0676,620000,0.45,16};
etab={0,2.7,0.025};
Export["diag1.txt","Params and Etab initialized"]
SiGeSuite[3,params,etab,4.89];
Export["suitediag.txt","Suite run complete"]
Quit[]
