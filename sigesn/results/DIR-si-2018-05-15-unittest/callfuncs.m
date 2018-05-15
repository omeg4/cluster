SetDirectory["/home/mbrunetti/cluster/sigesn/results/DIR-si-2018-05-15-unittest"]
params={50,0.1,650000,0.4,11.9};
etab={0,1,1};
Export["diag1.txt","Params and Etab initialized"]
SiGeSuite[3,params,etab,4.89];
Export["suitediag.txt","Suite run complete"]
Quit[]
