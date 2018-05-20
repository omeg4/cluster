SetDirectory["/home/mbrunetti/cluster/sigesn/results/DIR-sn-2018-05-15-imitate"]
params={101,0.085,550000,0.5,24};
etab={0,2.7,0.025};
Export["diag1.txt","Params and Etab initialized"]
SiGeSuite[3,params,etab,4.89];
Export["suitediag.txt","Suite run complete"]
Quit[]
