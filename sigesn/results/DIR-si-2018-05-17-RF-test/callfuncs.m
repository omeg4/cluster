SetDirectory["/home/mbrunetti/cluster/sigesn/results/DIR-si-2018-05-17-RF-test"]
params={1.9,0.046,650000,0.4,11.9};
etab={0,2.5,0.5};
Export["diag1.txt","Params and Etab initialized"]
SiGeSuite[3,params,etab,4.89];
Export["suitediag.txt","Suite run complete"]
Quit[]
