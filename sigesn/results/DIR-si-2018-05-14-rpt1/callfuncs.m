SetDirectory["/home/mbrunetti/cluster/sigesn/results/DIR-si-2018-05-14-rpt1"]
fpath=ToString[StringForm["/home/mbrunetti/cluster/sigesn/results/`1`/",DIR-si-2018-05-14-rpt1]];
params={3.95,0.047,542000,0.35,11.9};
etab={0,1,0.25};
Export["diag1.txt","Params and Etab initialized"]
suite=SiGeSuite[3,params,etab,4.89];
Export["suitediag.txt","Suite run complete"]
Export["results.m",suite];
Quit[]
