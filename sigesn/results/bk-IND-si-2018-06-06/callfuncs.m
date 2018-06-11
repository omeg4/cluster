SetDirectory["/home/mbrunetti/cluster/sigesn/results/bk-IND-si-2018-06-06"]
params={19,0.046,506000,0.4,11.9};
etab={0,2.75,0.25};
dtab={1,5,1};
Export["inp.m",{params,4.89,etab,dtab}];
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",SGSIndSuite[params,4.89,etab,dtab,IndKeld]];
"Suite run complete. Processing data.">>>"diag.txt"
labels={"min","max"};
Export["proc.m",ProcessSuite[]];
"Processing complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
