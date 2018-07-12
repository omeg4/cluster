SetDirectory["/home/mbrunetti/cluster/sigesn/results/hdencsi-DIR-si-2018-06-26"]
params={19,0.046,506000,0.333,11.9};
etab={0,0.2,0.01};
Export["inp.m",{params,4.89,etab,{0,0,1}}]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",SiGeSuite[3,params,etab,4.89]];
"Suite run complete. Processing data.">>>"diag.txt"
labels={"min","max"};
Export["proc.m",ProcessSuite[]];
"Processing complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
