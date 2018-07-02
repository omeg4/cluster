SetDirectory["/home/mbrunetti/cluster/sigesn/results/fssi-vf1-DIR-fssi-vf1-2018-07-01"]
params={0.95,0.046,532000,0.4,11.9};
etab={1,2,1};
Export["inp.m",{params,1,etab,{0,0,1}}]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",SiGeSuite[3,params,etab,1]];
"Suite run complete. Processing data.">>>"diag.txt"
labels={"min","max"};
Export["proc.m",ProcessSuite[]];
"Processing complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
