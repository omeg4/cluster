SetDirectory["/home/mbrunetti/cluster/sigesn/results/hdcrit-DIR-ge-2018-08-24"]
params={16.5,0.0676,620000,0.45,16};
etab={0.020,0.029,0.001};
Export["inp.m",{params,1,etab,{0,0,1}}]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",SiGeSuite[3,params,etab,1]];
"Suite run complete. Processing data.">>>"diag.txt"
labels={"min","max"};
Export["proc.m",ProcessSuite[]];
"Processing complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
