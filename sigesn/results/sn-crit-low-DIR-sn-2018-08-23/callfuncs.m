SetDirectory["/home/mbrunetti/cluster/sigesn/results/sn-crit-low-DIR-sn-2018-08-23"]
params={50.5,0.085,550000,0.5,24};
etab={0.58,0.62,0.02};
Export["inp.m",{params,1,etab,{0,0,1}}]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",SiGeSuite[3,params,etab,1]];
"Suite run complete. Processing data.">>>"diag.txt"
labels={"min","max"};
Export["proc.m",ProcessSuite[]];
"Processing complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
