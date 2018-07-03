SetDirectory["/home/mbrunetti/cluster/sigesn/results/vf4-fssi-DIR-si-2018-07-02"]
params={0.95,0.046,575000,0.4,11.9};
etab={0,2,0.5};
Export["inp.m",{params,1,etab,{0,0,1}}]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",SiGeSuite[3,params,etab,1]];
"Suite run complete. Processing data.">>>"diag.txt"
labels={"min","max"};
Export["proc.m",ProcessSuite[]];
"Processing complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
