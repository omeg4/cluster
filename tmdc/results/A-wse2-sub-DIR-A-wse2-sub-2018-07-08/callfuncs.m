SetDirectory["/home/mbrunetti/cluster/tmdc/results/A-wse2-sub-DIR-A-wse2-sub-2018-07-08"]
params={0.4,0.53,7.571};
Export["inp.m",params,2.3]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.3]];
Quit[]
