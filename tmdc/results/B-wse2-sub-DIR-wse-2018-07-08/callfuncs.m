SetDirectory["/home/mbrunetti/cluster/tmdc/results/B-wse2-sub-DIR-wse-2018-07-08"]
params={0.3,0.36,7.571};
Export["inp.m",params,2.3]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.3]];
Quit[]
