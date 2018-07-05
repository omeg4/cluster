SetDirectory["/home/mbrunetti/cluster/tmdc/results/B-vac-DIR-wse2-2018-07-04"]
params={0.3,0.36,7.571};
Export["inp.m",params,1]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,1]];
Quit[]
