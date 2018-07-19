SetDirectory["/home/mbrunetti/cluster/tmdc/results/donck-MB-vac-DIR-wse2-2018-07-18"]
params={0.34,0.3364,7.18};
Export["inp.m",params,1]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,1]];
Quit[]
