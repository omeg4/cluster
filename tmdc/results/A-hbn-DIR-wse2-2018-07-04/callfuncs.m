SetDirectory["/home/mbrunetti/cluster/tmdc/results/A-hbn-DIR-wse2-2018-07-04"]
params={0.4,0.53,7.571};
Export["inp.m",params,4.89]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,4.89]];
Quit[]
