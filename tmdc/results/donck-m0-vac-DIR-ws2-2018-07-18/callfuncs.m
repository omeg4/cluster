SetDirectory["/home/mbrunetti/cluster/tmdc/results/donck-m0-vac-DIR-ws2-2018-07-18"]
params={0.32,0.32,6.03};
Export["inp.m",params,1]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,1]];
Quit[]
