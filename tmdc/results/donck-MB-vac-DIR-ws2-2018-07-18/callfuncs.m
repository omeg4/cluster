SetDirectory["/home/mbrunetti/cluster/tmdc/results/donck-MB-vac-DIR-ws2-2018-07-18"]
params={0.32,0.3128,6.03};
Export["inp.m",params,1]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,1]];
Quit[]
