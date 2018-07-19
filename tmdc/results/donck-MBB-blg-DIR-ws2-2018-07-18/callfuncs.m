SetDirectory["/home/mbrunetti/cluster/tmdc/results/donck-MBB-blg-DIR-ws2-2018-07-18"]
params={0.32,0.3982,6.03};
Export["inp.m",params,2.79]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.79]];
Quit[]
