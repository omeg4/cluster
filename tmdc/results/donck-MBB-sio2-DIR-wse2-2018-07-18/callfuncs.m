SetDirectory["/home/mbrunetti/cluster/tmdc/results/donck-MBB-sio2-DIR-wse2-2018-07-18"]
params={0.34,0.4494,7.18};
Export["inp.m",params,2.4]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.4]];
Quit[]
