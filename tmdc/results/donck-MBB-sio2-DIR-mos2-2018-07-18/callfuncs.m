SetDirectory["/home/mbrunetti/cluster/tmdc/results/donck-MBB-sio2-DIR-mos2-2018-07-18"]
params={0.5,0.5358,6.6};
Export["inp.m",params,2.4]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.4]];
Quit[]
