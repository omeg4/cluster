SetDirectory["/home/mbrunetti/cluster/tmdc/results/donck-MB-sio2-DIR-mose2-2018-07-18"]
params={0.54,0.5421,8.23};
Export["inp.m",params,2.4]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.4]];
Quit[]
