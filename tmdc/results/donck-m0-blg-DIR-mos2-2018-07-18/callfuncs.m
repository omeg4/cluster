SetDirectory["/home/mbrunetti/cluster/tmdc/results/donck-m0-blg-DIR-mos2-2018-07-18"]
params={0.5,0.5,6.6};
Export["inp.m",params,2.79]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.79]];
Quit[]
