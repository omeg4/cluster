SetDirectory["/home/mbrunetti/cluster/tmdc/results/hanbicki-vac-DIR-ws2-2018-07-30"]
Emass = 0.44
Hmass = 0.45
params={Emass,Hmass,6.03};
Export["inp.m",{params,1}//Flatten]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,1]];
Quit[]
