SetDirectory["/home/mbrunetti/cluster/tmdc/results/hanbicki-DIR-ws2-2018-07-29"]
Emass = 0.44
Hmass = 0.45
params={Emass,Hmass,6.03};
Export["inp.m",{params,2.45}//Flatten]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.45]];
Quit[]
