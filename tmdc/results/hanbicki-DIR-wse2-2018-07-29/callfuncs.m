SetDirectory["/home/mbrunetti/cluster/tmdc/results/hanbicki-DIR-wse2-2018-07-29"]
Emass = 0.53
Hmass = 0.52
params={Emass,Hmass,7.18};
Export["inp.m",{params,2.45}//Flatten]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.45]];
Quit[]
