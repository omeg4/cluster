SetDirectory["/home/mbrunetti/cluster/tmdc/results/B-em-hw-vac-IND-2018-07-04"]
{me,mh,chi2d,kappa,d}=Import["../../inputs.m"]
Export["inp.m",{me,mh,chi2d,kappa}];
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",IndKeld[3,{me,mh,chi2d},kappa,d]
Export["suite.m",suite];
Quit[]
