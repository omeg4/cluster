SetDirectory["/home/mbrunetti/cluster/tmdc/results/B-em-hw-hbn-IND-2018-07-05"]
{me,mh,chi2d,kappa,d}=Import["../../inputs.m"]
Export["inp.m",{me,mh,chi2d,kappa, d}];
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",IndKeld[3,{me,mh,chi2d},kappa,d]];
Quit[]
