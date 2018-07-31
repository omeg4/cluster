SetDirectory["/home/mbrunetti/cluster/tmdc/results/A-em-hw-hbn-IND-2018-07-05"]
{me,mh,chi2d,kappa,d}=Import["inputs.m"]
Export["suite.m",IndKeld[3,{me,mh,chi2d},kappa,d]];
Quit[]
