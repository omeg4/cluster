SetDirectory["/home/mbrunetti/cluster/tmdc/results/B-em-hw-sub-keld-IND-2018-07-07"]
{me,mh,chi2d,kappa,d}=Import["inputs.m"]
Export["suite.m",IndKeld[3,{me,mh,chi2d},kappa,d]];
Quit[]
