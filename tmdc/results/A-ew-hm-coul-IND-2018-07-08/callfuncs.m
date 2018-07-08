SetDirectory["/home/mbrunetti/cluster/tmdc/results/A-ew-hm-coul-IND-2018-07-08"]
{me,mh,chi2d,kappa,d}=Import["inputs.m"]
Export["suite.m",IndCoul[3,{me,mh,chi2d},kappa,d]];
Quit[]
