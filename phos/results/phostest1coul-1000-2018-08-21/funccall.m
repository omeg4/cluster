SetDirectory["/home/mbrunetti/cluster/phos/results/phostest1coul-1000-2018-08-21"]
Export["inps.txt",{phostest1coul,"C",1,1,1,1,0,1000}];
result=CompPhos[4,{1,1,1,1},0,"C",1000];
Export["suite.m",result]
Quit[]
