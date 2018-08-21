SetDirectory["/home/mbrunetti/cluster/phos/results/phostest1coul-10000-2018-08-21"]
Export["inps.txt",{phostest1coul,"C",1,1,1,1,0,10000}];
result=CompPhos[4,{1,1,1,1},0,"C",10000];
Export["suite.m",result]
Quit[]
