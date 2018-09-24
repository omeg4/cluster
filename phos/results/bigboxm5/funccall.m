SetDirectory["/home/mbrunetti/cluster/phos/results/bigboxm5"]
result=CompPhos[20,{1/4,1/9,1,1},1,VCho[1,1]];
Export["suite.m",result];
Quit[]
