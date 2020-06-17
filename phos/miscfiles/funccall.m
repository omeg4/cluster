SetDirectory["/home/mbrunetti/cluster/phos"]
Export["inps.txt",{,,,,,,,}];
result=CompPhos[,{,,,},,,];
Export["suite.m",result]
Quit[]
