SetDirectory["/home/mbrunetti/cluster/sigesn/results/fsge-DIR-ge-2018-06-05"]
Export["proc.m",ProcessSuite[]];
"Reproc complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
