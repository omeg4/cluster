SetDirectory["/home/mbrunetti/cluster/sigesn/results/fssi-DIR-si-2018-06-05"]
Export["proc.m",ProcessSuite[]];
"Reproc complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
