SetDirectory["/home/mbrunetti/cluster/sigesn/results/small-sihbn-DIR-si-2018-05-27"]
Export["proc.m",ProcessSuite[]];
"Reproc complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
