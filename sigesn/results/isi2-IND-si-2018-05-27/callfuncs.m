params={0.95,0.046,650000,0.4,11.9};
etab={0,2.75,0.25};
dtab={1,5,1};
Export["diag1.txt","Params, Dtab, and Etab initialized"];
Export["inp.m",{params,4.89,etab,dtab}];
Export["suite.m",SGSIndSuite[params,4.89,etab,dtab,IndCoul]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, processing..."];
Export["proc.m",ProcessSuite[]];
Quit[]
