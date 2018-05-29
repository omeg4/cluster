params={16.5,0.0676,620000,0.45,16};
etab={0,2.75,0.25};
dtab={1,6,1};
Export["diag1.txt","Params, Dtab, and Etab initialized"];
Export["inp.m",{params,4.89,etab,dtab}];
Export["suite.m",SGSIndSuite[params,4.89,etab,dtab,IndCoul]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, processing..."];
Export["proc.m",ProcessSuite[]];
Quit[]
