params={1.9,0.046,650000,0.4,11.9};
etab={0,2.75,0.25};
dtab={1,9,1};
Export["diag1.txt","Params, Dtab, and Etab initialized"];
Export["inp.m",{params,4.89,etab,dtab}];
Export["suite.m",SGSIndSuite[params,4.89,etab,dtab,IndCoul]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, processing..."];
Export["proc.m",ProcessSuite[]];
prc=Import["proc.m"];
Export["muplt.pdf",mkmuplt[prc]];
Export["egapplt.pdf",mkegplt[prc]];
Export["ebplt.pdf",mkebplt[prc]];
Export["f0plt.pdf",mkf0plt[prc]];
Export["absplt.pdf",mkabsplt[prc]];
Export["afacplt.pdf",mkafacplt[prc]];
Quit[]
