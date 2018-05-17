params={1.9,0.046,650000,0.4,11.9};
etab={0,2.5,0.5};
dtab={1,5,};
Export["diag1.txt","Params, Dtab, and Etab initialized"]
SGSIndSuite[params,4.89,etab,dtab,IndCoul];
Quit[]
