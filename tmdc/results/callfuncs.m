SetDirectory["/home/mbrunetti/cluster/tmdc/results"]
Emass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"];
Hmass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[ - /2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"];
params={Emass,Hmass,};
Export["inp.m",params,]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,]];
Quit[]
