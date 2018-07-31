SetDirectory["/home/mbrunetti/cluster/tmdc/results/vdd-A-vac-DIR-mose2-2018-07-19"]
Emass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.47 - 0.18/2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.313,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[0.94,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
Hmass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.47 - 0.18/2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.313,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[0.94,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
params={Emass,Hmass,8.23};
Export["inp.m",{params,1}//Flatten]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,1]];
Quit[]
