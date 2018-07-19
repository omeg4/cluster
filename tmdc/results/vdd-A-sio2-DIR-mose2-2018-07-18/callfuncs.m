SetDirectory["/home/mbrunetti/cluster/tmdc/results/vdd-A-sio2-DIR-mose2-2018-07-18"]
Emass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.47,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.313,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[0.94,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
Hmass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.47 - 0.18/2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.313,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[0.94,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
params={Emass,Hmass,8.23};
Export["inp.m",{params,2.4}//Flatten]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.4]];
Quit[]
