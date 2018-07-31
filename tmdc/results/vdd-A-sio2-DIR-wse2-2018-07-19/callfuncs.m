SetDirectory["/home/mbrunetti/cluster/tmdc/results/vdd-A-sio2-DIR-wse2-2018-07-19"]
Emass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.60 - 0.46/2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.310,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[1.19,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
Hmass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.60 - 0.46/2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.310,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[1.19,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
params={Emass,Hmass,7.18};
Export["inp.m",{params,2.4}//Flatten]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.4]];
Quit[]
