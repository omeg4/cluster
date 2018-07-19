SetDirectory["/home/mbrunetti/cluster/tmdc/results/vdd-B-sio2-DIR-mos2-2018-07-18"]
Emass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.66,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.193,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[1.1,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
Hmass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.66 - -0.15/2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.193,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[1.1,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
params={Emass,Hmass,6.6};
Export["inp.m",{params,2.4}//Flatten]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.4]];
Quit[]
