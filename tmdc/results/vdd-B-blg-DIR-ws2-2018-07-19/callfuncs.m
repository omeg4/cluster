SetDirectory["/home/mbrunetti/cluster/tmdc/results/vdd-B-blg-DIR-ws2-2018-07-19"]
Emass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.79 - -0.43/2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.197,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[1.37,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
Hmass = UnitConvert[(Quantity["ReducedPlanckConstant"]^2)*UnitConvert[Quantity[1.79 - -0.43/2,"Electronvolts"],"Hartrees"]/(2*(UnitConvert[Quantity[3.197,"Angstroms"],"BohrRadius"]^2)*(UnitConvert[Quantity[1.37,"Electronvolts"],"Hartrees"]^2)),"ElectronMass"]//QuantityMagnitude;
params={Emass,Hmass,6.03};
Export["inp.m",{params,2.79}//Flatten]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,2.79]];
Quit[]
