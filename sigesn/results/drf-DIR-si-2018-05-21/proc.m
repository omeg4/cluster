(* Created with the Wolfram Language : www.wolfram.com *)
{{1.9, 0.046, 650000, 0.4, 11.9}, 
 {{"min", Table[{dfromsuite[data$9672[[Dn]][[1]][[type]]], 
     Table[{Quantity[eperpfromsuite[data$9672[[Dn]][[En]][[type]]], 
        "Volts"/"Angstroms"], Quantity[egapfromsuite[data$9672[[Dn]][[En]][[
          type]]], "Millielectronvolts"], 
       Quantity[mufromsuite[data$9672[[Dn]][[En]][[type]]], "ElectronMass"], 
       Table[{Quantity[1000*H2eV*EVfromsuite[data$9672[[Dn]][[En]][[type]], 
            n, l], "Millielectronvolts"], Quantity[
          B2nm*r2fromsuite[data$9672[[Dn]][[En]][[type]], n, l], 
          "BohrRadius"]}, {n, 3}, {l, 0, n - 1}], 
       Table[{Quantity[1000*H2eV*Etrfromsuite[data$9672[[Dn]][[En]][[type]], 
            1, nf, 0, 1], "Millielectronvolts"], f0fromsuite[
          data$9672[[Dn]][[En]][[type]], 1, nf, 0, 1], 
         UnitConvert[Quantity[absfromsuite[data$9672[[Dn]][[En]][[type]], 
            params$9672, kappa$9672, 1, nf, 0, 1], 1/"BohrRadius"], 
          1/"Meters"], afacfromsuite[data$9672[[Dn]][[En]][[type]], 
          params$9672, kappa$9672, 1, nf, 0, 1]}, {nf, {2, 3}}]}, {En, Ne}]}, 
    {Dn, {}[[1]]}]}, {"max", Table[{dfromsuite[data$9672[[Dn]][[1]][[type]]], 
     Table[{Quantity[eperpfromsuite[data$9672[[Dn]][[En]][[type]]], 
        "Volts"/"Angstroms"], Quantity[egapfromsuite[data$9672[[Dn]][[En]][[
          type]]], "Millielectronvolts"], 
       Quantity[mufromsuite[data$9672[[Dn]][[En]][[type]]], "ElectronMass"], 
       Table[{Quantity[1000*H2eV*EVfromsuite[data$9672[[Dn]][[En]][[type]], 
            n, l], "Millielectronvolts"], Quantity[
          B2nm*r2fromsuite[data$9672[[Dn]][[En]][[type]], n, l], 
          "BohrRadius"]}, {n, 3}, {l, 0, n - 1}], 
       Table[{Quantity[1000*H2eV*Etrfromsuite[data$9672[[Dn]][[En]][[type]], 
            1, nf, 0, 1], "Millielectronvolts"], f0fromsuite[
          data$9672[[Dn]][[En]][[type]], 1, nf, 0, 1], 
         UnitConvert[Quantity[absfromsuite[data$9672[[Dn]][[En]][[type]], 
            params$9672, kappa$9672, 1, nf, 0, 1], 1/"BohrRadius"], 
          1/"Meters"], afacfromsuite[data$9672[[Dn]][[En]][[type]], 
          params$9672, kappa$9672, 1, nf, 0, 1]}, {nf, {2, 3}}]}, {En, Ne}]}, 
    {Dn, {}[[1]]}]}}}
