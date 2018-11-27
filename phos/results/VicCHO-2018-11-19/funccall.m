SetDirectory["/home/mbrunetti/cluster/phos/results/VicCHO-2018-11-19"]
kappa = 4.89;
Export["inps.txt",{VicCHO,mus,chiphos,kappa,{1,5},VCho[4.89],10^-3}];
result=Table[
	CompPhos[10,{mus[[i]][[1]],mus[[i]][[2]],chiphos,kappa},nhbn,VCho[4.89],10^-3],
	{i,4},{nhbn,1,5}
	];
Export["suite.m",result];
Export["proc.m",ProcessPhosInd[result]];
Quit[]
