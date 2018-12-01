SetDirectory["/home/mbrunetti/cluster/phos/results/VicCxact-hr-2018-11-29"]
kappa = 4.89;
Export["inps.txt",{VicCxact-hr,mus,chiphos,kappa,{50,50},VCoul[4.89],10^-5}];
result=Table[
	CompPhos[10,{mus[[i]][[1]],mus[[i]][[2]],chiphos,kappa},nhbn,VCoul[4.89],10^-5],
	{i,1},{nhbn,50,50}
	];
Export["suite.m",result];
Export["proc.m",ProcessPhosInd[result]];
Quit[]
