SetDirectory["/home/mbrunetti/cluster/phos/results/VicCHO-lr-2018-11-20"]
kappa = 4.89;
Export["inps.txt",{VicCHO-lr,mus,chiphos,kappa,{1,20},VCho[4.89],10^-3}];
result=Table[
	CompPhos[10,{mus[[i]][[1]],mus[[i]][[2]],chiphos,kappa},nhbn,VCho[4.89],10^-3],
	{i,1},{nhbn,1,20}
	];
Export["suite.m",result];
Export["proc.m",ProcessPhosInd[result]];
Quit[]
