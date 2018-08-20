getEV[mux_,muy_,d_,n_]:=result[[mux]][[muy]][[3]][[d]][[2]][[1]][[n]]
getEF[mux_,muy_,d_,n_]:=result[[mux]][[muy]][[3]][[d]][[2]][[2]][[n]]
getmux[mux_]:=result[[mux]][[1]][[1]]
getmuy[muy_]:=result[[1]][[muy]][[2]]
pickmu[mux_,muy_,xy_]:=result[[mux]][[muy]][[xy]]
Lphos=1.117/B2nm
getf0[mux_,muy_,d_,n_,xy_]:=2*pickmu[mux,muy,xy]*(getEV[mux,muy,d,n]-getEV[mux,muy,d,1])*NIntegrate[getEF[mux,muy,d,n]*If[xy==1,x,y]*getEF[mux,muy,d,1],{x,-s,s},{y,-s,s},MaxRecursion->100,WorkingPrecision->100]
(* How to organize processed data? *)
phosprocess={
	chiphos*B2nm*Quantity["Nanometers"],
	Table[
	{
		result[[mux]][[muy]][[1]]*Quantity["ElectronMass"],
		result[[mux]][[muy]][[2]]*Quantity["ElectronMass"],
		{
			"orthonormcheck",
			Table[
				NIntegrate[getEF[mux,muy,d,ni]*getEF[mux,muy,d,nj],{x,-s,s},{y,-s,s},MaxRecursion->100,WorkingPrecision->100],
				{ni,10},{nj,10}
			]	
		},
		{
			"f0 matrix with Etr",
			Table[
				{
					(getEV[mux,muy,d,n]-getEV[mux,muy,d,1])*H2eV*1000,
					getf0[mux,muy,d,n,1],
					getf0[mux,muy,d,n,2]
				},
				{n,2,10}
			]
		}
	},
	{mux,2},{muy,2},{d,10}
	]
};
Export[ToString@StringForm["/home/mbrunetti/phos/results/`1`/processed.m",projdir],phosprocess]
Quit[]
