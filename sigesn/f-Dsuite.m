SiGeSuite[mmax_, matpars_, eztab_, kappa_] :=
  Module[
   {
    (* Unpack material parameters and get units straight *)
    Egap = matpars[[1]],
    (* Express bandgaps in Joules *)
    d0 = matpars[[2]]*Unm/Um,
	(* Express buckling constant in meters *)
    vF = matpars[[3]],
	(* Fermi velocity is already in m/s *)
    ds = matpars[[4]]*Unm/Um, (* Monolayer in thickness in meters *)
	erel = matpars[[5]] - 1,
    \[Kappa] = kappa,
	min,
	max
	MUTABMIN={};
	MUTABMAX={};
	EVTABMIN={};
	EVTABMAX={};
	ETRTABMIN={};
	ETRTABMAX={};
	F0TAB={};
    },
    Table[
	Export[ToString@StringForm["calcs_e`1`.m",Eperp],
	{
		Eperp,
		min={"min", DirKeld[mmax, matpars[[2 ;;]], Egap, Eperp, \[Kappa], 1]},
		max={"max", DirKeld[mmax, matpars[[2 ;;]], Egap, Eperp, \[Kappa], -1]}
	}],
	{Eperp, eztab[[1]], eztab[[2]], eztab[[3]]}
	];
	Export[ToString@StringForm["tabledone.txt","Table done, files saved. Starting analysis of data"]];
	resultfiles = FileNames["*calcs*"];
	(* What to do with these separated files?
		Analyze each and do:
		List of {{Eperps},{quantities}}
		Quantities would be EVs and Etrs and also f0, alpha, abs
		*)
	Table[
		{
			thismin=Import[ToString@resultfiles[[i]]],
			MUTABMIN=Append[MUTABMIN, {thismin[[1]],thismin[[2]][[3]]}],
			EVTABMIN=Append[EVTAB, {thismin[[1]],thismin[[2]][[4]]}],
			ETRTABMIN=Append[ETRTABMIN, {thismin[[1]],Table[thismin[[2]][[4]][[j]][[j]]-thismin[[2]][[4]][[j+1]][[j+1]],{j,Length[thismin]-1}]}]
		},
		{i, Length[resultfiles]}
	];
	Export["fin-mutab-min.m",MUTABMIN//Transpose];
	Export["fin-evtab-min.m",EVTABMIN//Transpose];
	Export["fin-etrtab-min.m",EVTABMIN//Transpose];
]
