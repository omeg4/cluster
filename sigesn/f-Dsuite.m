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
	max,
    },
	ne=0;
    Table[{ne = ne + 1,
	Export[ToString@StringForm["calcs_e`1`.m",IntegerString[ne,10,3]],
	{
		{"min", DirKeld[mmax, matpars[[2 ;;]], Egap, Eperp, \[Kappa], 1]},
		{"max", DirKeld[mmax, matpars[[2 ;;]], Egap, Eperp, \[Kappa], -1]}
	}]},
	{Eperp, eztab[[1]], eztab[[2]], eztab[[3]]}
	];
	Export["tabledone.txt","Table done, files saved. Starting analysis of data"];
]
