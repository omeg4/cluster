#!/bin/bash
#

cd results/

export startdate=$(date -I)

export projdir="rawRP"

mkdir $projdir
cd $projdir

cp ../../f-filemine.m ../../all-direct-defs.wl .

export fullpath=$(pwd)

cat <<EOF > test.m
SetDirectory["$(pwd)"]
Import["f-filemine.m"]
Import["all-direct-defs.wl"]
Export["rawRP.m",
ParallelTable[
	{
	  mc[[2]],
	  type,
	  nn,
	  FPeperp[mc[[1]]["prc"],ep],
	  QFreq[gamex*(10^13)],
	  de,
	  conc,
	  t,
	  {
		{
			mc[[1]]["exres"][[type]][[ep]],
			mc[[1]]["Eck"][[type]][[ep]][de,0],
			mc[[1]]["Eck"][[type]][[ep]][de,0]-mc[[1]]["exres"][[type]][[ep]]
		},
		{
			QFreq[gamex*(10^13)],
			mc[[1]]["gamcav"][[type]][[nn]][[ep]][de],
		},
		{

			mc[[1]]["Vsav"][[type]][[nn]][[ep]][de],
			Re[mc[[1]]["Elpup"][[type]][[nn]][[ep]][de,0][[1]]],
			Re[mc[[1]]["Elpup"][[type]][[nn]][[ep]][de,0][[3]]/2],
			Min[mc[[1]]["exres"][[type]][[ep]],mc[[1]]["Eck"][[type]][[ep]][de,0]]-Re[mc[[1]]["Elpup"][[type]][[nn]][[ep]][de,0][[1]]]
		},
		{
			mc[[1]]["TcKavVs"][[type]][[nn]][[ep]][QConc[conc*10^14], de],
			mc[[1]]["CritnKavVs"][[type]][[nn]][[ep]][t,de]
		},
		QuantityMagnitude[mc[[1]]["Vsav"][[type]][[nn]][[ep]][de]] > QuantityMagnitude[Abs[HB*(mc[[1]]["gamcav"][[type]][[nn]][[ep]][de] - QFreq[gamex*(10^13)])]/2,"Millielectronvolts"],
		t < QuantityMagnitude[mc[[1]]["TcKavVs"][[type]][[nn]][[ep]][QConc[conc*10^14], de], "Kelvins"] && t < QuantityMagnitude[Re[mc[[1]]["ElpupVs"][[type]][[nn]][[ep]][QFreq[gamex*(10^13)], de, 0][[3]]]/(2*KB), "Kelvins"],
		QuantityMagnitude[mc[[1]]["TcKavVs"][[type]][[nn]][[ep]][QConc[conc*10^14], de],"Kelvins"] < t < QuantityMagnitude[Re[mc[[1]]["ElpupVs"][[type]][[nn]][[ep]][QFreq[gamex*(10^13)], de, 0][[3]]]/(2*KB), "Kelvins"],
		t > QuantityMagnitude[Re[mc[[1]]["ElpupVs"][[type]][[nn]][[ep]][QFreq[gamex*(10^13)], de, 0][[3]]]/(2*KB), "Kelvins"]
		}
	},
	{mc, {{simc, 1}, {gemc, 2}, {snmc, 3}, {ssmc, 4}, {bsmc, 5}}},
	{type, 2},
	{nn, 10},
	{ep,mc[[1]]["NNEtr"][[type]],mc[[1]]["Emax"]},
	{gamex, 5},
	{de, -0.1, 0.1, 0.01},
	{conc, 1, 100},
	{t, 150, 350}
	]
]//AbsoluteTiming
Quit[]
EOF

cat <<EOF > submit_mathematica.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here:
#PBS -N "rawRP"

#DO NOT CHANGE THE NODE NUMBER:
#PBS -l nodes=node27:ppn=1

#Combine output and error files:
#PBS -j oe

#If you want specific log filename, use the following line:
#PBS -o pbs.log

export PBS_O_WORKDIR=$(pwd)
echo \$PBS_O_WORKDIR

echo "Starting Mathematica job"

cd \$PBS_O_WORKDIR

#Submit mathematica job
#NOTE: Your test.m file SHOULD include the command "Quit[]" as the last command in the file

math -script \$PBS_O_WORKDIR/test.m

echo "Job finished"

EOF
qsub submit_mathematica.pbs

