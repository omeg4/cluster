#!/bin/bash
# Goal of this script/function is to: 
#	Build up a Matematica script in parts
#	Run a mathematica script
#	be able to change inputs of mathematica script (or function)
#	also be able to dictate the name of the mathematica job AND the name of the output
# NOTE: can use mathematica Save[] to append to files... So I can log steps taken and parameters used at each step by using, e.g.
# Save['combinedlog.txt', "finished calculation of X at Y distance and Z other parameters"]
# Save[ToString@StringForm[log for this event,...], "same as above"]

# Set some constants and read inputs
export s=$1
export ns=$2
export jobname=$3
export todaysdate=$(date -I)_
export logsuffix="_output.log"
export projdir=/home/mbrunetti/phos/results/$todaysdate$jobname/

mkdir $projdir

# Build Mathematica script in parts: First the general constants, then phos constants, then phos functions (for NDEigensystem, etc.) and then finally run the function with the inputs from the script
cat <<EOF > funccall.m
s = $s
ns = $ns
projdir = ToString[$todaysdate$jobname]
filename = ToString[$jobname]
result=Table[
	{
		mux,
		muy,
		Table[{nhbn,CompPhosKeld[mux,muy,4.89,rho,nhbn*lBN,10,s,ns]},{nhbn,10}]
	},
	{mux,{0.062963,0.0910365}},{muy,{0.659901,0.967742}}
]//AbsoluteTiming;
Export[
	"$projdir$todaysdate$jobname.m",
	result
];
EOF

cat <<EOF > analyzedata.m
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
Export[ToString@StringForm["/home/mbrunetti/phos/results/\`1\`/processed.m",projdir],phosprocess]
Quit[]
EOF

cat mmaconsts.m phosconsts.m phosfuncs.m funccall.m analyzedata.m > test.m

cp /home/mbrunetti/phos/test.m $projdir/test.m

# Build the job submission script (read inputs and modify file names, etc.)
cat <<EOF > submit_mathematica.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here: 
#PBS -N $jobname

#DO NOT CHANGE THE NODE NUMBER:
#PBS -l nodes=node27:ppn=1

#Combine output and error files:
#PBS -j oe

#If you want specific log filename, use the following line:
#PBS -o $todaysdate$jobname$logsuffix 

echo $PBS_O_WORKDIR

echo "Starting Mathematica job"

cd $PBS_O_WORKDIR

#Submit mathematica job
#NOTE: Your test.m file SHOULD include the command "Quit[]" as the last command in the file

math -script phos/test.m

echo "Job finished"
EOF

qsub submit_mathematica.pbs
