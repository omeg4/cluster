#!/bin/bash
# Goal of this script/function is to: 
#	Build up a Matematica script in parts
#	Run a mathematica script
#	be able to change inputs of mathematica script (or function)
#	also be able to dictate the name of the mathematica job AND the name of the output
# NOTE: can use mathematica Save[] to append to files... So I can log steps taken and parameters used at each step by using, e.g.
# Save['combinedlog.txt', "finished calculation of X at Y distance and Z other parameters"]
# Save[ToString@StringForm['log for this event,...], "same as above"]

# Set some constants and read inputs
export s=$1
export ns=$2
export jobname=$3
export todaysdate=$(date -I)_
export logsuffix="_output.log"

# Build Mathematica script in parts: First the general constants, then phos constants, then phos functions (for NDEigensystem, etc.) and then finally run the function with the inputs from the script
cat <<EOF > funccall.m
s = $s
ns = $ns
filename = $jobname
result=Table[
	{
		{mus[[j]][[1]],mus[[j]][[2]],(mus[[j]][[1]]*mus[[j]][[2]]/(mus[[j]][[1]]+mus[[j]][[2]]))},
		CompPhosKeld[mus[[j]][[1]],mus[[j]][[2]],4.89,rho,i*lBN,10,s,ns]
	},
	{i,10},{j,4}
];
Export[
	ToString@StringForm[
		"/home/mbrunetti/phos/files/\`1\`-\`2\`_\`3\`.m",
		DateObject[][[1]][[2]],
		DateObject[][[1]][[3]],
		filename
	],
	result
];
Quit[]
EOF

cat mmaconsts.m phosconsts.m phosfuncs.m funccall.m > test.m

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
