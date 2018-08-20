#!/bin/bash
# Goal of this script/function is to: 
#	Build up a Matematica script in parts
#	Run a mathematica script
#	be able to change inputs of mathematica script (or function)
#	also be able to dictate the name of the mathematica job AND the name of the output
#########
# Better organization:
# mkdir for each job. In each folder, put:
# 	full script that started the job
#		log
#		raw data
#		analysis results

# Set some constants and read inputs
export jobname=$1
export todaysdate=$(date -I)_
export logsuffix="_output.log"
export outpath="/home/mbrunetti/phos/results/"$todaysdate$jobname

# Build Mathematica script in parts: First the general constants, then phos constants, then phos functions (for NDEigensystem, etc.) and then finally run the function with the inputs from the script
cat <<EOF > analysisscript.m
data = Import["/home/mbrunetti/phos/files/2-12_sns1k.m"];
Module[
	{
		EVs = data[[2]][[1]],
		EFs = data[[2]][[2]],
		mux = mus[[1]][[1]],
		muy = mus[[1]][[2]]
	},
	ornorm = Table[
			NIntegrate[EFs[[i]]*EFs[[j]],{x,-1000,1000},{y,-1000,1000},MaxRecursion->100,WorkingPrecision->100],
		{i,10},{j,10}
	];
	f0 = Table[
		{
			(EVs[[i]]-EVs[[1]])*H2eV,
			2*mux*(EVs[[i]]-EVs[[1]])*(NIntegrate[EFs[[i]]*x*EFs[[1]],{x,-1000,1000},{y,-1000,1000},MaxRecursion->100,WorkingPrecision->100])^2,
			2*muy*(EVs[[i]]-EVs[[1]])*(NIntegrate[EFs[[i]]*y*EFs[[1]],{x,-1000,1000},{y,-1000,1000},MaxRecursion->100,WorkingPrecision->100])^2
		},
		{i,10}
	];
	plotEn = ListPlot[
		Table[
			-1000*H2eV*EVs[[i]],
			{i,10}
		],
		PlotTheme->"Detailed"
	];
	plotEtr = ListPlot[
		Table[
			1000*H2eV*(EVs[[i]]-EVs[[1]]),
			{i,2,10}
		],
		PlotTheme->"Detailed"
	];
	EFplots = Table[
		{
			EVs[[i]],
			Plot[EFs[[i]]/.y->0,{x,-1000,1000},AxesLabel->{"x","\[psi]"}],
			Plot[EFs[[i]]/.x->0,{y,-1000,1000},AxesLabel->{"y","\[psi]"}],
		},
		{i,10}
	];
	Export["/home/mbrunetti/phos/files/2-12_sns1k_analysis.m",{ornorm,f0}];
	Export["/home/mbrunetti/phos/files/2-12_sns1k_Ebplot.pdf",plotEn];
	Export["/home/mbrunetti/phos/files/2-12_sns1k_Etrplot.pdf",plotEtr];
]
EOF

cat mmaconsts.m phosconsts.m phosfuncs.m analysisscript.m > test.m

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
#PBS -o $jobname$logsuffix 

echo $PBS_O_WORKDIR

echo "Starting Mathematica job"

cd $PBS_O_WORKDIR

#Submit mathematica job
#NOTE: Your test.m file SHOULD include the command "Quit[]" as the last command in the file

math -script phos/test.m

echo "Job finished"
EOF

mkdir /home/mbrunetti/phos/results/$todaysdate$jobname
cp /home/mbrunetti/phos/test.m /home/mbrunetti/phos/results/$todaysdate$jobname/test.m
cd /home/mbrunetti/phos/results/$todaysdate$jobname

qsub submit_mathematica.pbs
