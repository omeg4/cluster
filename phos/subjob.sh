#!/bin/bash
#

cd results/

export startdate=$(date -I)

export jobby=$1
export kappa=$2
export di=$3
export df=$4
export pot=$5
export eps=$6
export nmax=10

export projdir="$jobby-$startdate"

mkdir $projdir
cd $projdir

echo $(pwd)

export fullpath=$(pwd)

cat <<EOF > funccall.m
SetDirectory["$(pwd)"]
kappa = $kappa;
Export["inps.txt",{$jobby,mus,chiphos,kappa,{$di,$df},$pot,$eps}];
result=Table[
	CompPhos[$nmax,{mus[[i]][[1]],mus[[i]][[2]],chiphos,kappa},nhbn,$pot,$eps],
	{i,4},{nhbn,$di,$df}
	];
Export["suite.m",result];
Export["proc.m",ProcessPhosInd[result]];
Quit[]
EOF

cat ../../mmaconsts.m ../../phosconsts.m ../../phosfuncs.m funccall.m > test.m

cat <<EOF > submit_mathematica.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here:
#PBS -N $jobby

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
