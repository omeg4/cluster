#!/bin/bash
#

cd results/

export startdate=$(date -I)

export DorI="DIR"

export mater=$1
export jobby=$2

export projdir="$jobby-$DorI-$mater-$startdate"

mkdir $projdir
cd $projdir

cp ../../f-dirkeld.m ../../run_dir.sh ../../f-Dsuite.m .

export fullpath=$(pwd)

export ingap=$3
export buck=$4
export vF=$5
export thicc=$6
export eps=$7
export kappa=$8
export ezi=$9
export ezf=${10}
export ezstep=${11}

echo "Initializing direct calculations"
cat <<EOF > callfuncs.m
SetDirectory["$(pwd)"]
params={$ingap,$buck,$vF,$thicc,$eps};
etab={$ezi,$ezf,$ezstep};
Export["inp.m",{params,$kappa,etab,{0,0,1}}]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",SiGeSuite[3,params,etab,$kappa]];
"Suite run complete. Processing data.">>>"diag.txt"
labels={"min","max"};
Export["proc.m",ProcessSuite[]];
"Processing complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
EOF
cat /home/mbrunetti/cluster/sigesn/params.m /home/mbrunetti/cluster/sigesn/f-dirkeld.m /home/mbrunetti/cluster/sigesn/f-Dsuite.m /home/mbrunetti/cluster/sigesn/f-filemine.m callfuncs.m > test.m
cat <<EOF > submit_mathematica.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here:
#PBS -N $mater-$jobby

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

