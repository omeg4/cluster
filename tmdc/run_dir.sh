#!/bin/bash
#

cd results/

export startdate=$(date -I)

export DorI="DIR"

echo "Material?"

read mater

echo "What's the name of this job?"
read jobby

export projdir="$jobby-$DorI-$mater-$startdate"

mkdir $projdir
cd $projdir

export fullpath=$(pwd)

echo "m_e?"
read masse

echo "m_h?"
read massh

echo "chi_2D?"
read chi2d

echo "Kappa?"
read kappa

echo "Initializing direct calculations"
cat <<EOF > callfuncs.m
SetDirectory["$(pwd)"]
params={$masse,$massh,$chi2d};
Export["inp.m",{params,$kappa}]
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",DirKeld[3,params,kappa]];
Quit[]
EOF
cat /home/mbrunetti/cluster/tmdc/params.m /home/mbrunetti/cluster/tmdc/f-dirkeld.m /home/mbrunetti/cluster/tmdc/f-Dsuite.m /home/mbrunetti/cluster/tmdc/f-filemine.m callfuncs.m > test.m
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

