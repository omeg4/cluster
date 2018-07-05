#!/bin/bash
#

cd results/

export startdate=$(date -I)

export DorI="IND"

echo "What's the name of this job?"
read jobby

export projdir="$jobby-$DorI-$startdate"

mkdir $projdir
cd $projdir

export fullpath=$(pwd)

echo "Initializing indirect calculations"
cat <<EOF > callfuncs.m
SetDirectory["$(pwd)"]
{me,mh,chi2d,kappa,d}=Import["../../inputs.m"]
Export["inp.m",{me,mh,chi2d,kappa, d}];
"Inputs saved. Initializing suite.">>>"diag.txt"
Export["suite.m",IndKeld[3,{me,mh,chi2d},kappa,d]];
Quit[]
EOF

cat /home/mbrunetti/cluster/tmdc/params.m /home/mbrunetti/cluster/tmdc/f-indkeld.m /home/mbrunetti/cluster/tmdc/f-filemine.m callfuncs.m > test.m
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


