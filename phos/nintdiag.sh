#!/bin/bash
#
# Diagnostic submit

export DIR=~/cluster/phos/results/0408-nint
mkdir $DIR

cat <<EOF > end.m
Quit[]
EOF

cat neatphos.m nintdiag.m end.m > $DIR/test.m

cd $DIR

ln -s ../0323-short/0323-ndeout-2-eps0.9_maxi1000.m nde9.m
ln -s ../0323-short/0323-ndeout-2-eps0.10_maxi1000.m nde10.m

cat <<EOF > submit_mathematica.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here:
#PBS -N "nint"

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
