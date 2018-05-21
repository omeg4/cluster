#!/bin/bash
#

cat <<EOF > callanalysis.m
SetDirectory["$(pwd)"]
ProcessSuite[];
Quit[]
EOF
cat /home/mbrunetti/cluster/sigesn/params.m /home/mbrunetti/cluster/sigesn/f-dirkeld.m /home/mbrunetti/cluster/sigesn/f-Dsuite.m /home/mbrunetti/cluster/sigesn/f-normEF.m /home/mbrunetti/cluster/sigesn/f-filemine.m callanalysis.m > test2.m
cat <<EOF > submit_mathematica_analysis.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here:
#PBS -N $mater-$jobby-analysis

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

math -script \$PBS_O_WORKDIR/test2.m

echo "Job finished"

EOF
qsub submit_mathematica_analysis.pbs


