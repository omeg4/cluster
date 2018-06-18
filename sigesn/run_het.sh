#!/bin/bash
#

cd results/

export startdate=$(date -I)

export DorI="IND"

echo "HetCoul or HetKeld ?"
read funcname

echo "What's the name of this job?"
read jobby

export projdir="$jobby-$DorI-$startdate"

mkdir $projdir
cd $projdir

export fullpath=$(pwd)

echo "Initializing indirect calculations using "$funcname
cat <<EOF > callfuncs.m
SetDirectory["$(pwd)"]
{params1,params2,kappa,etab,dtab}=Import["../../inputs.m"]
Export["inp.m",{params1,param2,kappa,etab,dtab}];
"Inputs saved. Initializing suite.">>>"diag.txt"
{time,suite}=IndHeteroSuite[params1,params2,kappa,etab,dtab,$funcname];
Export["suite.m",suite];
"Suite run complete. Processing data.">>>"diag.txt"
ToString@StringForm["Time to run suite: `` hours",time/3600]>>>"diag.txt"
labels={"min","max"};
Export["proc.m",ProcessSuite[]];
"Processing complete. Checking normalization.">>>"diag.txt"
checksuitenorm[Import["suite.m"]];
Quit[]
EOF

cat /home/mbrunetti/cluster/sigesn/params.m /home/mbrunetti/cluster/sigesn/f-hetkeld.m /home/mbrunetti/cluster/sigesn/f-hetcoul.m /home/mbrunetti/cluster/sigesn/f-Isuite.m /home/mbrunetti/cluster/sigesn/f-filemine.m callfuncs.m > test.m
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


