#!/bin/bash
#

cd results/

export startdate=$(date -I)

export jobby=$1
export cork=$2
export dorr=${10}
export mux=$3
export muy=$4
export chi2d=$5
export kappa=$6
export dee=$7
export sss=$8
export nss=$9
export nmax=4

export projdir="$jobby-$sss-$startdate"

mkdir $projdir
cd $projdir

echo $(pwd)

export fullpath=$(pwd)

cat <<EOF > funccall.m
SetDirectory["$(pwd)"]
Export["inps.txt",{$jobby,$cork,$dorr,$mux,$muy,$chi2d,$kappa,$dee,$sss,$nss}];
result=CompPhos[$nmax,{$mux,$muy,$chi2d,$kappa},$dee,$cork,$dorr,$sss,$nss];
Export["suite.m",result]
Quit[]
EOF

cat ../../mmaconsts.m ../../phosconsts.m ../../phosfuncs.m funccall.m > test.m

cat <<EOF > submit_mathematica.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here:
#PBS -N $jobby-$sss

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
