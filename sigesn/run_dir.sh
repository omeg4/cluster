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

echo "Parameters:"
echo "Band gap? [meV]"
read ingap
echo "Buckling? [nm]"
read buck
echo "v_F? [m/s]"
read vF
echo "Monolayer Thickness? [nm]"
read thicc
echo "Epsilon of the material?"
read eps
echo "Kappa for the environment?"
read kappa

echo "E_z initial? [V/Ang]"
read ezi
echo "E_z final? [V/Ang]"
read ezf
echo "E_z step? [V/Ang]"
read ezstep

echo "Initializing direct calculations"
cat <<EOF > callfuncs.m
SetDirectory["$(pwd)"]
params={$ingap,$buck,$vF,$thicc,$eps};
etab={$ezi,$ezf,$ezstep};
Export["inp.m",{params,$kappa,etab,{0,0,1}}]
Export["diag1.txt","Params and Etab initialized"]
Export["suite.m",SiGeSuite[3,params,etab,$kappa]];
labels={"min","max"};
Export["suitediag.txt","Suite run complete, beginning analysis of data files."]
Export["proc.m",ProcessSuite[]];
Quit[]
EOF
cat /home/mbrunetti/cluster/sigesn/params.m /home/mbrunetti/cluster/sigesn/f-dirkeld.m /home/mbrunetti/cluster/sigesn/f-Dsuite.m /home/mbrunetti/cluster/sigesn/f-normEF.m /home/mbrunetti/cluster/sigesn/f-filemine.m callfuncs.m > test.m
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

