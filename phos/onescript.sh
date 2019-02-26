#!/bin/bash
#
# Since this script is in the Phos folder we don't need to worry about TMDC and Xene. Just do this for all the parameters and then we can adapt for other materials and place the scripts in the corresponding folder.

### Part 1: read series of inputs for parameters and run appropriate jobs
# What we want to do here is eliminate the need for a separate "joblist.txt" file and just keep everything in one place

# while IFS='' read -r line || [[ -n "$line" ]]; do
#     cd ~/cluster/phos/results
#     IFS="|" read -r -a element <<< $line
#     export startdate=$(date -I)
#     export jobby=${element[0]}
#     export kappa=${element[1]}
#     export di=${element[2]}
#     export df=${element[3]}
#     export pot=${element[4]}
#     export eps=${element[5]}
#     export nmax=10
cd ~/cluster/phos/results/

export jobby=$1
export startdate=$(date -I)

export projdir="$jobby-$startdate"

echo $projdir
mkdir $projdir
cd $projdir

echo $(pwd)

export fullpath=$(pwd)

cat <<EOF > funccall.m
SetDirectory["$(pwd)"];
time = AbsoluteTiming[ result = makeproc[]; ] [[1]];
Export["assoc.m", result];
Export["totaltime.txt", time];
Quit[]
EOF

cat ~/cluster/phos/neatphos.m funccall.m > test.m

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

# done < "$1"
