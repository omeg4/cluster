#!/bin/bash
#

cd results/

export startdate=$(date -I)

echo "Direct or indirect? (DIR or IND *only*)"

read DorI

echo "Material?"

read mater

echo "What's the name of this job?"
read jobby

export projdir="$DorI-$mater-$startdate-$jobby"

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
#################################################
#################################################
if [ $DorI = "DIR" ]
then
	#Direct
	echo "Initializing direct calculations"
	cat <<EOF > callfuncs.m
SetDirectory["$(pwd)"]
fpath=ToString[StringForm["/home/mbrunetti/cluster/sigesn/results/\`1\`/",$projdir]];
params={$ingap,$buck,$vF,$thicc,$eps};
etab={$ezi,$ezf,$ezstep};
Export["diag1.txt","Params and Etab initialized"]
suite=SiGeSuite[3,params,etab,$kappa];
Export["suitediag.txt","Suite run complete"]
Export["results.m",suite];
Quit[]
EOF
	cat /home/mbrunetti/cluster/sigesn/sgsinit.m callfuncs.m > test.m
	pwd
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
#################################################
#################################################
elif [ $DorI = "IND" ]
then

	#Indirect
	echo "D initial? [N_hBN]"
	read di
	echo "D final? [N_hBN]"
	read df
	echo "D step? [N_hBN]"
	read ds
	echo "Initializing indirect calculations"
	cat <<EOF > callfuncs.m
params={$ingap,$buck,$vF,$thicc,$eps};
etab={$ezi,$ezf,$ezstep};
dtab={$di,$df,$dstep};
suite=SGSIndSuite[params,etab,dtab,ToString[$mater-$jobby],ToString[$fullpath]];
fpath=ToString[$fullpath];
Export[fpath<>"results.m",suite];
Quit[]
EOF
#################################################
#################################################
elif [ $DorI = "PREV" ]
then
prevfile=$(cat prev.txt)
$DorI
$mater
$jobby
$ingap
$buck
$vF
$thicc
$eps
$kappa
$ezi
$ezf
$ezstep
else

	echo "Could not understand if you wanted DIR or IND calculation. Please try again"
	rm -rf "$mater$jobby"
fi

cat <<EOF > prev.txt
$DorI
$mater
$jobby
$ingap
$buck
$vF
$thicc
$eps
$kappa
$ezi
$ezf
$ezstep
EOF
