#!/bin/bash
# Goal of this script/function is to: 
#	Build up a Matematica script in parts
#	Run a mathematica script
#	be able to change inputs of mathematica script (or function)
#	also be able to dictate the name of the mathematica job AND the name of the output
# NOTE: can use mathematica Save[] to append to files... So I can log steps taken and parameters used at each step by using, e.g.
# Save['combinedlog.txt', "finished calculation of X at Y distance and Z other parameters"]
# Save[ToString@StringForm[log for this event,...], "same as above"]

# Set some constants and read inputs
export jobname=$1
export todaysdate=$(date -I)_
export logsuffix="_output.log"
export projdir=/home/mbrunetti/tmdc/results/$todaysdate$jobname

mkdir $projdir
export filename=$todaysdate$jobname

cat <<EOF > tmdc-make.m
(thedata=exactCOULsuite[EXmatrange,4.9,{1,10,1},3,na,damp,"$filename","pdf"]);//AbsoluteTiming
Export["$projdir/$filename.m",thedata]//AbsoluteTiming
EOF

cat <<EOF > tmdc-analysis.m
types={"MoS2-hi","MoS2-lo","MoSe2-hi","MoSe2-lo","WS2-hi","WS2-lo","WSe2-hi","WSe2-lo"};
Export["$projdir/$filename-eigenenergies.m",Table[{types[[t]],Table[{nhbn-2,thedata[[1]][[t]][[nhbn]][[2]]},{nhbn,3,11}]},{t,8}]]
(* Extract r^2 *)
Export["$projdir/$filename-radius.m",Table[{types[[t]],Table[{nhbn-2,extractradius[thedata,EXmatrange,t,nhbn,1,1]},{nhbn,3,11}]},{t,8}]]
(* Extract f0, 1\[Rule]2 *)
Export["$projdir/$filename-f012.m",Table[{types[[t]],Table[{nhbn-2,extracthiloF0[thedata,EXmatrange,t,nhbn,2,1]},{nhbn,3,11}]},{t,8}]]
(* Extract f0, 1\[Rule]3 *)
Export["$projdir/$filename-f013.m",Table[{types[[t]],Table[{nhbn-2,extracthiloF0[thedata,EXmatrange,t,nhbn,3,1]},{nhbn,3,11}]},{t,8}]]
(* Extract \[Alpha] 1\[Rule]2 *)
Export["$projdir/$filename-alpha12.m",Table[{types[[t]],Table[{nhbn-2,extractalpha[thedata,EXmatrange,t,nhbn,2,1]},{nhbn,3,11}]},{t,8}]]
(* Extract \[Alpha] 1\[Rule]3 *)
Export["$projdir/$filename-alpha13.m",Table[{types[[t]],Table[{nhbn-2,extractalpha[thedata,EXmatrange,t,nhbn,3,1]},{nhbn,3,11}]},{t,8}]]
(* Extract \[ScriptCapitalA] 1\[Rule]2 *)
Export["$projdir/$filename-afac12.m",Table[{types[[t]],Table[{nhbn-2,extractAFAC[thedata,EXmatrange,t,nhbn,2,1]},{nhbn,3,11}]},{t,8}]]
(* Extract \[ScriptCapitalA] 1\[Rule]3 *)
Export["$projdir/$filename-afac13.m",Table[{types[[t]],Table[{nhbn-2,extractAFAC[thedata,EXmatrange,t,nhbn,3,1]},{nhbn,3,11}]},{t,8}]]
EOF

cat tmdc-init.m tmdc-make.m tmdc-analysis.m > test.m

cp /home/mbrunetti/tmdc/test.m $projdir/test.m

# Build the job submission script (read inputs and modify file names, etc.)
cat <<EOF > submit_mathematica.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here: 
#PBS -N $jobname

#DO NOT CHANGE THE NODE NUMBER:
#PBS -l nodes=node27:ppn=1

#Combine output and error files:
#PBS -j oe

#If you want specific log filename, use the following line:
#PBS -o $todaysdate$jobname$logsuffix 

echo $PBS_O_WORKDIR

echo "Starting Mathematica job"

cd $PBS_O_WORKDIR

#Submit mathematica job
#NOTE: Your test.m file SHOULD include the command "Quit[]" as the last command in the file

math -script tmdc/test.m

echo "Job finished"
EOF

qsub submit_mathematica.pbs
