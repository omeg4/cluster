#!/bin/bash
#

cd results/

export startdate=$(date -I)

export projdir="BFTORP"

mkdir $projdir
cd $projdir

cp ../../f-filemine.m ../../all-direct-defs.m .

export fullpath=$(pwd)

cat <<EOF > test.m
SetDirectory["$(pwd)"]
Import["f-filemine.m"]
Import["all-direct-defs.m"]
Export["BFTORP.m",
  ParallelTable[
   RegionPlot[
    {
     t < QuantityMagnitude[
        mc["TcKavVs"][[type]][[nn]][[ep]][QConc[conc*10^14], de], 
        "Kelvins"] &&
      t < QuantityMagnitude[
        Re[mc["ElpupVs"][[type]][[nn]][[ep]][QFreq[10^13], de, 
            0][[3]]]/(2*KB), "Kelvins"],
     QuantityMagnitude[
       mc["TcKavVs"][[type]][[nn]][[ep]][QConc[conc*10^14], de], 
       "Kelvins"] < t < 
      QuantityMagnitude[
       Re[mc["ElpupVs"][[type]][[nn]][[ep]][QFreq[10^13], de, 
           0][[3]]]/(2*KB), "Kelvins"],
     t > QuantityMagnitude[
       Re[mc["ElpupVs"][[type]][[nn]][[ep]][QFreq[10^13], de, 
           0][[3]]]/(2*KB), "Kelvins"]
     }, {conc, 1, 50}, {t, 150, 350},
    PlotLegends -> {"T < (\!\(\*SubscriptBox[\(T\), \(pol\)]\) & \!\(\
\*SubscriptBox[\(T\), \(c\)]\)) - Stable BEC", 
      "(\!\(\*SubscriptBox[\(T\), \(c\)]\) < T < \
\!\(\*SubscriptBox[\(T\), \(pol\)]\)) - Stable normal LP", 
      "(T > \!\(\*SubscriptBox[\(T\), \(pol\)]\)) - No stable LP", ""},
    FrameLabel -> {"\!\(\*SubscriptBox[\(n\), \(pol\)]\) \
[\[Times]\!\(\*SuperscriptBox[\(10\), \(14\)]\) \!\(\*SuperscriptBox[\
\(m\), \(-2\)]\)]", "T [K]"},
    PerformanceGoal -> "Speed",
    ImageSize -> {1000, 1000},
    LabelStyle -> Directive[Black, 30]
    ],
   {mc, {simc, gemc, snmc, ssmc, bsmc}}, {type, 2}, {nn, 10}, {ep, 
    bsmc["Emax"]}, {de, -0.1, 0.1, 0.01}
   ]
  ] // AbsoluteTiming
Quit[]
EOF

cat <<EOF > submit_mathematica.pbs
#!/bin/sh

#Important: do not remove "#" symbol before PBS, keep it like that: "#PBS"


#You can set your job name here:
#PBS -N "BFTORP"

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

