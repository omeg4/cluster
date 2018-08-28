#!/bin/bash
#

while IFS='' read -r line || [[ -n "$line" ]]; do
	cd ~/cluster/sigesn
	IFS="|" read -r -a element <<< $line
	export mater=${element[0]}
	export jobby=${element[1]}
	export ingap=${element[2]}
	export buckl=${element[3]}
	export vferm=${element[4]}
	export thicc=${element[5]}
	export epsil=${element[6]}
	export kappa=${element[7]}
	export ezini=${element[8]}
	export ezfin=${element[9]}
	export ezstp=${element[10]}
	./run_dir.sh $mater $jobby $ingap $buckl $vferm $thicc $epsil $kappa $ezini $ezfin $ezstp
done < "$1"
