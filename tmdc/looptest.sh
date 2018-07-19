#!/bin/bash
#

while IFS='' read -r line || [[ -n "$line" ]]; do
	cd ~/cluster/tmdc
	IFS="|" read -r -a element <<< $line
	export mater=${element[0]}
	export jobby=${element[1]}
	export bandg=${element[2]}
	export lambd=${element[3]}
	export alatt=${element[4]}
	export thopp=${element[5]}
	export chi2d=${element[6]}
	export kappa=${element[7]}
	./run_dir.sh $mater $jobby $bandg $lambd $alatt $thopp $chi2d $kappa
done < "$1"
