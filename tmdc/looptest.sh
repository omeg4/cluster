#!/bin/bash
#

while IFS='' read -r line || [[ -n "$line" ]]; do
	cd ~/cluster/tmdc
	IFS="|" read -r -a element <<< $line
	export mater=${element[0]}
	export jobby=${element[1]}
	export emass=${element[2]}
	export hmass=${element[3]}
	export chi2d=${element[4]}
	export kappa=${element[5]}
	./run_dir.sh $mater $jobby $bandg $lambd $alatt $thopp $chi2d $kappa
done < "$1"
