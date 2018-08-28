#!/bin/bash
#

while IFS='' read -r line || [[ -n "$line" ]]; do
	cd ~/cluster/phos
	IFS="|" read -r -a element <<< $line
	export jobby=${element[0]}
	export cork=${element[1]}
	export mux=${element[2]}
	export muy=${element[3]}
	export chi2d=${element[4]}
	export kappa=${element[5]}
	export dee=${element[6]}
	export sss=${element[7]}
	export nss=${element[8]}
	./subjob.sh $jobby $cork $mux $muy $chi2d $kappa $dee $sss $nss
done < "$1"
