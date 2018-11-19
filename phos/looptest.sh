#!/bin/bash
#

while IFS='' read -r line || [[ -n "$line" ]]; do
	cd ~/cluster/phos
	IFS="|" read -r -a element <<< $line
	export jobby=${element[0]}
	export kappa=${element[1]}
	export di=${element[2]}
	export df=${element[3]}
	export pot=${element[4]}
	export eps=${element[5]}
	./subjob.sh $jobby $kappa $di $df $pot $eps
done < "$1"
