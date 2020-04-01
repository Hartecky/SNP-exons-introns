#!bin/bash

# Final file with coordinates of snps, exon start and exon end
coords='final_exons.txt'

# Listing all values
snp=( $( cat $coords | awk '{print $1}') )
exon_start=( $( cat $coords | awk '{print $2}') )
exon_end=( $( cat $coords | awk '{print $3}') )

# Working method but extremely slowly
# Iteration over every SNP value and checking 
# if each SNP is between exon_start and exon_end
# coordinates. Iteration through all file took
# too much time, so I figured out an awk method

i=0
counter=0
for value in ${exon_end[@]}; do
	new_val=$counter
	counter=0
	let "i++"
	for snps in ${snp[@]}; do
		if [[ $value > $snps ]]; then
			#statements
			let "counter++"
		else
			break
		fi
	done
	final=$(echo "scale=2; sqrt(($counter-$new_val)^2)" | bc)
	echo "Exon $i : $final SNPs"
done