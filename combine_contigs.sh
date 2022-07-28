#!/bin/bash

#This takes in a directory with a lot of .fna FASTA files that have a bunch of separate entries (like if there are >contig1...>contig2...etc) and makes it into one big fasta file with no separate >'s and just one header that is ">filename"

for f in ./*.fna; do
	sed -i '/^>/d' $f 
	sed -i '1i >'$f'' $f
	sed -e 's/.fna//g' -i $f
	sed -e 's_./__g' -i $f
done
