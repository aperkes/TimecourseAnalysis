#! /bin/bash

infile=$1
outfile="$1.info"

echo "Gene info" > $outfile

while read -r g; do
    echo $g
    #zgrep -h $g /data/sequencing/MollyGenome/GCF_000485575.1_Poecilia_formosa-5.1.2_rna.fna.gz >> $outfile
    grep -h $g ./gene_names.txt >> $outfile
    done <$infile

