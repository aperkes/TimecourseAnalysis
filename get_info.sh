#! /bin/bash

infile="A_genes.csv"
outfile="A_gene_info.txt"

echo "Gene info" > $outfile

while read -r g; do
    zgrep -h $g /data/sequencing/MollyGenome/GCF_000485575.1_Poecilia_formosa-5.1.2_rna.fna.gz >> $outfile
    done <$infile

