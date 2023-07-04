#! /bin/bash

echo Gene info > A_gene_info.txt

for g in 'cat A_genes.csv'; do
    zgrep -h $g /data/sequencing/MollyGenome/GCF_000485575.1_Poecilia_formosa-5.1.2_rna.fna.gz >> A_gene_info.txt
    done

