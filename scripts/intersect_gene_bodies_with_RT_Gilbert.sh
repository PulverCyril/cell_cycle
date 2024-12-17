#! /bin/bash

mkdir out/RTregions_vs_gene_bodies/

for BEDPATH in out/bed/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a data/1807_hg19_ens_coding_genes_body_symbol.bed -b $BEDPATH > out/RTregions_vs_gene_bodies/$BED;
done