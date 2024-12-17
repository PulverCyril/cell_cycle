#! /bin/bash

# downloading the RIF1 binding files
mkdir out/temp/rif1

cd out/temp/rif1
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160563/suppl/GSE160563_HCT116_RIF1_CutnRun_log2_FE_over_Input.bigwig

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160563/suppl/GSE160563_HCT116_RIF1_CutnRun_rep2_log2_FE_over_Input.bigwig


cd ../../../

# liftOver the hg19 gene bodies coordinates to hg38, sorted by acrophase
bin_external/liftOver -bedPlus=6 out/gene_bodies_all_genes_sorted_acrophase_tssstart.bed data/hg19ToHg38.over.chain.gz out/gene_bodies_all_genes_sorted_acrophase_tssstart_hg38.bed out/temp/unlifted.bed

# scoring over the bigwig
computeMatrix scale-regions \
 -S out/temp/rif1/GSE160563_HCT116_RIF1_CutnRun_rep2_log2_FE_over_Input.bigwig \
 -R out/gene_bodies_all_genes_sorted_acrophase_tssstart_hg38.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/gene_bodies_all_genes_sorted_acrophase_tssstart_hg38_rif1.tab.gz \
 --outFileNameMatrix out/temp/gene_bodies_all_genes_sorted_acrophase_tssstart_hg38_rif1.tab > out/temp/gene_bodies_all_genes_sorted_acrophase_tssstart_hg38_rif1.log

# sorted as in data/1807_hg19_ens_coding_genes_body_symbol.bed
bin_external/liftOver -bedPlus=6 data/1807_hg19_ens_coding_genes_body_symbol.bed data/hg19ToHg38.over.chain.gz out/1807_hg19_ens_coding_genes_body_symbol_hg38.bed out/temp/unlifted.bed

# scoring over the bigwig
computeMatrix scale-regions \
 -S out/temp/rif1/GSE160563_HCT116_RIF1_CutnRun_rep2_log2_FE_over_Input.bigwig \
 -R out/1807_hg19_ens_coding_genes_body_symbol_hg38.bed \
 -b 1000 -a 0 \
 --sortRegions no\
 -out out/temp/1807_hg19_ens_coding_genes_body_symbol_hg38_rif1_cutnrun_rep2_log2_FE_over_Input.tab.gz \
 --outFileNameMatrix out/temp/1807_hg19_ens_coding_genes_body_symbol_hg38_rif1_cutnrun_rep2_log2_FE_over_Input.tab \
 --outFileSortedRegions out/temp/1807_hg19_ens_coding_genes_body_symbol_hg38_rif1_cutnrun_rep2_log2_FE_over_Input_ID.tsv
