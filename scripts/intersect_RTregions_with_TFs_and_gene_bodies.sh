#! /bin/bash

# liftover of RTregions from hg38 to hg19
mkdir out/rif1ko_hct_RTregions/hg19
# TODO simplify as one file only, we don't need to have them separated
for BEDPATH in out/rif1ko_hct_RTregions/hg38/*; do
    BED=$(basename ${BEDPATH})
    bin_external/liftOver -bedPlus=4 $BEDPATH data/hg38ToHg19.over.chain.gz out/rif1ko_hct_RTregions/hg19/$BED out/temp/unlifted.bed;
done
bin_external/liftOver -bedPlus=4 out/rif1ko_hct_RTregions/rif1ko_hct_regions_hg38.bed data/hg38ToHg19.over.chain.gz out/rif1ko_hct_RTregions/rif1ko_hct_regions_hg19.bed out/temp/unlifted.bed

# ENCODE data x rif1ko HCT regions
mkdir out/encode_x_rif1ko_hct_RTregions/

for CHIPBEDPATH in data/ENCODE_k562_tf_peaks/bed_filtered_hg19/*; do
    CHIPBED=$(basename ${CHIPBEDPATH})
    CHIP="${CHIPBED%.*}"
    
    #bedtools intersect -wo -a out/rif1ko_hct_RTregions/rif1ko_hct_regions_hg19.bed -b $CHIPBEDPATH > out/encode_x_rif1ko_hct_RTregions/$CHIPBED;
    bedtools intersect -wo -a out/rif1ko_hct_RTbins.bed -b $CHIPBEDPATH > out/encode_x_rif1ko_hct_RTregions/$CHIPBED;
done


# same for najafabadi and schmitges chip-seq data, all ZNFs
mkdir out/schmitges_x_rif1ko_hct_RTregions/
for BEDPATH in data/schmitges_all_znfs_chips/hg19_peaks_chipatlas/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a out/rif1ko_hct_RTbins.bed -b $BEDPATH > out/schmitges_x_rif1ko_hct_RTregions/$BED;
done


# intersecting with promoters
mkdir out/najafabadi_x_rif1ko_hct_RTregions/
for BEDPATH in data/najafabadi_all_znfs_chips/hg19_peaks_80/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a out/rif1ko_hct_RTbins.bed -b $BEDPATH > out/najafabadi_x_rif1ko_hct_RTregions/$BED;
done


# intersecting KZFP Chip-seq with the RT regions
#bedtools intersect -wo -a out/rif1ko_hct_RTregions/rif1ko_hct_regions_hg19.bed -b data/hg19_peaks_un_filt_macs80.bed > out/rif1ko_hct_RTregions_vs_kzfp_peaks.bed
bedtools intersect -wo -a out/rif1ko_hct_RTbins.bed -b data/hg19_peaks_un_filt_macs80.bed > out/rif1ko_hct_RTbins_vs_kzfp_peaks.bed


# intersecting gene bodies with the RT regions
#bedtools intersect -wo -a out/rif1ko_hct_RTregions/rif1ko_hct_regions_hg19.bed -b data/1807_hg19_ens_coding_genes_body_symbol.bed > out/rif1ko_hct_RTregions_vs_gene_bodies.bed
bedtools intersect -wo -a out/rif1ko_hct_RTbins.bed -b data/1807_hg19_ens_coding_genes_body_symbol.bed > out/rif1ko_hct_RTbins_vs_gene_bodies.bed