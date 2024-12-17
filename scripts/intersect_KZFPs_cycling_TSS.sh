#! /bin/bash

bedtools intersect -wo -a out/tss_rhythmic_genes_sorted_acrophase.bed -b data/hg19_peaks_un_filt_macs80.bed > out/tss_rhythmic_genes_vs_kzfp_peaks.bed

bedtools intersect -wo -a data/tss_clustered.bed -b data/hg19_peaks_un_filt_macs80.bed > out/promoters_hg19_genes_vs_kzfp_peaks.bed


# KZFP ZNFs vs binding by KZFPs as measured by ChIP in the Trono lab: 
bedtools intersect -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b data/hg19_peaks_un_filt_macs80.bed > out/kzfp_znfs_vs_kzfp_peaks.bed
