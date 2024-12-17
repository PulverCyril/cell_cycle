#! /bin/bash

bedtools intersect -loj -wo -a out/OutOfFrameZFA_KZFPs_ZFPs.bed -b data/martina_processed_files/2301_H3K9wt_DB_ADJP005_FC2.bed > out/OutOfFrameZFA_KZFPs_KZFPsdeltaK9_274KO_293T.bed

bedtools intersect -loj -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b data/martina_processed_files/2301_H3K9wt_DB_ADJP005_FC2.bed > out/kzfp_znfs_deltaK9_274KO_293T_sorted_cluster_acrophase.bed

# ZNF274 ENCODE K562: ENCFF742DPO
bedtools intersect -loj -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b data/ENCODE_k562_tf_peaks/bed_filtered_hg19/ENCFF742DPO.bed > out/kzfp_znfs_znf274_encode_ENCFF742DPO_sorted_cluster_acrophase.bed


# ZNF274 ENCODE K562: ENCFF045CFL
bedtools intersect -loj -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b data/ENCODE_k562_tf_peaks/bed_filtered_hg19/ENCFF045CFL.bed > out/kzfp_znfs_znf274_encode_ENCFF045CFL_sorted_cluster_acrophase.bed
