#! /bin/bash

mkdir data/ENCODE_k562_tf_peaks/bed/
cd data/ENCODE_k562_tf_peaks/bed/
xargs -n 1 curl -O -L < ../files.txt

cd bed
gunzip *.gz

# filtering peaks to keep only those with q-value < 0.05, i.e. -log10(q) > 1.301

mkdir ../bed_filtered/

for BED in *.bed; do
    awk '{if ($9 > 1.301) print $0;}' $BED > ../bed_filtered/$BED;
done

# liftover to hg19. These are promoters, they should be well annotated
cd ../../../

mkdir data/ENCODE_k562_tf_peaks/bed_filtered_hg19

for BEDPATH in data/ENCODE_k562_tf_peaks/bed_filtered/*; do
    BED=$(basename ${BEDPATH})
    bin_external/liftOver -bedPlus=6 $BEDPATH data/hg38ToHg19.over.chain.gz data/ENCODE_k562_tf_peaks/bed_filtered_hg19/$BED out/temp/unlifted.bed;
done


for BED in *data/.bed; do
done

# computing the intersect with promoters.

mkdir out/encode_x_promoters/

for BEDPATH in data/ENCODE_k562_tf_peaks/bed_filtered_hg19/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a data/tss_clustered.bed -b $BEDPATH > out/encode_x_promoters/$BED;
done


# same for najafabadi and schmitges chip-seq data
mkdir out/schmitges_x_promoters/
for BEDPATH in data/schmitges_kzfp_chips/hg19_peaks_chipatlas/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a data/tss_clustered.bed -b $BEDPATH > out/schmitges_x_promoters/$BED;
done

mkdir out/najafabadi_x_promoters/
for BEDPATH in data/najafabadi_kzfp_chips/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a data/tss_clustered.bed -b $BEDPATH > out/najafabadi_x_promoters/$BED;
done

# same for schmitges chip-seq data, all ZNFs
mkdir out/schmitges_all_znfs_x_promoters/
for BEDPATH in data/schmitges_all_znfs_chips/hg19_peaks_chipatlas/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a data/tss_clustered.bed -b $BEDPATH > out/schmitges_all_znfs_x_promoters/$BED;
done

# extracting the tar archive
tar -xvf data/najafabadi_all_znfs_chips/GSE58341_RAW.tar -C data/najafabadi_all_znfs_chips/hg19_peaks_GEO/
rm data/najafabadi_all_znfs_chips/hg19_peaks_GEO/*summits*
gunzip data/najafabadi_all_znfs_chips/hg19_peaks_GEO/*.gz

# filtering out peaks with less than MACS 80
mkdir data/najafabadi_all_znfs_chips/hg19_peaks_80/
for BEDPATH in data/najafabadi_all_znfs_chips/hg19_peaks_GEO/*; do
    BED=$(basename ${BEDPATH})
    awk '{if($5 >= 80) {print; next}}' $BEDPATH > data/najafabadi_all_znfs_chips/hg19_peaks_80/$BED;
done

# intersecting with promoters
mkdir out/najafabadi_all_znfs_x_promoters/
for BEDPATH in data/najafabadi_all_znfs_chips/hg19_peaks_80/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a data/tss_clustered.bed -b $BEDPATH > out/najafabadi_all_znfs_x_promoters/$BED;
done



# Computing intersect with ZNFs of KZFPs:

mkdir out/encode_x_kzfp_znfs/

for BEDPATH in data/ENCODE_k562_tf_peaks/bed_filtered_hg19/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b $BEDPATH > out/encode_x_kzfp_znfs/$BED;
done

# same for najafabadi and schmitges chip-seq data
mkdir out/schmitges_x_kzfp_znfs/
for BEDPATH in data/schmitges_kzfp_chips/hg19_peaks_chipatlas/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b $BEDPATH > out/schmitges_x_kzfp_znfs/$BED;
done

mkdir out/najafabadi_x_kzfp_znfs/
for BEDPATH in data/najafabadi_kzfp_chips/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b $BEDPATH > out/najafabadi_x_kzfp_znfs/$BED;
done


# same for Schmitges chip-seq datasets, without subsetting on KZFPs:
mkdir out/schmitges_x_all_znfs/
for BEDPATH in data/schmitges_all_znfs_chips/hg19_peaks_chipatlas/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b $BEDPATH > out/schmitges_x_kzfp_znfs/$BED;
done


# intersecting with promoters
mkdir out/najafabadi_all_znfs_x_all_znfs/
for BEDPATH in data/najafabadi_all_znfs_chips/hg19_peaks_80/*; do
    BED=$(basename ${BEDPATH})
    bedtools intersect -wo -a out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed -b $BEDPATH > out/najafabadi_all_znfs_x_all_znfs/$BED;
done
