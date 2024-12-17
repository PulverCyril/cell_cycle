#! /bin/bash
mkdir out
mkdir out/temp
mkdir out/temp/bigwigs
cd out/temp/bigwigs/

#ZNF274 bigwig UC Davis in K562
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig

#ZNF274_M01 bigwig UC Davis in K562
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig

# ZNF274 ChIP-seq in 293T Michael
wget https://chip-atlas.dbcls.jp/data/hg19/eachData/bw/SRX2512747.bw

#H3K9me3 bigwig UC Davis in K562
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig

#H3K4me3 bigwig in K562
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeSydhHistoneK562H3k4me3bUcdSig.bigWig

#H3K4me1 bigwig in K562
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeSydhHistoneK562H3k4me1UcdSig.bigWig

#H3K27ac bigwig in K562
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig

#H3K27me3 bigwig in K562
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig

#h3k36me3 bigwig in K562:
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig

#ChromHMM state bigwig in K562
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeAwgSegmentationChromhmmK562.bed.gz

#ATAC-seq in K562 from ENCODE, replicate 1 of 3, mapped to hg19 by ChIP-ATLAS:
wget https://chip-atlas.dbcls.jp/data/hg19/eachData/bw/SRX10476524.bw

#CTCF ChIP-seq in K562 from ENCODE
wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeBroadHistoneK562CtcfStdSig.bigWig

#ZNF75A ChIP-exo in 293Ts
wget https://chip-atlas.dbcls.jp/data/hg19/eachData/bw/SRX14890453.bw

#ZNF75D replicate 2 ChIP-exo in 293Ts
wget https://chip-atlas.dbcls.jp/data/hg19/eachData/bw/SRX2512882.bw

#ATRX ChIP-seq in K562:
wget https://chip-atlas.dbcls.jp/data/hg19/eachData/bw/SRX1097312.bw

cd ../../../


# all gene bodies + 1kb 5' for the 5'-most promoter, sorted as in ../data/1807_hg19_ens_coding_genes_body_symbol: 
# H3K9me3. 
#Strand is apparently taken care automatically by deeptools.
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R data/1807_hg19_ens_coding_genes_body_symbol.bed \
 -b 1000 -a 0 \
 --sortRegions no\
 -out out/temp/1807_hg19_ens_coding_genes_body_symbol_H3K9me3.tab.gz \
 --outFileNameMatrix out/temp/1807_hg19_ens_coding_genes_body_symbol_H3K9me3.tab \
 --outFileSortedRegions out/temp/1807_hg19_ens_coding_genes_body_symbol_H3K9me3_ID.tsv
 
 
# ZNF274 in K562 
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R data/1807_hg19_ens_coding_genes_body_symbol.bed \
 -b 1000 -a 0 \
 --sortRegions no\
 -out out/temp/1807_hg19_ens_coding_genes_body_symbol_znf274k562m01.tab.gz \
 --outFileNameMatrix out/temp/1807_hg19_ens_coding_genes_body_symbol_znf274k562m01.tab \
 --outFileSortedRegions out/temp/1807_hg19_ens_coding_genes_body_symbol_znf274k562m01_ID.tsv

# not done here: RIF1 binding

# all TSS as sorted in tss_clustered.bed
# H3K4me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me3bUcdSig.bigWig \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_h3k4me3.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_h3k4me3.tab \
 --outFileSortedRegions out/temp/tss_clustered_h3k4me3_ID.tsv
 
 # ATAC on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX10476524.bw \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_atac.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_atac.tab \
 --outFileSortedRegions out/temp/tss_clustered_atac_ID.tsv

#H3K27ac
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_h3k27ac.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_h3k27ac.tab \
 --outFileSortedRegions out/temp/tss_clustered_h3k27ac_ID.tsv
 
#H3K4me1
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me1UcdSig.bigWig \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_h3k4me1.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_h3k4me1.tab \
 --outFileSortedRegions out/temp/tss_clustered_h3k4me1_ID.tsv

#H3K9me3
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_h3k9me3.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_h3k9me3.tab \
 --outFileSortedRegions out/temp/tss_clustered_h3k9me3_ID.tsv

# H3K9me3 high resolution centered on promoters +- 1.5kb. For the K9 trough plot in the heatmap
computeMatrix reference-point \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R data/tss_clustered.bed \
 -a 1500 -b 1500 \
 --referencePoint center \
 --binSize 10 \
 --averageTypeBins median \
 --sortRegions no \
 --missingDataAsZero \
-out out/temp/tss_clustered_h3k9me3_domains.tab.gz \
--outFileNameMatrix out/temp/tss_clustered_h3k9me3_domains.tab \
 --outFileSortedRegions out/temp/tss_clustered_h3k9me3_domains_ID.tsv

#H3K27me3
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_h3k27me3.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_h3k27me3.tab \
 --outFileSortedRegions out/temp/tss_clustered_h3k27me3_ID.tsv

# ZNF274 ENCODE replicate 2 (m01)
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_znf274m01.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_znf274m01.tab \
 --outFileSortedRegions out/temp/tss_clustered_znf274m01_ID.tsv

#ZNF75A
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX14890453.bw \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_znf75a.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_znf75a.tab \
 --outFileSortedRegions out/temp/tss_clustered_znf75a_ID.tsv

#ZNF75D
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX2512882.bw \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_znf75d.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_znf75d.tab \
 --outFileSortedRegions out/temp/tss_clustered_znf75d_ID.tsv

#CTCF
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562CtcfStdSig.bigWig \
 -R data/tss_clustered.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/tss_clustered_ctcf.tab.gz \
 --outFileNameMatrix out/temp/tss_clustered_ctcf.tab \
 --outFileSortedRegions out/temp/tss_clustered_ctcf_ID.tsv

# KZFPs and ZFPs: histone marks on Zinc finger encoding arrays
# H3K9me3
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/OutOfFrameZFA_KZFPs_ZFPs.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/OutOfFrameZFA_KZFPs_ZFPs_h3k9me3.tab.gz \
 --outFileNameMatrix out/temp/OutOfFrameZFA_KZFPs_ZFPs_h3k9me3.tab \
 --outFileSortedRegions out/temp/OutOfFrameZFA_KZFPs_ZFPs_h3k9me3_ID.tsv

#h3k27me3
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig \
 -R out/OutOfFrameZFA_KZFPs_ZFPs.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/OutOfFrameZFA_KZFPs_ZFPs_h3k27me3.tab.gz \
 --outFileNameMatrix out/temp/OutOfFrameZFA_KZFPs_ZFPs_h3k27me3.tab \
 --outFileSortedRegions out/temp/OutOfFrameZFA_KZFPs_ZFPs_h3k27me3_ID.tsv

#ATAC
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX10476524.bw \
 -R out/OutOfFrameZFA_KZFPs_ZFPs.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/OutOfFrameZFA_KZFPs_ZFPs_atac.tab.gz \
 --outFileNameMatrix out/temp/OutOfFrameZFA_KZFPs_ZFPs_atac.tab \
 --outFileSortedRegions out/temp/OutOfFrameZFA_KZFPs_ZFPs_atac_ID.tsv

#ZNF274 replicate m01
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R out/OutOfFrameZFA_KZFPs_ZFPs.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf274m01.tab.gz \
 --outFileNameMatrix out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf274m01.tab \
 --outFileSortedRegions out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf274m01_ID.tsv

#ZNF75A
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX14890453.bw \
 -R out/OutOfFrameZFA_KZFPs_ZFPs.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf75a.tab.gz \
 --outFileNameMatrix out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf75a.tab \
 --outFileSortedRegions out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf75a_ID.tsv

#ZNF75D
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX2512882.bw \
 -R out/OutOfFrameZFA_KZFPs_ZFPs.bed \
 -a 0 -b 0 \
 --sortRegions no\
 -out out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf75d.tab.gz \
 --outFileNameMatrix out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf75d.tab \
 --outFileSortedRegions out/temp/OutOfFrameZFA_KZFPs_ZFPs_znf75d_ID.tsv


# all genes: histone marks at the TSS
# H3K4me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me3bUcdSig.bigWig \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k4me3.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k4me3.tab
 
# ATAC-seq on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX10476524.bw \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_atac.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_atac.tab

# H3K27Ac on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k27ac.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_acrophase_h3k27ac.tab
 
# H3K4me1 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me1UcdSig.bigWig \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k4me1.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_acrophase_h3k4me1.tab
 
# H3K9me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k9me3.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_acrophase_h3k9me3.tab
 
 
# H3K9me3 wide domains centered on TSS, according to acrophase
# old params
# -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
# -out out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k9me3_domains.tab.gz \
# --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_acrophase_h3k9me3_domains.tab
########## Not too bad, yields a hole of zero signal which corresponds to zero K9 in the middle, which seems to become smaller in M-G1, for late replicating genes
computeMatrix reference-point \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/temp/tss_rhythmic_genes_sorted_acrophase.bed \
 -a 1500 -b 1500 \
 --referencePoint center \
 --binSize 10 \
 --averageTypeBins median \
 --sortRegions keep \
 --missingDataAsZero \
-out out/temp/tss_rhythmic_genes_sorted_acrophase_h3k9me3_domains.tab.gz \
--outFileNameMatrix out/temp/tss_rhythmic_genes_sorted_acrophase_h3k9me3_domains.tab

plotHeatmap \
 -m out/temp/tss_rhythmic_genes_sorted_acrophase_h3k9me3_domains.tab.gz\
 -out out/tss_rhythmic_genes_sorted_acrophase_h3k9me3_domains.png \
 --heatmapHeight 15  \
 -max 6 \
 --sortRegions keep \
 --plotTitle 'H3K9me3 domains rhythmic genes'
 
#H3K27me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k27me3.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k27me3.tab

# ZNF274 on TSS, ENCODE replicate 2, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_znf274_m01.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_znf274_m01.tab
 
#ZNF75A at TSS in 293Ts (Michael)
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX14890453.bw \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_znf75a.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_znf75a.tab
 
#ZNF75D at TSS in 293Ts (Michael, replicate 2)
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX2512882.bw \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_znf75d.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_znf75d.tab

# CTCF TSS all genes
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562CtcfStdSig.bigWig \
 -R out/tss_all_genes_sorted_acrophase_tssstart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_all_genes_sorted_acrophase_tssstart_ctcf.tab.gz \
 --outFileNameMatrix out/temp/tss_all_genes_sorted_acrophase_tssstart_ctcf.tab

# all rhythmic genes: selecting the top active TSS per gene based on H3K4me3, H3K4me1, H3K27Ac and ATAC-seq

# H3K4me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me3bUcdSig.bigWig \
 -R out/tss_rhythmic_genes_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_rhythmic_genes_strict_sorted_acrophase_h3k4me3.tab.gz \
 --outFileNameMatrix out/temp/tss_rhythmic_genes_strict_sorted_acrophase_h3k4me3.tab
 
# ATAC-seq on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX10476524.bw \
 -R out/tss_rhythmic_genes_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_rhythmic_genes_strict_sorted_acrophase_atac.tab.gz \
 --outFileNameMatrix out/temp/tss_rhythmic_genes_strict_sorted_acrophase_atac.tab

# H3K27Ac on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig \
 -R out/tss_rhythmic_genes_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_rhythmic_genes_strict_sorted_acrophase_h3k27ac.tab.gz \
 --outFileNameMatrix out/temp/tss_rhythmic_genes_strict_sorted_acrophase_h3k27ac.tab
 
# H3K9me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/tss_rhythmic_genes_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_rhythmic_genes_strict_sorted_acrophase_h3k9me3.tab.gz \
 --outFileNameMatrix out/temp/tss_rhythmic_genes_strict_sorted_acrophase_h3k9me3.tab
 
#H3K27me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig \
 -R out/tss_rhythmic_genes_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/tss_rhythmic_genes_strict_sorted_acrophase_h3k27me3.tab.gz \
 --outFileNameMatrix out/temp/tss_rhythmic_genes_strict_sorted_acrophase_h3k27me3.tab


# KZFPS: TSS: selecting the top active TSS per gene based on H3K4me3, H3K4me1 and H3K27Ac #

# H3K4me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me3bUcdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_h3k4me3.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_h3k4me3.tab
 
# H3K4me1 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me1UcdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_h3k4me1.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_h3k4me1.tab

# H3K27Ac on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_h3k27ac.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_h3k27ac.tab
 
# ATAC-seq on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX10476524.bw \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_atac.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_atac.tab
 
 
# CTCF on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562CtcfStdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_ctcf.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_ctcf.tab
 
#ZNF75A at TSS in 293Ts (Michael)
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX14890453.bw \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_znf75a.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_znf75a.tab
 
#ZNF75D at TSS in 293Ts (Michael, replicate 2)
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX2512882.bw \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_znf75d.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_znf75d.tab
 
 
# H3K9me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_h3k9me3.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_h3k9me3.tab

# ZNF274 on TSS, ENCODE replicate 2, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_znf274_m01.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_znf274_m01.tab
 
# KZFPs: ZNFs

# ZNF274 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_strict_sorted_acrophase_znf274.tab.gz \
 --outFileNameMatrix out/temp/kzfp_tss_strict_sorted_acrophase_znf274.tab

 
# H3K9me3 on ZNFs, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_acrophase_h3k9me3.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_acrophase_h3k9me3.tab
 
 
# ZNF274 on ZNFs, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_acrophase_znf274.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_acrophase_znf274.tab
 
 
#ZNF75A at ZNFs, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX14890453.bw \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_acrophase_znf75a.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_acrophase_znf75a.tab
 
#ZNF75D at ZNFs, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX2512882.bw \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_acrophase_znf75d.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_acrophase_znf75d.tab

# ZNF274 on ZNFs, using the UCSC utilitary tool 
bin_external/bigWigAverageOverBed out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig out/kzfp_zfcoordinates_sorted_acrophase.bed out/temp/kzfp_znfs_strict_sorted_acrophase_ZNF274_bigWigAverageOverBed.tab

# ZNF274 on ZNFs, according to coordinates
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_coordinates.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_coordinates_znf274.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_coordinates_znf274.tab 

# Signal at ZNFs for plotting (with 1kb margins on either sides of the ZNFs)
 
# ZNF274 in K562 replicate 1 on ZNFs, according to cluster and acrophase within cluster: 
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed \
 -a 1000 -b 1000 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274.tab

# ZNF274 in K562 replicate 2 on ZNFs, according to cluster and acrophase within cluster: 

computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed \
 -a 1000 -b 1000 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_m01.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_m01.tab


#ZNF274 at ZNFs in 293T (Imbeault) 

computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX2512747.bw \
 -R out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed \
 -a 1000 -b 1000 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_imbeault.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_imbeault.tab

#ZNF274 at ZNFs in 293T (Begnis) 
computeMatrix scale-regions \
 -S data/martina_bigwigs/ZNF274_HA_mean.bw \
 -R out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed \
 -a 1000 -b 1000 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_begnis.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_begnis.tab


#ZNF75A at ZNFs in 293Ts (Michael)
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX14890453.bw \
 -R out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed \
 -a 1000 -b 1000 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf75a.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf75a.tab
 
#ZNF75D at ZNFs in 293Ts (Michael, replicate 2)
computeMatrix scale-regions \
 -S out/temp/bigwigs/SRX2512882.bw \
 -R out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed \
 -a 1000 -b 1000 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf75d.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf75d.tab



#H3K36 on ZNFs, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_acrophase_h3k36me3.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_acrophase_h3k36me3.tab


#H3K9me3 on ZNFs, according to acrophase and cluster
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_h3k9me3.tab.gz \
 --outFileNameMatrix out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_h3k9me3.tab
 

# KRAB-less C2H2 ZNFs, according to acrophase and DNAStart of the ZNF array
# kept strictly at 0 as it serves to pick the ZNF array the with the highest ZNF274 binding per gene, where necessary.
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R out/c2h2s_nokrab_zfcoordinates_sorted_acrophase_DNAStart.bed \
 -a 0 -b 0 \
 --sortRegions keep\
 -out out/temp/c2h2s_nokrab_zfcoordinates_strict_sorted_acrophase_dnastart_znf274m01.tab.gz \
 --outFileNameMatrix out/temp/c2h2s_nokrab_zfcoordinates_strict_sorted_acrophase_dnastart_znf274m01.tab
 
 
# With 1kb margins, as this contains exatctly one ZF array per C2H2-krabless ZNF, for plotting purposes
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R out/c2h2s_nokrab_zfcoordinates_top274_sorted_acrophase.bed \
 -a 1000 -b 1000 \
 --sortRegions keep\
 -out out/temp/c2h2s_nokrab_zfcoordinates_top274_sorted_acrophase_znf274m01.tab.gz \
 --outFileNameMatrix out/temp/c2h2s_nokrab_zfcoordinates_top274_sorted_acrophase_znf274m01.tab

###########################################################

# ZNF274 ChIP-seq at ZNFs, according to acrophase K562 replicate 1
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 4000 -b 4000 \
 --sortRegions keep\
 -out out/temp/znfs_ZNF274_K562_acrophase_signal.tab.gz
 
plotHeatmap \
 -m out/temp/znfs_ZNF274_K562_acrophase_signal.tab.gz\
 -out out/znfs_ZNF274_K562_acrophase_signal.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'ZNF274 ChIP-seq in K562'
 
 
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_signif_cycling_zfcoordinates_sorted_acrophase.bed \
 -a 1000 -b 1000 \
 --sortRegions keep\
 -out out/temp/znfs_ZNF274_K562_signif_cycling_signal_acrophase.tab.gz
 
plotHeatmap \
 -m out/temp/znfs_ZNF274_K562_signif_cycling_signal_acrophase.tab.gz\
 -out out/znfs_ZNF274_K562_signif_cycling_signal_acrophase.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10\
 --plotTitle 'ZNF274 ChIP-seq in K562'
 
# ZNF274 ChIP-seq in K562 replicate 1 at ZNFs, according to MESOR

computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_MESOR.bed \
 -a 4000 -b 4000 \
 --sortRegions keep\
 -out out/temp/znfs_ZNF274_K562_signal_MESOR.tab.gz
 
plotHeatmap \
 -m out/temp/znfs_ZNF274_K562_signal_MESOR.tab.gz\
 -out out/znfs_ZNF274_K562_signal_MESOR.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'ZNF274 ChIP-seq in K562'
 
# ZNF274 ChIP-seq replicate 2 at ZNFs in K562 replicate, according to acrophase 

computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274m01UcdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 4000 -b 4000 \
 --sortRegions keep\
 -out out/temp/znfs_ZNF274_K562_M01_signal_acrophase.tab.gz
 
plotHeatmap \
 -m out/temp/znfs_ZNF274_K562_M01_signal_acrophase.tab.gz\
 -out out/znfs_ZNF274_K562_M01_signal_acrophase.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 20 \
 --plotTitle 'ZNF274 ChIP-seq in K562'
 
# ZNF274 ChIP-seq in K562 replicate 1 at TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/tss_ZNF274_K562_signal_acrophase.tab.gz
 
plotHeatmap \
 -m out/temp/tss_ZNF274_K562_signal_acrophase.tab.gz\
 -out out/tss_ZNF274_K562_signal_acrophase.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'ZNF274 ChIP-seq in K562'
 
# ZNF274 at the KRAB domain, as a negative control for the effect of cell cycle on chromatin opening

# replicate 1
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhTfbsK562Znf274UcdSig.bigWig \
 -R out/kzfp_krabcoordinates_sorted_acrophase.bed \
 -a 4000 -b 4000 \
 --sortRegions keep\
 -out out/temp/KRABs_ZNF274_K562_signal_acrophase.tab.gz
 
plotHeatmap \
 -m out/temp/KRABs_ZNF274_K562_signal_acrophase.tab.gz\
 -out out/KRABs_ZNF274_K562_signal_acrophase.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'ZNF274 ChIP-seq in K562 at KRAB domains'

# H3K9me3 K652 at ZNFs according to acrophase

computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 4000 -b 4000 \
 --sortRegions keep\
 -out out/temp/znfs_h3k9me3_K562_signal_acrophase.tab.gz
 
plotHeatmap \
 -m out/temp/znfs_h3k9me3_K562_signal_acrophase.tab.gz\
 -out out/znfs_h3k9me3_K562_signal_acrophase.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K9me3 ChIP-seq in K562'
 
# H3K9me3 on TSS (for now about twice as many TSS as KZFPs, we should maybe filter by H3K4me1 signal)
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/tss_h3k9me3_K562_signal_acrophase.tab.gz
 
plotHeatmap \
 -m out/temp/tss_h3k9me3_K562_signal_acrophase.tab.gz\
 -out out/tss_h3k9me3_K562_signal_acrophase.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K9me3 at KZFP TSS in K562'
 
# H3K9me3 on TSS, according to MESOR
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k9me3StdSig.bigWig \
 -R out/kzfp_tss_sorted_MESOR.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/tss_h3k9me3_sorted_MESOR.tab.gz
 
plotHeatmap \
 -m out/temp/tss_h3k9me3_sorted_MESOR.tab.gz\
 -out out/tss_h3k9me3_sorted_MESOR.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K9me3 at KZFP TSS in K562'
 
 
# H3K36me3 at ZNFs, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_acrophase.bed \
 -a 4000 -b 4000 \
 --sortRegions keep\
 -out out/temp/znfs_h3k36me3_K562_signal_acrophase.tab.gz
 
plotHeatmap \
 -m out/temp/znfs_h3k36me3_K562_signal_acrophase.tab.gz\
 -out out/znfs_h3k36me3_K562_signal_acrophase.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K36me3 ChIP-seq in K562'
 
 
 # H3K36me3 at ZNFs, according to MESOR
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k36me3StdSig.bigWig \
 -R out/kzfp_zfcoordinates_sorted_MESOR.bed \
 -a 4000 -b 4000 \
 --sortRegions keep\
 -out out/temp/znfs_h3k36me3_K562_signal_MESOR.tab.gz
 
plotHeatmap \
 -m out/temp/znfs_h3k36me3_K562_signal_MESOR.tab.gz\
 -out out/znfs_h3k36me3_K562_signal_MESOR.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K36me3 ChIP-seq in K562'
 

# H3K4me3 on TSS, according to MESOR
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me3bUcdSig.bigWig \
 -R out/kzfp_tss_sorted_MESOR.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_sorted_MESOR_h3k4me3.tab.gz
 
plotHeatmap \
 -m out/temp/kzfp_tss_sorted_MESOR_h3k4me3.tab.gz\
 -out out/kzfp_tss_sorted_MESOR_h3k4me3.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K4me3 at KZFP TSS in K562'
 
# H3K4me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me3bUcdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_sorted_acrophase_h3k4me3.tab.gz
 
plotHeatmap \
 -m out/temp/kzfp_tss_sorted_acrophase_h3k4me3.tab.gz\
 -out out/kzfp_tss_sorted_acrophase_h3k4me3.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K4me3 at KZFP TSS in K562'
 
# H3K4me1 on TSS, according to MESOR
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me1UcdSig.bigWig \
 -R out/kzfp_tss_sorted_MESOR.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_sorted_MESOR_h3k4me1.tab.gz
 
plotHeatmap \
 -m out/temp/kzfp_tss_sorted_MESOR_h3k4me1.tab.gz\
 -out out/kzfp_tss_sorted_MESOR_h3k4me1.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K4me1 at KZFP TSS in K562'
 
 # H3K4me1 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k4me1UcdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_sorted_acrophase_h3k4me1.tab.gz
 
plotHeatmap \
 -m out/temp/kzfp_tss_sorted_acrophase_h3k4me1.tab.gz\
 -out out/kzfp_tss_sorted_acrophase_h3k4me1.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K4me1 at KZFP TSS in K562'
 
 # H3K27ac on TSS, according to MESOR
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig \
 -R out/kzfp_tss_sorted_MESOR.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_sorted_MESOR_h3k27ac.tab.gz
 
plotHeatmap \
 -m out/temp/kzfp_tss_sorted_MESOR_h3k27ac.tab.gz\
 -out out/kzfp_tss_sorted_MESOR_h3k27ac.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K27ac at KZFP TSS in K562'
 
 # H3K27ac on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeBroadHistoneK562H3k27acStdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_sorted_acrophase_h3k27ac.tab.gz
 
plotHeatmap \
 -m out/temp/kzfp_tss_sorted_acrophase_h3k27ac.tab.gz\
 -out out/kzfp_tss_sorted_acrophase_h3k27ac.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K27Ac at KZFP TSS in K562'
 
# H3K27me3 on TSS, according to MESOR
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig \
 -R out/kzfp_tss_sorted_MESOR.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_sorted_MESOR_h3k27me3.tab.gz
 
plotHeatmap \
 -m out/temp/kzfp_tss_sorted_MESOR_h3k27me3.tab.gz\
 -out out/kzfp_tss_sorted_MESOR_h3k27me3.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K27me3 at KZFP TSS in K562'
 
# H3K27me3 on TSS, according to acrophase
computeMatrix scale-regions \
 -S out/temp/bigwigs/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig \
 -R out/kzfp_tss_sorted_acrophase.bed \
 -a 2000 -b 2000 \
 --sortRegions keep\
 -out out/temp/kzfp_tss_sorted_acrophase_h3k27me3.tab.gz
 
plotHeatmap \
 -m out/temp/kzfp_tss_sorted_acrophase_h3k27me3.tab.gz\
 -out out/kzfp_tss_sorted_acrophase_h3k27me3.png \
 --heatmapHeight 15  \
 --sortRegions keep \
 -max 10 \
 --plotTitle 'H3K27me3 at KZFP TSS in K562'
