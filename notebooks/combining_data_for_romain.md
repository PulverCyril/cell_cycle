---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.7
  kernelspec:
    display_name: R
    language: R
    name: ir
---

# Combining data for Romain

This notebook only uses processed data in other scripts and notebooks to gather it for Romain

```R
library(tidyr)
```

For all genes:

-> KZNFs / C2H2 / TF / Other
-> RT score (continuous)
-> H3K9me3 level
-> Rif1 level
-> ZNF274 level
-> H3K9me3 FC upon 274KO
-> RNA FC upon 274KO


loading the gene metadata

```R
gene_metadata = read.table('../data/1807_hg19_ens_coding_genes_body_symbol.bed', 
                            sep = '\t',
                           header = F)
gene_metadata %>% head()
```

```R
gene_metadata = gene_metadata %>% 
    dplyr::rename(chr = 1, 
                  start = 2,
                    end = 3,
                    ensembl = 4,
                    strand = 6,
                    symbol = 7) %>%
    dplyr::select(-V5)

```

```R
gene_metadata %>% head()
```

```R
gene_metadata %>% dim()
```

```R
gene_metadata %>% dplyr::filter(symbol == "ZNF724P")
```

Renaming ZNF724P as ZNF724. Verified that the coordinates matched in the GrCh37 ENSEMBL browser

```R
gene_metadata$symbol[which(gene_metadata$symbol == "ZNF724P")] = "ZNF724"
```

```R
# TF information
kzfps_jonas = read.table('../data/kzfps_jonas.tsv', 
                         header = T, 
                         sep = '\t',
                         quote = "")
kzfps_jonas %>% head()
kzfps_jonas %>% dim()                 
```

```R
kzfps_olga = read.table('../data/from_Olga//Table4Cyril.txt',
                        sep = '\t',
                        header = T)
kzfps_olga %>% head()
```

```R
kzfps_olga$summaryKZFP %>% table()
```

```R
kzfps_olga = kzfps_olga %>% 
    dplyr::filter(Common.Name == "Human",
                  !summaryKZFP %in% c("no.KZFP"))
kzfps_olga %>% dim()
```

```R
# is_KZFP
gene_metadata = gene_metadata %>% 
    dplyr::mutate(isKZFP = ifelse(test = ((symbol %in% kzfps_jonas$assigned_gene) |
                                         symbol %in% kzfps_olga$KRAB.Gene.id),
                                 yes = T,
                                 no = F))
```

```R
gene_metadata %>% dplyr::filter(symbol == "ZNF724P")
```

```R
gene_metadata$isKZFP %>% summary()
```

```R
# is ZFP
tfs = read.csv('../data/DatabaseExtract_v_1.01.csv', 
                 sep = ',',
                header = T,
                check.names = T)
tfs %>% head()
tfs %>% dim()
```

```R
tfs = tfs %>% dplyr::filter(Is.TF. == "Yes")
tfs %>% dim()
```

```R
tfs %>% head()
```

```R
# getting C2H2:
gene_metadata = gene_metadata %>% 
    dplyr::mutate(isC2H2 = ifelse(test = (symbol %in% (tfs %>% dplyr::filter(grepl('C2H2 ZF', DBD)) %>% dplyr::pull(HGNC.symbol))) |
                                      (ensembl %in% (tfs %>% dplyr::filter(grepl('C2H2 ZF', DBD)) %>% dplyr::pull(Ensembl.ID))),
                                 yes = T,
                                 no = F))
```

```R
gene_metadata$isC2H2 %>% table()
```

```R
# getting all TFs

gene_metadata = gene_metadata %>% 
    dplyr::mutate(isTF = ifelse(test = (symbol %in% (tfs %>% dplyr::pull(HGNC.symbol))) |
                                      (ensembl %in% (tfs %>% dplyr::pull(Ensembl.ID))),
                                 yes = T,
                                 no = F))
```

```R
gene_metadata$isTF %>% table()
```

```R
gene_metadata$tf_category = ifelse(gene_metadata$isKZFP, yes = 'KZFP', 
                                    no = ifelse(gene_metadata$isC2H2, yes = 'C2H2 ZF', 
                                                no = ifelse(gene_metadata$isTF, yes = 'TF', no = 'not TF')))
table(gene_metadata$tf_category)
```

```R
# continuous replication timing
replitiming_k562_cont = read.table('../out/RTregions_vs_gene_bodies/RTindex_K562.bed', sep = '\t', header = F)
replitiming_k562_cont %>% head()
```

```R
read_RTindex_over_gene_bodies <- function(bed_path) {
    RTindex = read.table(bed_path, sep = '\t', header = F)
    RTindex_gene_bodies = RTindex %>% 
    dplyr::select(V4, V11, V12) %>% 
    dplyr::rename(ensembl = V4, RTindex = V11, overlap = V12) %>% 
    unique() %>% dplyr::group_by(ensembl) %>% 
    dplyr::mutate(RTindex_avg = (as.double((RTindex %*% overlap)))/sum(overlap)) %>%
    dplyr::select(ensembl, RTindex_avg) %>% unique()    
}
```

```R
RTindexK562_gene_bodies = read_RTindex_over_gene_bodies('../out/RTregions_vs_gene_bodies/RTindex_K562.bed')
```

```R
RTindexK562_gene_bodies %>% head()
```

```R
gene_metadata = gene_metadata %>% 
    dplyr::left_join(., RTindexK562_gene_bodies) %>%
    dplyr::rename(RTindex_K562 = RTindex_avg)
```

```R
RTdiff_RIF1KO = read_RTindex_over_gene_bodies('../out/RTregions_vs_gene_bodies/RTdiff_HCT_RIF1KO_vs_WT.bed')
```

```R
gene_metadata = gene_metadata %>% 
    dplyr::left_join(., RTdiff_RIF1KO) %>%
    dplyr::rename(RTdiff_RIF1KO = RTindex_avg)
```

```R
RTdiff_ZNF274KO = read_RTindex_over_gene_bodies('../out/RTregions_vs_gene_bodies/RTdiffRT_diff_293T_ZNF274KO_vs_WT.bed')
```

```R
gene_metadata = gene_metadata %>% 
    dplyr::left_join(., RTdiff_ZNF274KO) %>%
    dplyr::rename(RTdiff_ZNF274KO = RTindex_avg)
```

```R
gene_metadata %>% head()
```

```R
# genomic marks
sum_signal_computeMatrix <- function(matrix_file, id_file, feature_name, FUN = "sum") {
    id = read.table(id_file, comment.char = "", 
            header = T,
            sep = '\t')
    mat = read.table(matrix_file, sep = '\t', skip=3)
    if(FUN == "sum") {
        id[, feature_name] = mat %>% rowSums(na.rm = T)
    } else if (FUN == "mean") {
            id[, feature_name] = mat %>% rowMeans(na.rm = T)
        }
    return(id)
    }
```

```R
h3k9me3_tss_body = sum_signal_computeMatrix('../out/temp/1807_hg19_ens_coding_genes_body_symbol_H3K9me3.tab',
                                           '../out/temp/1807_hg19_ens_coding_genes_body_symbol_H3K9me3_ID.tsv',
                                           "H3K9me3_tss_body")
```

```R
znf274_tss_body = sum_signal_computeMatrix('../out/temp/1807_hg19_ens_coding_genes_body_symbol_znf274k562m01.tab',
                                           '../out/temp/1807_hg19_ens_coding_genes_body_symbol_znf274k562m01_ID.tsv',
                                           "ZNF274_tss_body")
```

```R
rif1_tss_body = sum_signal_computeMatrix('../out/temp/1807_hg19_ens_coding_genes_body_symbol_hg38_rif1_cutnrun_rep2_log2_FE_over_Input.tab',
                                           '../out/temp/1807_hg19_ens_coding_genes_body_symbol_hg38_rif1_cutnrun_rep2_log2_FE_over_Input_ID.tsv',
                                           "RIF1_tss_body",
                                        FUN = "mean")
```

```R
# merging
gene_metadata = gene_metadata %>% 
    dplyr::left_join(., h3k9me3_tss_body %>% dplyr::select(name, H3K9me3_tss_body), by = c("ensembl"="name")) %>%
    dplyr::left_join(., znf274_tss_body %>% dplyr::select(name, ZNF274_tss_body), by = c("ensembl"="name")) %>%
    dplyr::left_join(., rif1_tss_body %>% dplyr::select(name, RIF1_tss_body), by = c("ensembl"="name"))
```

```R
gene_metadata %>% head()
```

```R
# merging martina's H3K9me3 and RNAseq fold changes and pvals upon ZNF274KO
```

```R
# adding Martina's expression data
ZNF274KO_vs_WT_HEK293T = read.table('../data/martina_processed_files/DE_ZNF274KO_vs_WT_HEK293T.tsv', sep = '\t', header = 1, quote = "")
ZNF274KO_vs_WT_HEK293T %>% head()
ZNF274KO_vs_WT_HEK293T %>% dim()
```

```R
gene_metadata = gene_metadata %>% 
    dplyr::left_join(., 
           ZNF274KO_vs_WT_HEK293T %>% 
            dplyr::select(ensembl, pval, p_adj, foldChange) %>%
            dplyr::rename(p_val_ZNF274KO_vs_WT_HEK293T = pval,
                         p_adj_ZNF274KO_vs_WT_HEK293T = p_adj,
                         foldChange_ZNF274KO_vs_WT_HEK293T = foldChange))
```

## Genomic marks at promoters

```R
# promoters
promoters = read.csv('../data/tss_clustered.bed', sep = '\t', header = F) %>%
    dplyr::select(-V5) %>%
    dplyr::rename(chrom = V1, start = V2, end = V3, ensembl = V4, strand = V6)

promoters %>% head()
```

```R
# loading epigenomic data at the tss:
marks = list(CTCF = "ctcf", ZNF75D = "znf75d", ZNF75A = "znf75a", ZNF274 = "znf274m01", H3K27me3 = "h3k27me3", H3K9me3 = "h3k9me3", H3K27ac = "h3k27ac", ATAC = "atac", H3K4me3 = "h3k4me3")
```

```R
for(m in names(marks)) {
    mark_df = read.table(paste0('../out/temp/tss_clustered_', marks[[m]], '_ID.tsv'), sep = '\t', skip=1, header = F) %>%
        dplyr::select(c(1, 2, 3, 4, 6)) %>%
        dplyr::rename(chrom = V1, start = V2, end = V3, ensembl = V4, strand = V6)
    mark_at_promoters = read.table(paste0('../out/temp/tss_clustered_', marks[[m]], '.tab'), sep = '\t', skip=3, header = F)
    mark_df[paste0(m, '_tss')] = mark_at_promoters %>% rowSums(na.rm = T)
    promoters = promoters %>% dplyr::left_join(., mark_df)
    }
```

```R
promoters_top_h3k4me3_h3k27ac = promoters %>% dplyr::group_by(ensembl) %>% dplyr::slice_max(order_by = H3K4me3_tss, na_rm = T, n = 1) %>% dplyr::ungroup()
promoters_top_h3k4me3_h3k27ac = promoters_top_h3k4me3_h3k27ac %>% dplyr::group_by(ensembl) %>% dplyr::slice_max(order_by = H3K27ac_tss, na_rm = T, n = 1) %>% dplyr::ungroup()
promoters_top_h3k4me3_h3k27ac %>% dim()
```

```R
promoters_top_h3k4me3_h3k27ac %>% head()
```

```R
# saving top_h3k9me3_h3k27ac promoters
promoters_top_h3k4me3_h3k27ac %>% 
    dplyr::select(chrom, start, end, ensembl, strand) %>% 
    write.table('../out/tables/promoters_top_h3k4me3_h3k27ac.bed', sep = '\t', col.names = F, row.names = F, quote = F)
```

```R
# merging into the gene metadata table
gene_metadata = gene_metadata %>% dplyr::left_join(., 
                                                   promoters_top_h3k4me3_h3k27ac %>% dplyr::rename(start_TSS = start, end_TSS = end))
```

```R
gene_metadata %>% head()
```

## Chromatin marks at ZNFs for KZFPs and C2H2 ZFs

```R
# getting the ZNF coordinates from my KZFP scanning script
zfps = read.csv('../data/human_KZFPTable.csv', sep = '\t')
zfps %>% head(1)
zfps %>% dim()
zfps %>% colnames()
```

```R
# TODO later
#zfps = zfps %>% dplyr::filter(!ensembl %in% gene_metadata$ensembl)
#zfps %>% dim()
```

```R
# TODO later zfps = zfps %>% dplyr::left_join(., gene_metadata %>% dplyr::select(ensembl, tf_category, is_C2H2, is_KZFP, is_TF))
```

```R
chr_conversion = read.table("../data/hg19.chromAlias.txt", sep = '\t', header = T)
chr_conversion %>% head()
```

```R
zfps = zfps %>% dplyr::left_join(., chr_conversion %>% dplyr::select(ucsc, refseq), by = c("chromosome" = "refseq")) %>% dplyr::rename(chr = ucsc)
zfps %>% head()
```

```R
# adding ZNF812, which is NC_000019.9_KZFP_261
zfps$GeneName[which(zfps$X == "NC_000019.9_KZFP_261")] = "ZNF812"
```

```R
zfp_symbols = c(zfps$GeneName %>% unique(), "ZNF496", "ZNF79")
zfp_symbols %>% length()
```

```R
# importing KZFPs from Jonas's list: doi: 10.1101/gr.277722.123
kzfps_jonas = read.table('../data/kzfps_jonas.csv', header = T, sep = ';', quote = "")
kzfps_jonas %>% dim()
```

```R
# looking for duplicated ZFPs
zfp_symbols = zfps %>% dplyr::filter(chr %in% c(paste0("chr", 1:22), "chrX", "chrY")) %>% dplyr::select(GeneName) %>% dplyr::arrange(GeneName) %>% unlist() %>% as.character()
zfp_symbols
```

```R
zfps %>% dplyr::filter(GeneName == "CHAMP1")
```

```R
zfp_symbols[which(zfp_symbols %>% duplicated)]
```

Of course we can drop all "notAnnotated" zfs

```R
zfps = zfps %>% dplyr::filter(GeneName != "notAnnotated")
zfps %>% dim()
```

```R
zfps = zfps %>% dplyr::filter(
                        chr %in% c(paste0("chr", 1:22), "chrX", "chrY")) 
zfps %>% dim()
```

```R
zfp_symbols = zfps %>% dplyr::pull(GeneName)
zfp_symbols[which(zfp_symbols %>% duplicated)]
```

```R
# taking into account the fact that some have the same ZNF but different KRAB domains
zfps = zfps %>% 
    dplyr::select(chr, 
                  Strand, 
                  GeneName, 
                  OutOfFrameZFAStart, 
                  OutOfFrameZFAEnd, 
                  OutOfFrameZFAid,
                  numberOfSingleZFs, 
                  CanonicalScore) %>% 
    unique() %>%
    dplyr::rename(chrom = chr, strand = Strand, symbol = GeneName)
zfps %>% dim()
```

```R
# adding missing KZFPs
znf496_zfcoordinates = list("chr1", "-", "ZNF496", 247463860, 247464369, "ZNF496_znf", 4, NA)
znf79_zfcoordinates = list("chr9", "+", "ZNF79", 130206556, 130207464, "ZNF79_znf", 11, NA)

zfps = rbind(zfps, znf496_zfcoordinates, znf79_zfcoordinates)
zfps %>% head()
```

```R
zfp_symbols = zfps %>% dplyr::pull(symbol)
zfp_symbols[which(zfp_symbols %>% duplicated)]
```

```R
zfps %>% 
    dplyr::filter(symbol == "PRDM15") %>% 
    dplyr::select(symbol, OutOfFrameZFAStart, OutOfFrameZFAEnd, OutOfFrameZFAid, numberOfSingleZFs, CanonicalScore)
```

```R
# exporting all ZNF coordinates as bedfiles for computematrix
znf_arrays = zfps %>% 
    dplyr::select(chrom, 
                  OutOfFrameZFAStart, 
                  OutOfFrameZFAEnd, 
                  OutOfFrameZFAid, 
                  numberOfSingleZFs, 
                  strand)
```

```R
znf_arrays %>% write.table('../out/OutOfFrameZFA_KZFPs_ZFPs.bed', sep = '\t', col.names = F, row.names = F, quote = F)
```

```R
# loading epigenomic data at the tss:
marks = list(ZNF75D = "znf75d", ZNF75A = "znf75a", ZNF274 = "znf274m01", H3K27me3 = "h3k27me3", H3K9me3 = "h3k9me3", ATAC = "atac")
```

```R
for(m in names(marks)) {
    mark_df = read.table(paste0('../out/temp/OutOfFrameZFA_KZFPs_ZFPs_', marks[[m]], '_ID.tsv'), sep = '\t', skip=1, header = F) %>%
        dplyr::select(c(1, 2, 3, 4, 6)) %>%
        dplyr::rename(chrom = V1, OutOfFrameZFAStart = V2, OutOfFrameZFAEnd = V3, OutOfFrameZFAid = V4, strand = V6)
    mark_at_znfs = read.table(paste0('../out/temp/OutOfFrameZFA_KZFPs_ZFPs_', marks[[m]], '.tab'), sep = '\t', skip=3, header = F)
    mark_df[paste0(m, '_znfs')] = mark_at_znfs %>% rowSums(na.rm = T)
    znf_arrays = znf_arrays %>% dplyr::left_join(., mark_df)
    }
```

```R
# importing the histone marks over ZNFs:
znf_arrays %>% head()
```

```R
znf_arrays = znf_arrays %>% dplyr::mutate(OutOfFrameZFALength = OutOfFrameZFAEnd-OutOfFrameZFAStart)
```

```R
znf_arrays %>% head()
```

```R
# adding deltaK9
deltaK9 = read.table('../out/OutOfFrameZFA_KZFPs_KZFPsdeltaK9_274KO_293T.bed', sep = '\t', header = F) %>%
    dplyr::select(-V5, -V7, -V6) %>%
    dplyr::rename(chr = V1, 
                  OutOfFrameZFAStart = V2, 
                  OutOfFrameZFAEnd = V3,
                  OutOfFrameZFAid = V4,
                  peak_start = V8,
                  peak_end = V9,
                  peak_foldChange = V10,
                  dk9_padj = V11) %>%
    dplyr::mutate(peak_foldChange = as.double(peak_foldChange), dk9_padj = as.double(dk9_padj))
deltaK9 %>% head()
```

```R
deltaK9 %>% dplyr::filter(OutOfFrameZFAid == "NC_000012.11_RF_OOFZFA_636")
```

```R
znf_arrays$deltaK9_ZNF274KO_znf = znf_arrays$OutOfFrameZFAid %in% (deltaK9 %>% dplyr::filter(dk9_padj < 0.05, peak_foldChange < 0))$OutOfFrameZFAid
```

```R
znf_arrays$deltaK9_ZNF274KO_znf %>% summary()
```

```R
# computing a single epigenomic score per ZFP and KZFP
zfp_epigenomic_at_znfs = zfps %>% 
    dplyr::left_join(., znf_arrays) %>% 
    dplyr::group_by(symbol) %>%
    dplyr::mutate(OutOfFrameZFALengthTotal = sum(OutOfFrameZFALength)) %>%
    dplyr::mutate(OutOfFrameZFALengthFrac = OutOfFrameZFALength/OutOfFrameZFALengthTotal) %>%
    dplyr::summarize(ZNF75D_znfs = sum(ZNF75D_znfs*OutOfFrameZFALengthFrac),
                ZNF75A_znfs = sum(ZNF75A_znfs*OutOfFrameZFALengthFrac),
                  ZNF274_znfs = sum(ZNF274_znfs*OutOfFrameZFALengthFrac),
                   H3K27me3_znfs = sum(H3K27me3_znfs*OutOfFrameZFALengthFrac),
                   H3K9me3_znfs = sum(H3K9me3_znfs*OutOfFrameZFALengthFrac),
                  ATAC_znfs = sum(ATAC_znfs*OutOfFrameZFALengthFrac),
                     deltaK9_ZNF274KO_znf = sum(deltaK9_ZNF274KO_znf)>0
                 )
```

```R
zfp_epigenomic_at_znfs %>% head()
zfp_epigenomic_at_znfs %>% dim()
```

```R
gene_metadata = gene_metadata %>% dplyr::left_join(., zfp_epigenomic_at_znfs)
```

```R
gene_metadata %>% head()
```

```R
gene_metadata %>% dim()
```

Adding the cluster information

```R
# importing KZFPs from Jonas's list: doi: 10.1101/gr.277722.123
kzfps_jonas = read.table('../data/kzfps_jonas.csv', header = T, sep = ';', quote = "")
kzfps_jonas %>% dim()

```

```R
# KZFPs that are present in the scanning list, but not in Jonas's list
gene_metadata %>% dplyr::filter(isKZFP) %>% dplyr::filter(!symbol %in% kzfps_jonas$assigned_gene) %>% dplyr::pull(symbol)
```

Looks good.

<!-- #region -->


We will have to add cluster information for ZNF496 and ZNF79 ourselves.
<!-- #endregion -->

```R
# adding cluster information
kzfps_jonas %>% colnames()
gene_metadata %>% colnames()
```

```R
# any duplicates in jonas KZFP gene names?
kzfps_jonas %>% dplyr::filter(duplicated(assigned_gene))

```

```R
gene_metadata = gene_metadata %>% 
    dplyr::left_join(., kzfps_jonas %>% 
                         dplyr::select(assigned_gene, age_MA, species, z_C2H2_miss, cluster) %>%
                         dplyr::rename(age.DeTribolet = age_MA), by = c("symbol"="assigned_gene"))
```

```R
gene_metadata %>% head()
```

ZNF79 is not in any cluster, it is isolated on chromosome 9
ZNF496 is in a cluster together with ZNF669, ZNF124 and others.

```R
znf79_idx = which(gene_metadata$symbol == "ZNF79")
znf496_idx = which(gene_metadata$symbol == "ZNF496")
cluster_idx = which(colnames(gene_metadata) == "cluster")
```

```R
gene_metadata[znf79_idx, cluster_idx] = "noCluster"
gene_metadata[znf496_idx, cluster_idx] = "chr1.2"

```

```R
gene_metadata$cluster = ifelse((gene_metadata$isKZFP & is.na(gene_metadata$cluster)), yes = "noCluster", no = gene_metadata$cluster)
```

```R
cluster_names = gene_metadata$cluster %>% unique() %>% sort()
cluster_names
```

```R
kzfp_cluster_order = c(cluster_names[length(cluster_names)],
                       cluster_names[1:2],
                       cluster_names[20:30],
                       cluster_names[3:19],
                       cluster_names[31])
kzfp_cluster_order                  
```

## Adding Olga's KZFP information

```R
kzfps_olga = read.table('../data/from_Olga/Table4Cyril.txt', sep = '\t', header = T)
kzfps_olga$summaryKZFP %>% table()
kzfps_olga %>% head(1)
```

```R
kzfps_olga = kzfps_olga %>% dplyr::filter(Common.Name == "Human", summaryKZFP != "No.KZFP")
kzfps_olga %>% dim()
kzfps_olga %>% head()
```

```R
kzfps_olga %>% dplyr::select(KRAB.Gene.id, KRAB.Tycko.Score.Min) %>% unique() %>% dim()
```

```R
kzfps_olga %>% dplyr::select(KRAB.Gene.id, KRAB.Tycko.Score.Min) %>% unique() %>% dplyr::filter(duplicated(KRAB.Gene.id))
```

There are many KZFPs without gene symbols. We filter them out. What about the others?

```R
kzfps_olga %>% dplyr::select(KRAB.Gene.id, KRAB.Tycko.Score.Min) %>% dplyr::filter(KRAB.Gene.id=="ZNF589")
```

For those with a symbol, we can clearly take the mean tycko score.

```R
kzfps_olga = kzfps_olga %>% 
    dplyr::filter(KRAB.Gene.id != '.') %>%
    dplyr::select(KRAB.Gene.id, KRAB.Tycko.Score.Min) %>% 
    dplyr::group_by(KRAB.Gene.id) %>%
    dplyr::summarize(Tycko_score_avg = mean(KRAB.Tycko.Score.Min))
```

```R
kzfps_olga %>% head()
```

Let's expand the KRAB.Gene.id column by comma

```R
kzfps_olga = kzfps_olga %>% separate_rows(KRAB.Gene.id, sep=",")
kzfps_olga %>% dim()
```

```R
kzfps_olga %>% dplyr::filter(duplicated(KRAB.Gene.id))
```

```R
kzfps_olga = kzfps_olga %>%
    dplyr::group_by(KRAB.Gene.id) %>%
    dplyr::summarize(Tycko_score_avg = mean(Tycko_score_avg))
```

```R
gene_metadata %>% dim()
```

```R
gene_metadata = gene_metadata %>% dplyr::left_join(., kzfps_olga %>% dplyr::rename(symbol = KRAB.Gene.id))
```

### Adding the KZFP ages from corrected file (KZFP_Info.txt) given by Romain. Is this Alex Coudray??

```R
kzfp_info = read.table('../data/KZFP_Info.txt.tsv', header = T, sep = '\t', quote = '')
kzfp_info %>% dim()
kzfp_info %>% head()
```

```R
# renaming ZNF724P to ZNF724
kzfp_info$Name[which(kzfp_info$Name=="ZNF724P")] = "ZNF724"
```

```R
kzfp_info %>% dplyr::filter(duplicated(Name))
```

```R
kzfp_info %>% dplyr::select(Name, Age.Years) %>% unique() %>% dim()
```

```R
kzfp_info %>% dplyr::select(Name, Age.Years) %>% unique() %>% dplyr::select(Name) %>% unique() %>%dim()
```

```R
gene_metadata = gene_metadata %>% dplyr::left_join(., kzfp_info %>% 
                                                       dplyr::select(Name, Age.Years, Age.Species) %>% 
                                                       dplyr::rename(symbol = Name, Age.Years.Coudray = Age.Years, Age.Species.Coudray = Age.Species) %>%
                                                       unique())
```

### Adding the GenOrigin ages


```R
genorigin = read.table('../data/GenOrigin/Homo_sapiens.csv', sep = ',', header = T)
genorigin %>% head()
```

```R
genorigin %>% dplyr::filter(external_gene_name == "ADK")  
```

```R
genorigin %>% dplyr::filter(external_gene_name == "ADK")  
```

```R
genorigin = genorigin %>% dplyr::mutate(gene_age_integer = as.integer(gene_age))
```

```R
genorigin %>% dplyr::filter(is.na(gene_age_integer)) %>% dim()
```

```R
genorigin %>% dplyr::filter(is.na(gene_age_integer), gene_age == ">4290") %>% dim()
```

We can safely reassign those to 4290

```R
genorigin %>% dplyr::pull(gene_age_integer) %>% max(na.rm = T)
```

```R
genorigin$gene_age_integer = ifelse(genorigin$gene_age == ">4290", yes = 4290, no = genorigin$gene_age_integer)
```

```R
genorigin %>% dplyr::filter(is.na(gene_age)) %>% dim()
```

```R
genorigin$gene_age = NULL
genorigin = genorigin %>% dplyr::rename(gene_age = gene_age_integer)
```

```R
genorigin%>% dim()
```

```R
genorigin %>% dplyr::select(external_gene_name) %>% unique() %>% dim()
```

```R
genorigin %>% dplyr::select(ensembl_gene_id) %>% unique() %>% dim()
```

```R
genorigin %>% dplyr::filter(ensembl_gene_id %in% gene_metadata$ensembl) %>% dim()
```

```R
genorigin %>% dplyr::filter(ensembl_gene_id %in% gene_metadata$ensembl) %>% dplyr::pull(gene_age) %>% as.integer() %>% hist()
```

```R
genorigin %>% dplyr::filter(!ensembl_gene_id %in% gene_metadata$ensembl) %>% dplyr::pull(gene_age) %>% as.integer() %>% hist()
```

Problem: most of the ones we don't find are young genes. Why???

```R
genorigin %>% dplyr::filter(!external_gene_name %in% gene_metadata$symbol) %>% dplyr::pull(gene_age) %>% length()
genorigin %>% dplyr::filter(!external_gene_name %in% gene_metadata$symbol) %>% dplyr::pull(gene_age) %>% as.integer() %>% hist()
```

Sort of the same problem...

```R
genorigin %>% dplyr::filter(!ensembl_gene_id %in% gene_metadata$ensembl) %>% dplyr::pull(gene_age) %>% length()

```

```R
genorigin %>% dplyr::filter(duplicated(external_gene_name)) %>% head()
```

```R
genorigin %>% dplyr::filter(external_gene_name=="PDE11A") %>% head()
```

Many of the duplicated gene symbols are because of different ensembl gene id for the same gene, which still receive the exact same evolutionary age.

```R
genorigin %>% dplyr::select(-ensembl_gene_id) %>% unique() %>% dim()
```

```R
genorigin %>% dplyr::select(-ensembl_gene_id) %>% unique() %>%
dplyr::filter(external_gene_name %in% gene_metadata$symbol) %>% dim()
```

```R
genorigin %>% dplyr::select(-ensembl_gene_id) %>% unique() %>%
dplyr::filter(!external_gene_name %in% gene_metadata$symbol) %>% dplyr::arrange(external_gene_name) %>% tail(100)
```

```R
genorigin %>% dplyr::select(-ensembl_gene_id) %>% unique() %>%
dplyr::filter(external_gene_name %in% gene_metadata$symbol, duplicated(external_gene_name))
```

```R
genorigin %>% dplyr::select(-ensembl_gene_id) %>% unique() %>%
dplyr::filter(external_gene_name %in% gene_metadata$symbol) %>%
dplyr::filter(external_gene_name == "WEE2")
```

Many of the genes have several values for their age. We should keep the max value to be conservative.


For instance, how could WEE2 be 3 MYO??

```R
genorigin_for_join = genorigin %>% dplyr::select(-ensembl_gene_id) %>% unique() %>%
dplyr::group_by(external_gene_name) %>% dplyr::slice_max(order_by = gene_age, n = 1)
```

```R
genorigin_for_join %>% dim()
```

```R
genorigin_for_join$external_gene_name %>% unique() %>% length()
```

```R
gene_metadata = gene_metadata %>% dplyr::left_join(., genorigin_for_join %>% dplyr::rename(gene_age.GenOrigin = gene_age,
                                                                                          gene_Interval.GenOrigin = gene_Interval,
                                                                                          gene_branch.GenOrigin = gene_branch), by = c("symbol" = "external_gene_name"))
```

## Adding GenTree ages



```R
genTree = read.table('../data/hg19_ver73_age.tsv', sep = '\t', header = T)
genTree %>% dim()
genTree %>% head()
```

```R
genTree_branches = read.csv('../data/GenTree_branches.csv', header = F)
colnames(genTree_branches) = c("branch", "split.clade", "split.age")
genTree = genTree %>% dplyr::left_join(., genTree_branches)
```

```R
genTree %>% head()
```

```R
genTree %>% dplyr::filter(duplicated(gene))
```

```R
gene_metadata = gene_metadata %>% dplyr::left_join(., genTree %>% dplyr::select(gene, branch, split.clade, split.age) %>%
                                       dplyr::rename(branch.GenTree = branch,
                                                    clade.GenTree = split.clade,
                                                    age.GenTree = split.age),
                                   by = c("ensembl" = "gene"))
```

### Combining age data:
1) Imbeault, 2) genTree (vertebrates), 3) GenOrigin (goes to the origin of life)

```R
genTree_age_max = gene_metadata$age.GenTree %>% max(na.rm = T)
genTree_age_max
```

```R
gene_metadata$age_combined = ifelse(is.na(gene_metadata$Age.Years.Coudray), no = gene_metadata$Age.Years.Coudray, yes = gene_metadata$age.GenTree)

gene_metadata$age_combined = ifelse(is.na(gene_metadata$age_combined), no = gene_metadata$age_combined, yes = gene_metadata$gene_age.GenOrigin)

gene_metadata$age_combined = ifelse((gene_metadata$age_combined == genTree_age_max) & (gene_metadata$gene_age.GenOrigin > genTree_age_max), no = gene_metadata$age_combined, yes = gene_metadata$gene_age.GenOrigin)
```

```R
gene_metadata %>% head()
```

```R
gene_metadata %>% write.table('../out/tables/gene_metadata_augmented.tsv',
                             sep = '\t',
                             quote = F,
                             col.names = T,
                             row.names = F)
```

## Assigning 100kb-wide windows around promoters to rhythmic genes for enrichment in pyTEnrich

```R
rhythm = read.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', header = T)

windows = read.table('../data/windows_merged.bed', sep = '\t', header = F) %>% dplyr::rename(chrom = V1, window_start = V2, window_end = V3, ensembl = V4)

windows %>% head()
windows %>% dim()
```

```R
dir.create('../out/windows_100kb_wide_up_and_downstream_of_rhythmic_promoters/')
temp_df = rhythm %>% dplyr::left_join(windows)

for (p in (rhythm$phase_assigned %>% unique())) {
    temp_df %>% 
        dplyr::filter(phase_assigned == p, padj < 0.05, !is.na(window_start)) %>%
        dplyr::select(chrom, window_start, window_end, ensembl) %>%
        dplyr::arrange(chrom, window_start, window_end, ensembl) %>%
        write.table(paste0('../out/windows_100kb_wide_up_and_downstream_of_rhythmic_promoters/windows_100kb_', gsub("/", "_", p), '.bed'), sep = '\t', col.names = F, row.names = F, quote = F)
    }

```

## Sanity checks

```R
# quick qc check:
gene_metadata %>% head()
```

```R
gene_metadata %>% dplyr::filter(symbol == "ZNF724")
gene_metadata %>% dplyr::filter(symbol == "ZNF519")
gene_metadata %>% dplyr::filter(symbol == "ZNF274")
```

```R
library(ggplot2)
```

```R
gene_metadata %>% dplyr::filter(isKZFP) %>% ggplot(aes(x=Age.Years.Coudray, y = age.DeTribolet))+ geom_point() + ggrepel::geom_label_repel(aes(x = Age.Years.Coudray, y = age.DeTribolet, label = symbol))
```

```R
gene_metadata %>% dplyr::filter(isKZFP) %>% ggplot(aes(x=Age.Years.Coudray, y = age.GenTree))+ geom_point() + ggrepel::geom_label_repel(aes(x = Age.Years.Coudray, y = age.GenTree, label = symbol))
```

```R
gene_metadata %>% dplyr::filter(isKZFP) %>% ggplot(aes(x=Age.Years.Coudray, y = gene_age.GenOrigin))+ geom_point() + ggrepel::geom_label_repel(aes(x = Age.Years.Coudray, y = gene_age.GenOrigin, label = symbol))
```

```R
gene_metadata %>% ggplot(aes(x=age.GenTree, y = gene_age.GenOrigin))+ geom_point()
```

```R
gene_metadata %>% ggplot(aes(x=factor(age.GenTree), y = gene_age.GenOrigin))+ geom_boxplot()
```

```R
gene_metadata %>% dplyr::filter(isTF) %>% ggplot(aes(x=age.GenTree, y = gene_age.GenOrigin))+ geom_point()
```

Even for TFs, the agreement is often not that good between GenTree and GenOrigin

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = gene_age.GenOrigin)) + geom_violin()
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = age.GenTree)) + geom_violin()
```

Gentree stops at the Euteleostome (90% of vertebrates) split, while GenOrigin goes up to the origin of life.

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = log(ZNF274_tss_body+1))) + geom_violin()
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = log(ZNF274_znfs+1))) + geom_violin()
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = RTdiff_ZNF274KO)) + geom_violin()
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = RTdiff_RIF1KO)) + geom_violin()
```

```R
gene_metadata %>% dplyr::filter(isKZFP) %>% ggplot(aes(x = log(ATAC_znfs+1), y = log(ZNF274_znfs+1))) + geom_point()
```

```R
gene_metadata %>% dplyr::filter(isKZFP) %>% ggplot(aes(x = log(H3K9me3_znfs+1), y = log(ZNF274_znfs+1))) + geom_point()
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = log(H3K9me3_tss_body+1))) + geom_violin()
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = RIF1_tss_body)) + geom_violin()
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = RTindex_K562)) + geom_violin()
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = RTdiff_RIF1KO)) + geom_violin() + ggpubr::stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(2, 3), c(2,4)))
```

```R
gene_metadata %>% colnames()
```

```R
library(ggpubr)
```

```R
gene_metadata %>% 
    ggplot(aes(x = RIF1_tss_body, y = RTdiff_RIF1KO)) + 
    geom_point() + facet_wrap(facets = "tf_category", nrow = 2) + 
    xlim(c(-200, 200)) +
    ggpubr::stat_cor(method = "spearman")
```

```R
gene_metadata %>% ggplot(aes(x = tf_category, y = RTdiff_ZNF274KO)) + geom_violin() + ggpubr::stat_compare_means(comparisons = list(c(1, 2), c(1, 3), c(2, 3), c(2,4)))
```

```R
gene_metadata %>% ggplot(aes(x = RIF1_tss_body, y = RTindex_K562)) + geom_point() + facet_wrap("tf_category")
```

```R
gene_metadata %>% dplyr::filter(isKZFP) %>% ggplot(aes(x = log(ZNF274_tss_body+1), y = RTdiff_ZNF274KO)) + geom_point() + ggpubr::stat_cor(method = "spearman")
```

## Tycko score vs rhythmicity?

```R
rhythm = read.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', header = T)
```

```R
to_plot = rhythm %>% dplyr::left_join(., gene_metadata)
```

```R
to_plot %>% dim()
```

```R
to_plot %>% 
    dplyr::filter(isKZFP) %>% 
    ggplot(aes(x = phase_assigned %in% c("S2", "G2", "G2/M"), 
                                                 y = Tycko_score_avg))+
    geom_violin() +
    ggpubr::stat_compare_means()
```

```R
to_plot %>% 
    dplyr::filter(isKZFP, padj < 0.05) %>% 
    ggplot(aes(x = phase_assigned %in% c("S2", "G2", "G2/M"), 
                                                 y = Tycko_score_avg))+
    geom_violin() +
    ggpubr::stat_compare_means()
```

```R
to_plot %>% 
    dplyr::filter(padj < 0.05, !phase_assigned %in% c("S1", "S2", "M")) %>% 
    ggplot(aes(x = phase_assigned %in% c("eG1", "lG1", "G1/S"), y = log(H3K9me3_tss_body +1))) + 
               geom_violin() + 
    geom_boxplot() +
    ggpubr::stat_compare_means()
```

## Exporting the data for Evarist and the shiny app

```R
rhythm = read.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', header = T)
```

```R
rhythm %>% head()
rhythm %>% dim()
```

```R
# merging with the gene_metadata information
rhythm %>% colnames()
```

```R
gene_metadata %>% colnames()
```

```R
rhythm_for_plot = rhythm %>% 
    dplyr::select(ensembl, entrez, symbol, genename, chr, start, end, strand, coding,
                        mesor, amplitude, acrophase, rsquared, pvalue, phase, padj, stars, phase_rounded, delta, phase_assigned,
                        RTindex_avg, rif1) %>%
    dplyr::left_join(., gene_metadata %>% dplyr::select(ensembl, tf_category))
```

```R
rhythm_for_plot %>% dim()
```

```R
expr_corrected_long = read.table('../out/tables/expr_corrected_long.tsv', sep = '\t', header = T)
```

```R
expr_corrected_long %>% head()
```

```R
rhythm_for_plot = rhythm_for_plot %>% dplyr::left_join(., expr_corrected_long %>% dplyr::select(ensembl, sample_id, norm_adj_counts, gate, pca_outlier, phase_corrected))
```

```R
rhythm_for_plot %>% head()
```

```R
rhythm_for_plot %>% colnames()
```

```R
norm.counts_no_outliers_coding = read.table('../out/tables/norm.counts_no_outliers_coding.tsv', sep = '\t', header = T)
sample_order = colnames(norm.counts_no_outliers_coding)
norm.counts_no_outliers_coding %>% head()
```

```R
norm.counts_no_outliers_coding = norm.counts_no_outliers_coding %>% 
    tibble::rownames_to_column(var = "ensembl")
```

```R
norm.counts_no_outliers_coding = norm.counts_no_outliers_coding  %>%
    tidyr::pivot_longer(., cols = 2:ncol(norm.counts_no_outliers_coding)) %>%
    dplyr::rename(sample_id = name, norm.counts = value)
```

```R
norm.counts_no_outliers_coding = norm.counts_no_outliers_coding %>% dplyr::left_join(., expr_corrected_long %>% dplyr::select(-norm_adj_counts))
```

```R
norm.counts_no_outliers_coding %>% head()
```

```R
# adding the fitted curve
```

Expr = MESOR + Amplitude * cos(2*pi/T + acrophase), T = 8, acrophase in radians

```R
rhythm_for_plot %>% dplyr::select(delta) %>% summary()
```

```R
batch = readRDS('../out/tables/bulk_RNAseq_chronogram_batches.RDS')
batch %>% length()
sample_order %>% length()
```

```R
batch_df = data.frame(batch, sample_order) %>% dplyr::rename(sample_id = sample_order)
```

```R
norm.counts_no_outliers_coding = norm.counts_no_outliers_coding %>% dplyr::left_join(., batch_df)
```

```R
rhythm_for_plot = rhythm_for_plot %>% 
    dplyr::left_join(norm.counts_no_outliers_coding) %>%
    dplyr::select(-RTindex_avg, -rif1, -norm_adj_counts)
```

```R
rhythm_for_plot %>% colnames()
rhythm_for_plot %>% dim()
```

```R
rhythm_for_plot %>% write.table('../out/tables/rhythm_for_plot.tsv', col.names = T, row.names = F, quote = F, sep = '\t')
```

```R
g = "ZNF274"
to_plot = rhythm_for_plot %>% 
    dplyr::filter(!pca_outlier) %>%
    dplyr::filter(symbol == g) %>% 
    dplyr::mutate(batch_corrected = norm.counts - delta)

    fit.phase = to_plot$phase %>% unique()
    fit.acr = to_plot$acrophase %>% unique()
    fit.mesor = to_plot$mesor %>% unique()
    fit.ampl = to_plot$amplitude %>% unique()
    fit.padj = to_plot$padj %>% unique()
    fit.phase_ass = to_plot$phase_assigned %>% unique()
    x_estimation = seq(from = -0.5, to = 7.5, by = 0.1)
    estimated_expr = sapply(x_estimation, function(x) fit.mesor + fit.ampl*(cos(2*pi*x/8+fit.phase)))
    fit.df = data.frame(x_estimation, estimated_expr)

    p = to_plot %>% 
    ggplot(., aes(x = factor(phase_corrected, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")), 
           y = norm.counts)) + 
    geom_point(size = 2, aes(col = factor(batch, levels = c(0, 1))), show.legend = F) + 
    scale_color_manual(values = c("black", "grey")) + 
    theme_classic() + 
    geom_curve(data = to_plot %>% dplyr::filter(batch == 1),
             aes(x = factor(phase_corrected, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")),
                                xend = factor(phase_corrected, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")), 
                                                  y = norm.counts, 
            yend = batch_corrected),
             lty = 3,             
            col = "grey") +
    geom_point(data = to_plot %>% dplyr::filter(batch == 1),
    aes(x = factor(phase_corrected, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")),
        y = batch_corrected),
        col = "black") + 
    labs(x = "Phase", y = "expression [log2]") +
    ggtitle(label = paste0(g, 
                           " | peak phase = ", fit.phase_ass,
                           " | adj. p = ", format(fit.padj, scientific = T, digits = 2))) +
    theme(plot.title = element_text(size=11)) +
    geom_line(data = fit.df, aes(x = x_estimation + 1, estimated_expr)) + 
    geom_hline(yintercept = fit.mesor, lty = 2) + 
    geom_segment(data = data.frame(x = fit.acr + 1, xend = fit.acr + 1, 
                                   y = fit.mesor, yend = fit.mesor + fit.ampl),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 lty = 2
                )
p
dev.copy(svg, filename = paste0('../out/figures/cell_cycle_marker_genes/lognormadjcounts/', g, '_lognormadjcounts.svg'), width = 4, height = 4)
dev.off()
```

```R
# example gene: E2F1
g = "E2F1"
to_plot = rhythm_for_plot %>% 
    dplyr::filter(!pca_outlier) %>%
    dplyr::filter(symbol == g) %>% 
    dplyr::mutate(batch_corrected = norm.counts - delta)

    fit.phase = to_plot$phase %>% unique()
    fit.acr = to_plot$acrophase %>% unique()
    fit.mesor = to_plot$mesor %>% unique()
    fit.ampl = to_plot$amplitude %>% unique()
    fit.padj = to_plot$padj %>% unique()
    fit.phase_ass = to_plot$phase_assigned %>% unique()
    x_estimation = seq(from = -0.5, to = 7.5, by = 0.1)
    estimated_expr = sapply(x_estimation, function(x) fit.mesor + fit.ampl*(cos(2*pi*x/8+fit.phase)))
    fit.df = data.frame(x_estimation, estimated_expr)

    p = to_plot %>% 
    ggplot(., aes(x = factor(phase_corrected, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")), 
           y = norm.counts)) + 
    geom_point(size = 2, aes(col = factor(batch, levels = c(0, 1))), show.legend = F) + 
    scale_color_manual(values = c("black", "grey")) + 
    theme_classic() + 
    geom_curve(data = to_plot %>% dplyr::filter(batch == 1),
             aes(x = factor(phase_corrected, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")),
                                xend = factor(phase_corrected, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")), 
                                                  y = norm.counts, 
            yend = batch_corrected),
             lty = 3,             
            col = "grey") +
    geom_point(data = to_plot %>% dplyr::filter(batch == 1),
    aes(x = factor(phase_corrected, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")),
        y = batch_corrected),
        col = "black") + 
    labs(x = "Phase", y = "expression [log2]") +
    ggtitle(label = paste0(g, 
                           " | peak phase = ", fit.phase_ass,
                           " | adj. p = ", format(fit.padj, scientific = T, digits = 2))) +
    theme(plot.title = element_text(size=11)) +
    geom_line(data = fit.df, aes(x = x_estimation + 1, estimated_expr)) + 
    geom_hline(yintercept = fit.mesor, lty = 2) + 
    geom_segment(data = data.frame(x = fit.acr + 1, xend = fit.acr + 1, 
                                   y = fit.mesor, yend = fit.mesor + fit.ampl),
                 aes(x = x, xend = xend, y = y, yend = yend),
                 lty = 2
                )
p
dev.copy(svg, filename = paste0('../out/figures/cell_cycle_marker_genes/lognormadjcounts/', g, '_lognormadjcounts.svg'), width = 4, height = 4)
dev.off()
```

## age of specific KZFPs

```R
#Romain's shortlist: 
romain_list = c("ZNF311", "ZNF334", "ZNF514", "ZNF519", "ZNF586", "ZNF587", "ZNF641", "ZNF669", "ZNF84")
```

```R
gene_metadata %>% dplyr::filter(symbol %in% romain_list) %>% dplyr::select(symbol, age_combined)
```

## Delta H3K9me3

```R
gene_metadata %>% colnames()
```

```R
library(ggplot2)
```

```R
gene_metadata %>% ggplot(aes(x = p_adj_ZNF274KO_vs_WT_HEK293T < 0.05, fill = isKZFP)) + geom_bar(position = "fill")
```

```R
gene_metadata %>% dplyr::filter(isKZFP) %>% ggplot(aes(x = p_adj_ZNF274KO_vs_WT_HEK293T < 0.05, fill = deltaK9_ZNF274KO_znf)) + geom_bar(position = "fill")
```
