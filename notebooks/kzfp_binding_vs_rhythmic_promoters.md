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

# Are there evolutionarily recent TFs enriching at rhyhtmic promoters, and in a phase-specific way?

<!-- #region -->
Promoters were defined as clusters of TSS annotated in ensembl (Pulver et al., Genome Biology 2023)


## Conventions
We create column prefixes in the following way:
- `promoter_` indicates that we have something to do with the promoter overlapping the peak
- `TF_` indicates that we have something to do with the TF that was ChIP-ed
- `peak_` indicates that we have something to do with the peak itself


Scripts `download_Schmitges_peaks_all_ZNFs.sh` (Schmitges data), `download_ENCODE_TF_peaks.sh` (ENCODE data, intersect peaks vs promoters for ENCODE, Schmitges and Najafabadi), `intersect_KZFPs_cycling_TSS.sh` (Tronolab data) must be run prior to this notebook.

<!-- #endregion -->

```R
library(tidyverse)
library(R.utils) # to count lines of the filtered and liftovered bedfiles, to know the total number of peaks
```

## Loading the ENCODE TFs x promoters bedtools intersect outputs

```R
df_list = list()
```

```R
encode_tfs_metadata = read.table('../data/ENCODE_k562_tf_peaks/experiment_report_2024_2_19_10h_10m.tsv', skip = 1, header = T, sep = '\t')
```

```R
i = 0
for (bed in list.files('../out/encode_x_promoters/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    p = paste0('../out/encode_x_promoters/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        temp_df = read.table(paste0('../out/encode_x_promoters/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        temp_df$TF_symbol = encode_tfs_metadata[grepl(filename, encode_tfs_metadata$Files) %>% which(), "Target.gene.symbol"]
        temp_df$peaks_total = R.utils::countLines(paste0('../data/ENCODE_k562_tf_peaks/bed_filtered_hg19/', bed))
        i = i+1

        df_list[[i]] = temp_df}
    }
```

```R
peaks_x_promoters_encode = do.call(rbind, df_list)
```

```R
peaks_x_promoters_encode$cell_line = "K562"
peaks_x_promoters_encode$publication = "ENCODE"
```

```R
peaks_x_promoters_encode %>% head() %>% print()
```

```R
# identifying the overlap column: 
min(897329, 896245) - max(895467, 895935)
```

```R
peaks_x_promoters_encode = peaks_x_promoters_encode %>% dplyr::select(-V5, -V10, -V11, -V12, -V14, -V15, -V16)
peaks_x_promoters_encode %>% dim()
```

```R
colnames(peaks_x_promoters_encode) = c("promoter_chr", "promoter_start", "promoter_end", "promoter_ensembl", "promoter_strand", 
                                      "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "filename", "TF_symbol", "peaks_total", "cell_line", "publication")
peaks_x_promoters_encode %>% head()
peaks_x_promoters_encode %>% dim()
```

## Loading the Tronolab ChIP-exo and ChIP-seq peaks x promoters bedtools intersect output

```R
peaks = read.table('../data/hg19_peaks_un_filt_macs80.bed', sep = '\t', header = F, ) %>% dplyr::select(-V6)
```

```R
colnames(peaks) = c("peak_chr", "peak_start", "peak_end", "peak_ID", "peak_score")
```

```R
# extracting KZFP names
peaks$TF_symbol = sapply(peaks$peak_ID, function(x) return((str_split(x, "_", simplify = T))[1])) %>% as.vector()
```

```R
peaks %>% head()
```

```R
peaks$TF_symbol %>% unique() %>% length()
```

```R
all_tss_hg19 = read.table('../data/tss_clustered.bed', sep = '\t')
all_tss_hg19 %>% head()
```

```R
# counting all peaks per KZFP
peaks = peaks %>% dplyr::group_by(TF_symbol) %>% dplyr::mutate(TF_peaks_total=n()) %>% ungroup()
```

```R
peaks %>% head()
```

```R
peaks %>% dim()
```

```R
# annotating low confidence ChIP-seq
peaks$peak_lowconf = grepl("lowconf", peaks$peak_ID)
```

```R
peaks_x_promoters_tronolab = read.table('../out/promoters_hg19_genes_vs_kzfp_peaks.bed', sep = '\t', header = F)
```

```R
peaks_x_promoters_tronolab %>% head()
```

```R
# identifying the overlap column:
min(861618, 859795) - max(859760, 859146)
# V13
```

```R
peaks_x_promoters_tronolab %>% dim()
```

```R
peaks_x_promoters_tronolab = peaks_x_promoters_tronolab %>% dplyr::select(-V5, -V12)
```

```R
colnames(peaks_x_promoters_tronolab) = c("promoter_chr", "promoter_start", "promoter_end", "promoter_ensembl", "promoter_strand", "peak_chr", "peak_start", "peak_end", "peak_ID", "peak_score", "peak_overlap")
```

```R
peaks_x_promoters_tronolab %>% head()
```

```R
any(!grepl("sampled_peak", peaks_x_promoters_tronolab$peak_ID))
```

```R
# adding the peak information
peaks_x_promoters_tronolab = peaks_x_promoters_tronolab %>% dplyr::left_join(., peaks %>% dplyr::select(peak_ID, TF_symbol, peak_lowconf, TF_peaks_total))
```

```R
# recovering the filename (i.e name of the chip itself)
peaks_x_promoters_tronolab$filename = sapply(peaks_x_promoters_tronolab$peak_ID, function(x) str_split(x, "_sampled_peak_")[[1]][[1]])
```

```R
peaks_x_promoters_tronolab %>% head()
```

## Loading the Schmitges ChIP-seq peaks x promoters bedtools intersect output

```R
# adding the schmitges peaks at promoters:
#schmitges_metadata = read.table('../data/schmitges_kzfp_chips//schmitges_kzfps.txt', sep = '', header = F, col.names = c("sample", "filename"))
schmitges_metadata = read.table('../data/schmitges_all_znfs_chips/schmitges_all_znfs.txt', sep = '', header = F, col.names = c("sample", "filename"))
schmitges_metadata$TF_symbol = schmitges_metadata$sample %>% sapply(., FUN = function(x) str_split(x, "_")[[1]][[1]])
schmitges_metadata %>% head()
```

```R
i = 0
df_list = list()
#for (bed in list.files('../out/schmitges_x_promoters/', pattern = "*.bed")) {
for (bed in list.files('../out/schmitges_all_znfs_x_promoters/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    #p = paste0('../out/schmitges_x_promoters/', bed)
    p = paste0('../out/schmitges_all_znfs_x_promoters/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        #temp_df = read.table(paste0('../out/schmitges_x_promoters/', bed), sep = '\t', header = F)
        temp_df = read.table(paste0('../out/schmitges_all_znfs_x_promoters/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        #temp_df$TF_symbol = schmitges_metadata %>% can be done afterwards
        #temp_df$peaks_total = R.utils::countLines(paste0('../data/schmitges_kzfp_chips/hg19_peaks_chipatlas/', bed))
        temp_df$peaks_total = R.utils::countLines(paste0('../data/schmitges_all_znfs_chips/hg19_peaks_chipatlas/', bed))
        i = i+1

        df_list[[i]] = temp_df}
    }
```

```R
peaks_x_promoters_schmitges = do.call(rbind, df_list)
```

```R
peaks_x_promoters_schmitges %>% head() %>% print()
```

```R
peaks_x_promoters_schmitges$cell_line = "HEK293T"
peaks_x_promoters_schmitges$publication = "schmitges"
peaks_x_promoters_schmitges = peaks_x_promoters_schmitges %>% dplyr::left_join(., schmitges_metadata %>% dplyr::select(filename, TF_symbol))
peaks_x_promoters_schmitges %>% head() %>% print()
```

```R
peaks_x_promoters_schmitges = peaks_x_promoters_schmitges %>% dplyr::select(-V5, -V10, -V12, -V13, -V14, -V15, -V16)
```

```R
colnames(peaks_x_promoters_schmitges) = c("promoter_chr", "promoter_start", "promoter_end", "promoter_ensembl", "promoter_strand", 
                                      "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "filename", "peaks_total", "cell_line", "publication", "TF_symbol")
peaks_x_promoters_schmitges %>% head()
peaks_x_promoters_schmitges %>% dim()
```

## Loading the Najafabadi ChIP-seq peaks x promoters bedtools intersect output

```R
# adding the najafabadi peaks at promoters
i = 0
df_list = list()
#for (bed in list.files('../out/najafabadi_x_promoters/', pattern = "*.bed")) {
for (bed in list.files('../out/najafabadi_all_znfs_x_promoters/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    #p = paste0('../out/najafabadi_x_promoters/', bed)
    p = paste0('../out/najafabadi_all_znfs_x_promoters/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        #temp_df = read.table(paste0('../out/najafabadi_x_promoters/', bed), sep = '\t', header = F)
        temp_df = read.table(paste0('../out/najafabadi_all_znfs_x_promoters/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        temp_df$sample = str_split(filename, "_")[[1]][[1]]
        #temp_df$peaks_total = R.utils::countLines(paste0('../data/najafabadi_kzfp_chips/', bed))
        temp_df$peaks_total = R.utils::countLines(paste0('../data/najafabadi_all_znfs_chips/hg19_peaks_80/', bed))
        i = i+1

        df_list[[i]] = temp_df}
    }
```

```R
peaks_x_promoters_najafabadi = do.call(rbind, df_list)
```

```R
peaks_x_promoters_najafabadi %>% head()
```

```R
peaks_x_promoters_najafabadi$cell_line = "HEK293T"
peaks_x_promoters_najafabadi$publication = "najafabadi"
peaks_x_promoters_najafabadi %>% head() %>% print()
```

```R
peaks_x_promoters_najafabadi = peaks_x_promoters_najafabadi %>% dplyr::select(-V5, -V10)
peaks_x_promoters_najafabadi %>% head()
```

```R
najafabadi_metadata = read.table('../data/najafabadi_all_znfs_chips/najafabadi_metadata.txt', sep = ' ', header = F) %>% 
    dplyr::rename(symbol = V1, sample = V2)
najafabadi_metadata %>% head()
```

```R
peaks_x_promoters_najafabadi = peaks_x_promoters_najafabadi %>% 
    dplyr::left_join(., najafabadi_metadata)
```

```R
peaks_x_promoters_najafabadi %>% head()
```

```R
colnames(peaks_x_promoters_najafabadi) = c("promoter_chr", "promoter_start", "promoter_end", "promoter_ensembl", "promoter_strand", 
                                      "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "filename", "sample", "peaks_total", "cell_line", "publication", "TF_symbol")
peaks_x_promoters_najafabadi %>% head()
peaks_x_promoters_najafabadi %>% dim()
```

## Concatenating the bedtools intersect outputs from all datasets

```R
# concatenating with the encode data
peaks_x_promoters_tronolab$cell_line = "HEK293T"
peaks_x_promoters_tronolab$publication = "Tronolab"
peaks_x_promoters = rbind(peaks_x_promoters_encode %>% dplyr::select("promoter_chr", "promoter_start", "promoter_end", "promoter_ensembl", "promoter_strand", 
                                             "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "TF_symbol", "cell_line", "publication", "peaks_total", "filename"),
      peaks_x_promoters_tronolab %>% dplyr::rename(peaks_total = TF_peaks_total) %>% dplyr::select("promoter_chr", "promoter_start", "promoter_end", "promoter_ensembl", "promoter_strand", 
                                             "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "TF_symbol", "cell_line", "publication", "peaks_total", "filename"),
                         peaks_x_promoters_najafabadi %>% dplyr::select("promoter_chr", "promoter_start", "promoter_end", "promoter_ensembl", "promoter_strand", 
                                             "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "TF_symbol", "cell_line", "publication", "peaks_total", "filename"),
                         peaks_x_promoters_schmitges %>% dplyr::select("promoter_chr", "promoter_start", "promoter_end", "promoter_ensembl", "promoter_strand", 
                                             "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "TF_symbol", "cell_line", "publication", "peaks_total", "filename"))
```

```R
# what are TF symbols without associated gene?
# RBM14,RBM14-RBM4 -> RBM14
# U2AF1L5,U2AF1  -> U2AF1
# COMMD3-BMI1,BMI1 -> BMI1
# ZNF724P -> ZNF724

peaks_x_promoters$TF_symbol[which(peaks_x_promoters$TF_symbol == "RBM14,RBM14-RBM4")] = "RBM4"
peaks_x_promoters$TF_symbol[which(peaks_x_promoters$TF_symbol == "U2AF1L5,U2AF1")] = "U2AF1"
peaks_x_promoters$TF_symbol[which(peaks_x_promoters$TF_symbol == "COMMD3-BMI1,BMI1")] = "BMI1"
peaks_x_promoters$TF_symbol[which(peaks_x_promoters$TF_symbol == "ZNF724P")] = "ZNF724"
```

```R
peaks_x_promoters %>% dim()
peaks_x_promoters$TF_symbol %>% unique() %>% length()
```

```R
# computing promoters and peak length
peaks_x_promoters = peaks_x_promoters %>% dplyr::mutate(promoter_length = promoter_end-promoter_start, peak_length = peak_end-peak_start)
```

```R
peaks_x_promoters %>% head()
```

Counting TFs for which we have binding data, and the number of genes where they bind


```R
peaks_x_promoters %>% dplyr::pull(TF_symbol) %>% unique() %>% length()
```

```R
peaks_x_promoters %>% dplyr::pull(promoter_ensembl) %>% unique() %>% length()
```

```R
# adding the gene rhythmicity information, both for KZFPs and the target genes
rhythm_all_genes = read.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', header = 1)
rhythm_all_genes %>% head()
```

```R
df_to_join = rhythm_all_genes %>% 
                         dplyr::rename(gene_chr = chr, 
                                       gene_start = start, 
                                       gene_end = end, 
                                       gene_strand = strand)
colnames(df_to_join) = paste0("promoter_", colnames(df_to_join))

peaks_x_promoters = peaks_x_promoters %>% 
    dplyr::left_join(x=., y=df_to_join)
peaks_x_promoters %>% colnames()
```

```R
peaks_x_promoters %>% dplyr::pull(promoter_ensembl) %>% unique() %>% length()
peaks_x_promoters %>% dplyr::pull(promoter_symbol) %>% unique() %>% length()
```

```R
peaks_x_promoters %>% dplyr::filter(is.na(promoter_symbol)) %>% head()
```

This is due to many ensembl ID having no associated HGNC symbol. The reverse is not true. 


Here already we have a huge discrepancy between promoter symbols and promoter ensembl ID...

```R
rhythm_all_genes %>% dplyr::pull(ensembl) %>% unique() %>% length()
rhythm_all_genes %>% dplyr::pull(symbol) %>% unique %>% length()
```

```R
rhythm_all_genes %>% dplyr::filter(is.na(symbol))
```

In contrast, there are no genes with NA symbol and/or HGNC symbol in the rhythm data.

```R
peaks_x_promoters %>% dplyr::filter(is.na(promoter_symbol)) %>% head()
```

```R
peaks_x_promoters %>% ncol()
```

```R
# adding the TF cycling information
df_to_join = rhythm_all_genes
colnames(df_to_join) = paste0("TF_", colnames(df_to_join))
peaks_x_promoters = peaks_x_promoters %>% 
    dplyr::left_join(x =., y = df_to_join)
peaks_x_promoters %>% colnames()
```

```R
# loading KZFP information
kzfps_cycling_info = read.csv('../out/tables/kzfps_cycling_info.tsv', header = T, sep = '\t')
kzfps_cycling_info %>% head()
```

Building a matrix of TF peaks per promoter. Note that we use the promoter_ensembl column to group_by, meaning that we take the background number of genes to be ~18000

```R
npeaks_per_promoter = peaks_x_promoters %>% dplyr::group_by(promoter_ensembl, TF_symbol) %>% dplyr::summarize(nPeaks = n()) %>% spread(., TF_symbol, nPeaks)
```

```R
npeaks_per_promoter %>% head()
```

```R
npeaks_per_promoter %>% dim()
```

```R
npeaks_per_promoter_mat = as.matrix(npeaks_per_promoter %>%tibble::column_to_rownames("promoter_ensembl") %>% as.data.frame())
```

```R
npeaks_per_promoter_mat_cleaned = ifelse(is.na(npeaks_per_promoter_mat), yes = 0, no = npeaks_per_promoter_mat)
```

```R
npeaks_per_promoter_mat_cleaned %>% head()
```

```R
npeaks_per_promoter_mat_cleaned %>% write.table('../out/tables/gene_promoters_x_ChIPseq_ENCODE_Trono_hughes.tsv', row.names = T, col.names = T, quote = F, sep = '\t')
```

```R
# TFs to targets
peaks_x_promoters %>% dplyr::select(TF_symbol, promoter_ensembl, promoter_symbol) %>% unique() %>% write.table('../out/tables/TFs_to_target_promoters.tsv', sep = '\t', row.names = F, col.names = T, quote = F)
```

### Exporting rhythmic promoters per phase, for enrichment computation in pyTEenrich

One bed file per phase

```R
peaks_x_promoters %>% colnames()
```

```R
peaks_x_promoters$promoter_ensembl %>% unique() %>% length()
```

```R
peaks_x_promoters %>% dplyr::select(1:4) %>% unique() %>% nrow()
```

```R
dir.create('../out/promoters_of_rhythmic_genes_per_phase/')
for (p in peaks_x_promoters$promoter_phase_assigned %>% unique()) {
    print(p)
    peaks_x_promoters %>% 
    dplyr::filter(promoter_phase_assigned == p, promoter_padj < 0.05) %>% 
    dplyr::select(promoter_chr, promoter_start, promoter_end, promoter_ensembl) %>%
    unique() %>%
    dplyr::arrange(promoter_chr, promoter_start, promoter_end) %>%
    write.table(paste0('../out/promoters_of_rhythmic_genes_per_phase/promoters_', gsub("/", "_", p), '.bed'), sep = '\t', col.names = F, row.names = F, quote = F)
}
```

## Intersections between TFs and promoters of rhythmic genes


## Binding profile of TFs w.r.t. rhythmic genes

For each TF, we compute:
- Total number of peaks
- Total number of peaks in promoters
- Total number of peaks in rhythmic promoters
- Total number of peaks in signif. promoters peaking in each phase

```R
peaks_x_promoters %>% head()
```

Note: peaks_total is simply the number of lines in the bedfile reporting the peak coordinates for each TF.

```R
# separated total peaks per TF AND publication
total_peaks_per_TF = peaks_x_promoters %>% 
    dplyr::select(TF_symbol, filename, cell_line, publication, peaks_total) %>% 
    unique() %>%
    dplyr::rename(TF_peaks_publication = peaks_total)
```

```R
total_peaks_per_TF %>% head()
```

```R
total_peaks_per_TF %>% arrange(TF_symbol) %>% head()
```

```R
# total peaks for each TF, summing over publications
total_peaks_per_TF = total_peaks_per_TF %>% dplyr::group_by(TF_symbol) %>% dplyr::mutate(TF_peaks_total = sum(TF_peaks_publication))
```

```R
total_peaks_per_TF %>% arrange(TF_symbol) %>% head()
```

Number of ChIP-seq experiments considered

```R
total_peaks_per_TF %>% dim()
```

Number of TFs considered

```R
total_peaks_per_TF$TF_symbol %>% unique() %>% length()
```

```R
total_peaks_per_TF %>% write.table('../out/tables/total_peaks_per_TF_and_publication.tsv', row.names = F, col.names = T, sep = '\t', quote = F)
```

```R
peaks_x_promoters %>% dplyr::pull(promoter_phase_assigned) %>% unique()
```

```R
peaks_x_promoters %>% dplyr::filter(is.na(promoter_phase_assigned)) %>% nrow()
```

```R
# filtering out peaks falling within promoters of genes that are not in the rythmicity dataset
peaks_x_promoters = peaks_x_promoters %>% dplyr::filter(!is.na(promoter_phase_assigned))
```

```R
peaks_x_promoters %>% dplyr::pull(promoter_ensembl) %>% unique() %>% length()
```

```R
peaks_x_promoters %>% dplyr::pull(promoter_symbol) %>% unique() %>% length()
```

```R
rhythm_all_genes %>% dplyr::pull(ensembl) %>% unique() %>% length()
```

```R
ChIP_centric_metrics = peaks_x_promoters %>% 
    dplyr::group_by(TF_symbol) %>% 
    dplyr::summarize( 
        TF_peaks_in_promoters = n(), 
        TF_promoters_bound = length(unique(promoter_ensembl)),
        TF_peaks_in_rhythmic_promoters = sum(promoter_padj < 0.05, na.rm = T),
        TF_rhythmic_promoters_bound = length(unique(promoter_ensembl)),
        TF_peaks_in_eG1 = sum((promoter_padj < 0.05) & (promoter_phase_assigned == "eG1"), na.rm = T),
        TF_peaks_in_lG1 = sum((promoter_padj < 0.05) & (promoter_phase_assigned == "lG1"), na.rm = T),
        TF_peaks_in_G1S = sum((promoter_padj < 0.05) & (promoter_phase_assigned == "G1/S"), na.rm = T),
        TF_peaks_in_S1 = sum((promoter_padj < 0.05) & (promoter_phase_assigned == "S1"), na.rm = T),
        TF_peaks_in_S2 = sum((promoter_padj < 0.05) & (promoter_phase_assigned == "S2"), na.rm = T),
        TF_peaks_in_SG2 = sum((promoter_padj < 0.05) & (promoter_phase_assigned == "S/G2"), na.rm = T),
        TF_peaks_in_G2M = sum((promoter_padj < 0.05) & (promoter_phase_assigned == "G2/M"), na.rm = T),
        TF_peaks_in_M = sum((promoter_padj < 0.05) & (promoter_phase_assigned == "M"), na.rm = T))
```

```R
ChIP_centric_metrics %>% dplyr::filter(TF_symbol == "ZNF519")
```

```R
ChIP_centric_metrics %>% dim()
```

```R
# merging with total peaks computed separately
ChIP_centric_metrics = ChIP_centric_metrics %>% dplyr::left_join(., (total_peaks_per_TF %>% dplyr::select(TF_symbol, TF_peaks_total) %>% unique()))
```

```R
ChIP_centric_metrics %>% dim()
```

```R
ChIP_centric_metrics %>% head()
```

```R
ChIP_centric_metrics %>% write.table('../out/tables/ChIP_centric_metrics_promoters.tsv', row.names = F, col.names = T, sep = '\t', quote = F)
```

```R
gene_centric_metrics = peaks_x_promoters %>% dplyr::select(promoter_padj, TF_symbol, promoter_ensembl) %>% unique()
```

```R
gene_centric_metrics %>% head()
```

Enrichment tests: are any TFs enriched in binding to rhythmic vs. non-rhythmic promoters? 


| | genes with prom. binding by TF | genes without prom. binding by TF | total |
|---|:---:|:---:|:---:|
| rhythmic genes | x | m-x | m |
| arhythmic genes | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
TF_vs_gene_metrics = gene_centric_metrics %>% dplyr::group_by(TF_symbol) %>% dplyr::summarize(rhythmic_genes_bound = sum(promoter_padj < 0.05, na.rm = T),
                                                                        arhythmic_genes_bound = sum(promoter_padj > 0.05, na.rm = T))
TF_vs_gene_metrics %>% head()
```

```R
# computing constants
n_rhythmic_genes = rhythm_all_genes %>% dplyr::filter(padj < 0.05) %>% nrow()
total_genes = rhythm_all_genes %>% nrow()
ratio_expected = n_rhythmic_genes/total_genes
ratio_expected
```

```R
enrichment_test <- function(row) {
    x_observed = as.integer(row[["rhythmic_genes_bound"]])
    m = n_rhythmic_genes
    n = total_genes - n_rhythmic_genes
    k = as.integer(row[["rhythmic_genes_bound"]]) + as.integer(row[["arhythmic_genes_bound"]])
    
    x = 0:m # is the variable tested
    probs <- dhyper(x, m, n, k, log = FALSE)
    
    
    # we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
    pval_one_sided = sum(probs[x>=x_observed])

    # we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
    pval_two_sided = sum(probs[probs <= pval_one_sided])
    
    return(pval_two_sided)
    }
```

```R
TF_vs_gene_metrics %>% dplyr::filter(TF_symbol == "E2F1")
```

```R
TF_vs_gene_metrics %>% head()
```

```R
# adding Olga's KZFPs just to make sure we don't miss any, beyond the ones expressed in the K562 data.
```

```R
kzfps_olga = read.table('../data/from_Olga/Table4Cyril.txt', sep = '\t', header = T)
kzfps_olga$summaryKZFP %>% table()
kzfps_olga %>% head(1)
```

```R
olga_kzfps_symbols = (kzfps_olga %>% dplyr::filter(Common.Name == "Human", summaryKZFP != "No.KZFP"))$KRAB.Gene.id %>% unique() %>% as.character()
olga_kzfps_symbols %>% head()
```

```R
rhythmic_vs_arhythmic_enrichment = TF_vs_gene_metrics %>% apply(., 1, enrichment_test) 
names(rhythmic_vs_arhythmic_enrichment) = TF_vs_gene_metrics$TF_symbol
rhythmic_vs_arhythmic_enrichment = as.data.frame(rhythmic_vs_arhythmic_enrichment) %>% tibble::rownames_to_column("TF_symbol")
colnames(rhythmic_vs_arhythmic_enrichment)[2] = 'p_enrich'
rhythmic_vs_arhythmic_enrichment$is_kzfp = ifelse(rhythmic_vs_arhythmic_enrichment$TF_symbol %in% c(kzfps_cycling_info$symbol, olga_kzfps_symbols), T, F)
rhythmic_vs_arhythmic_enrichment$padj_enrich = p.adjust(rhythmic_vs_arhythmic_enrichment$p_enrich, method = "BH")
rhythmic_vs_arhythmic_enrichment$ratio_observed = TF_vs_gene_metrics$rhythmic_genes_bound/(TF_vs_gene_metrics$arhythmic_genes_bound+TF_vs_gene_metrics$rhythmic_genes_bound)
rhythmic_vs_arhythmic_enrichment$ratio_expected = ratio_expected
rhythmic_vs_arhythmic_enrichment$ratio_foldChange = rhythmic_vs_arhythmic_enrichment$ratio_observed/rhythmic_vs_arhythmic_enrichment$ratio_expected

rhythmic_vs_arhythmic_enrichment = rhythmic_vs_arhythmic_enrichment %>% 
    dplyr::left_join(., rhythm_all_genes %>% dplyr::select(symbol, padj, phase_assigned, stars, RepliTiming) %>% dplyr::rename(TF_symbol = symbol, padj_rhythm = padj, stars_rhythm = stars))
```

```R
rhythmic_vs_arhythmic_enrichment %>% head()
```

```R
rhythmic_vs_arhythmic_enrichment %>% arrange(padj_enrich) %>% head()
```

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(padj_enrich < 0.05) %>% dplyr::arrange(desc(ratio_foldChange)) %>% head(20)
```

```R
rhythmic_vs_arhythmic_enrichment %>% write.table('../out/tables/rhythmic_vs_arythmic_binding_enrichment_ENCODE_tronolabKZFPs_hughesKZFPs.tsv', col.names = T, row.names = F, sep = '\t')
```

How many TFs, how many KZFPs, how many from ENCODE etc?

```R
peaks_x_promoters %>% dplyr::select(TF_symbol) %>% unique() %>% dim()
```

```R
peaks_x_promoters %>% dplyr::filter(TF_symbol %in% c(kzfps_cycling_info$symbol, olga_kzfps_symbols)) %>% dplyr::select(TF_symbol) %>% unique() %>% dim()
```

```R
peaks_x_promoters %>% dplyr::filter(!TF_symbol %in% c(kzfps_cycling_info$symbol, olga_kzfps_symbols)) %>% dplyr::select(TF_symbol) %>% unique() %>% dim()
```

Number of rhythmic promoters bound by KZFPs:

```R
peaks_x_promoters %>% colnames()
```

```R
peaks_x_promoters %>% 
    dplyr::filter(promoter_padj < 0.05, 
                  TF_symbol %in% c(kzfps_cycling_info$symbol, olga_kzfps_symbols)) %>%
    dplyr::select(promoter_ensembl) %>%
    unique() %>%
    dim()
```

```R
peaks_x_promoters %>% 
    dplyr::filter(promoter_padj < 0.05, 
                  TF_symbol %in% c(kzfps_cycling_info$symbol, olga_kzfps_symbols)) %>%
    dplyr::select(promoter_symbol) %>%
    unique() %>%
    dim()
```

All TFs

```R
peaks_x_promoters %>% 
    dplyr::filter(promoter_padj < 0.05) %>%
    dplyr::select(promoter_ensembl) %>%
    unique() %>%
    dim()
```

```R
peaks_x_promoters %>% 
    dplyr::filter(promoter_padj < 0.05) %>%
    dplyr::select(promoter_symbol) %>%
    unique() %>%
    dim()
```

Are genes coding for rhythmic TFs generally more enriched at promoters of rhythmic genes?

```R
rhythmic_vs_arhythmic_enrichment %>% ggplot(aes(x = padj_rhythm < 0.05, y = -log10(padj_enrich))) + geom_boxplot()
```

Not striking enrichment of rhythmic TFs amongst rhythmic promoter binders, though they seem to have a slightly higher mean


```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(!is_kzfp) %>% ggplot(aes(x = padj_rhythm < 0.05, y = -log10(padj_enrich))) + geom_boxplot()
```

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(!is_kzfp) %>% ggplot(aes(x = padj_rhythm < 0.05)) + geom_bar()
```

For non KZFP genes, the enrichment at rhythmic promoters is the same whether the gene is rhythmic or not.

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(is_kzfp) %>% ggplot(aes(x = padj_rhythm < 0.05)) + geom_bar()
```

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(is_kzfp) %>% ggplot(aes(x = padj_rhythm < 0.05, y = -log10(padj_enrich))) + geom_boxplot()
```

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(is_kzfp, padj_enrich < 0.05) %>% arrange(p_enrich) 
```

```R
rhythmic_vs_arhythmic_enrichment %>% 
    dplyr::filter(is_kzfp) %>% 
    arrange(p_enrich) %>% 
    ggplot(aes(x = ratio_foldChange, y = -log10(p_enrich), col = padj_enrich < 0.05)) +
    geom_point()
```

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(is.na(phase_assigned))
```

A lot of the TFs for which we have binding data, but no expression data in Romain's bulk phase-sorted RNAseq are ZNFs. Probaly due to low expression again...

```R
rhythmic_vs_arhythmic_enrichment %>% arrange(padj_enrich)
```

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(grepl("E2F", TF_symbol)) %>% arrange(padj_enrich)
```

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(grepl("MCM", TF_symbol)) %>% arrange(padj_enrich)
```

Hypothesis TFs which are highly rhythmic, yet do not enrich at rhythmic promoters are likely involved in replication, or in general transcription. They are DNA binding, but do not control transcription, rather replication. In contrast, TFs that are highly rhythmic and enrich at the promoters of rhythmic genes are likely TFs controlling cell cycling transcriptionally, that is pushing the cells towards via gene products.

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(padj_rhythm < 0.05) %>% arrange(desc(padj_enrich), padj_rhythm) %>% head(40)
```

Lots of KZFPs near the end of this list... is this to suggest that rhythmic KZFPs do not promote gene expression (most, at least) but rather control replication in some whay? TODO: cross this score with the perturb-seq data, and possibly a DNA damage signature, or replication stress signature.

```R
rhythmic_vs_arhythmic_enrichment %>% dplyr::filter(!is_kzfp, padj_rhythm < 0.05) %>% arrange(padj_enrich) %>% head()
```

```R
(rhythmic_vs_arhythmic_enrichment$padj_enrich < 0.05) %>% sum()
```

```R
rhythmic_vs_arhythmic_enrichment %>% nrow()
```

Most of the assayed TFs enrich within rhythmic vs non rhythmic promoters. But this may be because we're looking at highly transcribed genes as well.

```R
rhythmic_vs_arhythmic_enrichment %>% ggplot(aes(y = -log10(padj_enrich), x = is_kzfp)) + geom_boxplot()
```

Compared to other TFs with available ChIP-seq data, KZFPs are generally excluded from the promoter of rhythmic genes. But there are still exceptions with elevated enrichment.

The trend holds when excluding the Tronolab chip-seq


## Enrichment of TF binding to promoters of rhythmic genes within each phase of the cell cycle.
- Each phase vs. all other promoters

- Each phase vs. all rhythmic promoters


```R
gene_centric_metrics %>% head()
```

```R
# adding the phase information for the genes
phase_centric_metrics = gene_centric_metrics %>% dplyr::left_join(., rhythm_all_genes %>% dplyr::select(ensembl, symbol, phase_assigned), by = c("promoter_ensembl" = "ensembl"))
phase_centric_metrics %>% head()
```

```R
phase_centric_metrics %>% dim()
```

```R
phase_centric_metrics %>% dplyr::select(promoter_padj, TF_symbol, promoter_ensembl) %>% unique() %>% dim()
```

```R
phase_centric_metrics %>% dplyr::select(promoter_padj, TF_symbol, symbol) %>% unique() %>% dim()
```

No problem of gene duplication here.


| | genes with prom. binding by TF | genes without prom. binding by TF | total |
|---|:---:|:---:|:---:|
| Genes in phase | x | m-x | m |
| Genes not in phase | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
total_genes = rhythm_all_genes %>% nrow()
total_genes
```

```R
# result list, to be merged afterwards
phase_selected = "eG1"
phase_list = rhythm_all_genes$phase_assigned %>% unique()
df_list = list()
```

```R
enrichment_test_phase_vs_all_genes <- function(row) {
    x_observed = as.integer(row[["rhythmic_genes_in_phase_bound"]])
    m = n_genes_in_phase
    n = total_genes - n_genes_in_phase
    k = as.integer(row[["rhythmic_genes_in_phase_bound"]]) + as.integer(row[["genes_not_in_phase_bound"]])
    
    x = 0:m # is the variable tested
    probs <- dhyper(x, m, n, k, log = FALSE)
    
    
    # we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
    pval_one_sided = sum(probs[x>=x_observed])

    # we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
    pval_two_sided = sum(probs[probs <= pval_one_sided])
    
    return(pval_two_sided)
}
```

```R
rhythm_all_genes %>% head()
rhythm_all_genes %>% dim()
rhythm_all_genes %>% colnames()
```

```R
i = 1
for (phase_selected in rhythm_all_genes$phase_assigned %>% unique()) {
    print(phase_selected)
    n_genes_in_phase = rhythm_all_genes %>% 
        dplyr::filter(phase_assigned == phase_selected, padj < 0.05) %>% 
        nrow()
    TF_vs_genephase_metrics = phase_centric_metrics %>% 
        dplyr::group_by(TF_symbol) %>% 
        dplyr::summarize(rhythmic_genes_in_phase_bound = sum((promoter_padj < 0.05) & (phase_assigned == phase_selected), na.rm = T),
                                                                        genes_not_in_phase_bound = sum(!((promoter_padj < 0.05) & (phase_assigned == phase_selected)),na.rm = T))
    
    df_list[[phase_selected]] = TF_vs_genephase_metrics %>% 
        apply(., 1, enrichment_test_phase_vs_all_genes)# %>% 
       #p.adjust(., "BH", n = 8*nrow(TF_vs_genephase_metrics)) # we opt to correct for multiple testing later
    i = i+1
    }
```

```R
df_list %>% names()
```

```R
in_phase_vs_all_genes_enrichment = df_list %>% as.data.frame()
in_phase_vs_all_genes_enrichment %>% dim()
in_phase_vs_all_genes_enrichment %>% head()
```

```R
in_phase_vs_all_genes_enrichment$TF_symbol = TF_vs_genephase_metrics$TF_symbol
in_phase_vs_all_genes_enrichment$is_kzfp = ifelse(in_phase_vs_all_genes_enrichment$TF_symbol %in% c(kzfps_cycling_info$symbol, olga_kzfps_symbols), T, F)
in_phase_vs_all_genes_enrichment = in_phase_vs_all_genes_enrichment %>%
    dplyr::left_join(., rhythm_all_genes %>% 
                     dplyr::select(symbol, padj, phase_assigned, stars, RepliTiming) %>% 
                     dplyr::rename(TF_symbol = symbol, padj_rhythm = padj, stars_rhythm = stars))
```

```R
in_phase_vs_all_genes_enrichment %>% head()
```

```R
in_phase_vs_all_genes_enrichment %>% dplyr::select(TF_symbol) %>% unique() %>% nrow()
in_phase_vs_all_genes_enrichment %>% dplyr::filter(is_kzfp) %>% dplyr::select(TF_symbol) %>% unique() %>% nrow()
in_phase_vs_all_genes_enrichment %>% dplyr::filter(!is_kzfp) %>% dplyr::select(TF_symbol) %>% unique() %>% nrow()
```

```R
phase_selected = "G1.S"
in_phase_vs_all_genes_enrichment %>% arrange(!!as.symbol(phase_selected)) %>% head(50) %>% dplyr::select(all_of(phase_selected), TF_symbol, stars_rhythm, phase_assigned)
```

```R
in_phase_vs_all_genes_enrichment %>% arrange(!!as.symbol(phase_selected)) %>% dplyr::filter(TF_symbol == "E2F1") %>% dplyr::select(all_of(phase_selected), TF_symbol, stars_rhythm, phase_assigned)
```

ZNF519 is top 10, higher than E2F1

```R
in_phase_vs_all_genes_enrichment %>% dplyr::filter(TF_symbol == "ZNF519")
```

What if we compute an average ranking on the p_enrich values, for rhythmic KZFPs and TFs?

```R
n_TFs = in_phase_vs_all_genes_enrichment %>% dplyr::select(TF_symbol) %>% unique() %>% nrow()
n_TFs
```

```R
in_phase_vs_all_genes_enrichment$median_rank_signif_pval =  in_phase_vs_all_genes_enrichment %>% 
    dplyr::select(1:8) %>% 
    apply(., 2, function(x) ifelse(x > 0.05, yes = 1, no = x)) %>%
    apply(., 2, rank, ties="average") %>% 
    apply(., 1, median)

in_phase_vs_all_genes_enrichment$best_rank_signif_pval =  in_phase_vs_all_genes_enrichment %>% 
    dplyr::select(1:8) %>% 
    apply(., 2, function(x) ifelse(x > 0.05, yes = 1, no = x)) %>%
    apply(., 2, rank, ties="average") %>% 
    apply(., 1, min)
          
          
in_phase_vs_all_genes_enrichment$median_rank_signif_padj =  in_phase_vs_all_genes_enrichment %>% 
    dplyr::select(1:8) %>% 
    apply(., 2, function(x) p.adjust(x, method = "BH", n = 8*n_TFs)) %>%
    apply(., 2, function(x) ifelse(x > 0.05, yes = 1, no = x)) %>%
    apply(., 2, rank, ties="average") %>% 
    apply(., 1, median)

in_phase_vs_all_genes_enrichment$best_rank_signif_padj =  in_phase_vs_all_genes_enrichment %>% 
    dplyr::select(1:8) %>% 
    apply(., 2, function(x) p.adjust(x, method = "BH", n = 8*n_TFs)) %>%
    apply(., 2, function(x) ifelse(x > 0.05, yes = 1, no = x)) %>%
    apply(., 2, rank, ties="average") %>% 
    apply(., 1, min)
```

```R
in_phase_vs_all_genes_enrichment %>% arrange(median_rank_signif_pval) %>% head(10)
```

```R
in_phase_vs_all_genes_enrichment %>% arrange(best_rank_signif_pval) %>% head(10)
```

This yields very similar results to doing the enrichment on all rhythmic promoters at once. The genes coming up are those involved in transcription, in general. By contrast, focusing on each phase separately yields phase-specific TFs, such as the E2F TFs.

```R
in_phase_vs_all_genes_enrichment %>% dplyr::filter(padj_rhythm < 0.05, is_kzfp) %>% arrange(median_rank_signif_pval) %>% head(10)
```

```R
in_phase_vs_all_genes_enrichment %>% dplyr::filter(padj_rhythm < 0.05, is_kzfp) %>% arrange(median_rank_signif_padj) %>% head(10)
```

```R
in_phase_vs_all_genes_enrichment %>% dplyr::filter(padj_rhythm < 0.05, is_kzfp) %>% arrange(best_rank_signif_pval) %>% head(10)
```

```R
in_phase_vs_all_genes_enrichment %>% dplyr::filter(padj_rhythm < 0.05, is_kzfp) %>% arrange(best_rank_signif_padj) %>% head(10)
```

```R
in_phase_vs_all_genes_enrichment %>% 
    dplyr::filter(is_kzfp) %>%
    dplyr::filter_at(vars(1:8), any_vars(.<0.05)) %>%
    arrange(best_rank_signif_pval)
```

```R
in_phase_vs_all_genes_enrichment %>% 
    dplyr::filter(!is_kzfp) %>%
    dplyr::filter_at(vars(1:8), any_vars(.<0.05)) %>%
    arrange(best_rank_signif_pval)
```

```R
in_phase_vs_all_genes_enrichment %>% dplyr::select(TF_symbol, is_kzfp) %>% unique() %>% dplyr::pull(is_kzfp) %>% table()
```

```R
in_phase_vs_all_genes_enrichment$pval_enrich_min = apply(in_phase_vs_all_genes_enrichment %>% dplyr::select(1:8), 1, min, simplify = T)
```

## Plotting the enrichment over rhythmic promoters results

```R
# non-adjusted p
p = in_phase_vs_all_genes_enrichment %>% 
    dplyr::select(1:8) %>% 
    apply(., 2, function(x) sum(x< 0.05)) %>% 
              as.data.frame() %>% 
          tibble::rownames_to_column("phase") %>%
          dplyr::mutate(phase = factor(gsub("\\.", "/", phase), 
                                       levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M"))) %>%
        dplyr::rename(n_TFs = ".") %>% 
            ggplot(aes(x = phase, y = n_TFs)) + 
          geom_bar(stat = "identity") + 
            theme_classic() + 
            ylab("# TFs") + 
          xlab("") + 
          theme(axis.line.x = element_blank(),
                axis.ticks.x = element_blank())
p
```

```R
ggsave("../out/cell_cycle_figures/n_tfs_enriched_in_rhythmic_proms_per_phase_pval.svg", p, svg, height = 1.5, width = 3)
```

```R
# adjusted p
p = in_phase_vs_all_genes_enrichment %>% 
    dplyr::select(1:8) %>% 
    apply(., 2, function(x) p.adjust(x, "BH", 8*n_TFs)) %>% 
    apply(., 2, function(x) sum(x< 0.05)) %>% 
              as.data.frame() %>% 
          tibble::rownames_to_column("phase") %>%
          dplyr::mutate(phase = factor(gsub("\\.", "/", phase), 
                                       levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M"))) %>%
        dplyr::rename(n_TFs = ".") %>% 
            ggplot(aes(x = phase, y = n_TFs)) + 
          geom_bar(stat = "identity") + 
            theme_classic() + 
            ylab("# TFs") + 
          xlab("") + 
          theme(axis.line.x = element_blank(),
                axis.ticks.x = element_blank())
p
```

Over-representation of TFs enriching at the promoters of S2-to-G2/M, but that's expected since that's when the most expressed genes peak.

```R
ggsave("../out/cell_cycle_figures/n_tfs_enriched_in_rhythmic_proms_per_phase_padj.svg", p, svg, height = 1.5, width = 3)
```

```R
in_phase_vs_all_genes_enrichment %>% dplyr::select(1:8) %>% apply(., 2, function(x) sum(x< 0.05))
in_phase_vs_all_genes_enrichment %>% dplyr::select(1:8) %>% apply(., 1, function(x) sum(x< 0.05) > 0) %>% sum()
```

```R
in_phase_vs_all_genes_enrichment %>% dplyr::select(1:8) %>% apply(., 2, function(x) p.adjust(x, "BH", 8*n_TFs)) %>% apply(., 2, function(x) sum(x< 0.05))
in_phase_vs_all_genes_enrichment %>% dplyr::select(1:8) %>% apply(., 2, function(x) p.adjust(x, "BH", 8*n_TFs)) %>% apply(., 1, function(x) sum(x< 0.05) > 0) %>% sum()
```

We find 380 out of 744 TFs that are enriched in at least one phase. What if we enrich over them?

```R
in_phase_vs_all_genes_enrichment %>% head()
```

```R
in_phase_vs_all_genes_enrichment_long = in_phase_vs_all_genes_enrichment %>% 
    tidyr::pivot_longer(cols = 1:8, names_to = "promoters_active_in_phase", values_to = "pval_enrich") %>% 
    dplyr::mutate(promoters_active_in_phase = gsub("\\.", "/", promoters_active_in_phase))
in_phase_vs_all_genes_enrichment_long = in_phase_vs_all_genes_enrichment_long %>% dplyr::mutate(padj_enrich = p.adjust(pval_enrich, "BH", 8*n_TFs))
in_phase_vs_all_genes_enrichment_long %>% head()
```

Let's try just with the adjusted p values for the enrichment plots

```R
enriched_once = in_phase_vs_all_genes_enrichment %>% 
    dplyr::select(1:8) %>% 
    apply(., 2, function(x) p.adjust(x, "BH", 8*n_TFs)) %>%
    apply(., 1, function(x) sum(x<0.05)>0)
```

```R
in_phase_vs_all_genes_enrichment_long = in_phase_vs_all_genes_enrichment_long %>% 
    dplyr::mutate(promoters_active_in_phase = factor(promoters_active_in_phase, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")))
```

## Plotting example genes and ZNF519

Ranking based on signif. padj enrich

```R
in_phase_vs_all_genes_enrichment_long$rank = in_phase_vs_all_genes_enrichment_long %>% 
    dplyr::mutate(padj_enrich = ifelse(padj_enrich > 0.05, yes = 1, no = padj_enrich)) %>%
    group_by(promoters_active_in_phase) %>%
    dplyr::mutate(rank = rank(padj_enrich, ties = "average")) %>% 
    dplyr::pull(rank)
```

Focusing on rhythmic TFs does yield more specialized TFs, otherwise we are flooded by things such as general TFs, RNA pol II subunits etc.

There are TFs which score very high in one phase, but not in the others. These appear to be "specialized", as in they enrich particularly high in a specific phase. We could find them by dividing their median_rank by their best rank. The top ones should be specialized.

```R
in_phase_vs_all_genes_enrichment_long = in_phase_vs_all_genes_enrichment_long %>% dplyr::mutate(specificity = median_rank_signif_padj / rank)
```

```R
in_phase_vs_all_genes_enrichment_long %>% 
    dplyr::group_by(promoters_active_in_phase) %>% 
    dplyr::arrange(desc(specificity)) %>% pull(specificity) %>% hist()
```

Saving the data for external plotting:


```R
n_promoters_bound_per_TF = peaks_x_promoters %>% dplyr::select(TF_symbol, promoter_ensembl) %>%unique() %>% dplyr::group_by(TF_symbol) %>% summarize(n_promoters_bound = n())
n_promoters_bound_per_TF %>% head()
```

```R
peaks_x_promoters %>% colnames()
```

```R
in_phase_vs_all_genes_enrichment_long = in_phase_vs_all_genes_enrichment_long %>% dplyr::left_join(., n_promoters_bound_per_TF)
```

```R
in_phase_vs_all_genes_enrichment_long %>% write.table('../out/tables/TF_binding_enrichment_in_rhythmic_genes_from_each_phase.tsv',
                                                     sep = '\t',
                                                     col.names = T,
                                                     row.names = F,
                                                     quote = F)
```

```R
in_phase_vs_all_genes_enrichment_long %>% head()
```

```R
phase_centric_metrics %>% dplyr::filter(TF_symbol == "ZNF682") %>% dplyr::select(symbol) %>% write.table('../out/temp/znf682_bound.tsv', sep = '\t', quote = F, col.names = F, row.names = F)
```

Shows an enrichment for ZNC-C2H2 containing proteins.


## Enrichment of TF binding at promoters and ZNFs of M-G1 vs G1/S-G2 KZFPs 

Before running this part, the script ``

Reguation of KZFP expression may take place at the TSS. But it may also take place at their 3' ZNF-encoding region, which is rich in epigenomic signal (H3K9me3, ZNF274, ATRX, H3K36me3 etc).

We'll aggregate the binding enrichments computed over these two regions and adjust the p-values for multiple testing. Will ZNF274 binding on ZNFs be the top signal?


| | KZFP genes with prom. binding by TF | KZFP genes without prom. binding by TF | total |
|---|:---:|:---:|:---:|
| KZFP genes in M-lG1 | x | m-x | m |
| KZFP genes not in M-lG1 | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |


### Loading the ENCODE x ZNF bedtools intersect outputs

```R
# adding the peaks of all TFs that overlap with KZFP ZNFs

i = 0
df_list_encode_znfs = list()
for (bed in list.files('../out/encode_x_kzfp_znfs/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    p = paste0('../out/encode_x_kzfp_znfs/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        temp_df = read.table(paste0('../out/encode_x_kzfp_znfs/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        temp_df$TF_symbol = encode_tfs_metadata[grepl(filename, encode_tfs_metadata$Files) %>% which(), "Target.gene.symbol"]
        temp_df$peaks_total = R.utils::countLines(paste0('../data/ENCODE_k562_tf_peaks/bed_filtered_hg19/', bed))
        i = i+1

        df_list_encode_znfs[[i]] = temp_df
        }
}
```

```R
peaks_x_kzfp_znfs_encode = do.call(rbind, df_list_encode_znfs) %>% dplyr::filter(V7 != ".")
peaks_x_kzfp_znfs_encode$cell_line = "K562"
peaks_x_kzfp_znfs_encode$publication = "ENCODE"
peaks_x_kzfp_znfs_encode %>% head() %>% print()
```

```R
min(23689462 , 23689419   ) - max(23688470, 23688783)
#V17
```

```R
peaks_x_kzfp_znfs_encode = peaks_x_kzfp_znfs_encode %>% dplyr::select(-V5, -V10, -V11, -V12, -V14, -V15, -V16)
peaks_x_kzfp_znfs_encode %>% dim()
```

```R
colnames(peaks_x_kzfp_znfs_encode) = c("znf_chr", "znf_start", "znf_end", "znf_symbol", "znf_strand", 
                                      "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "filename", "TF_symbol", "peaks_total", "cell_line", "publication")
peaks_x_kzfp_znfs_encode %>% head()
peaks_x_kzfp_znfs_encode %>% dim()
```

```R
peaks_x_kzfp_znfs_encode %>% head()
```

```R
peaks_x_kzfp_znfs_encode %>% dplyr::group_by(TF_symbol) %>% dplyr::summarize(n = n()) %>% arrange(desc(n)) %>% head(10)
```

This is looking really good: TRIM28, ZNF274 and SETDB1 are the three top TFs with absolute number of peaks at ZNFs. Now, whether that constitutes enrichment, and if there is enrichment between M-lG1 and S2-G2 is another question. CBX3 is also involved in binding heterochromatin. RBM25 binds mRNA ?? for mRNA splicing??


### Loading the Tronolab x ZNF bedtools intersect outputs

```R
peaks_x_kzfp_znfs_tronolab = read.table('../out/kzfp_znfs_vs_kzfp_peaks.bed', sep = '\t', header = F)
```

```R
peaks_x_kzfp_znfs_tronolab %>% head()
```

```R
# identifying the overlap column:
min(3381507, 3380958) - max(3380353, 3380319)
# V13
```

```R
peaks_x_kzfp_znfs_tronolab %>% dim()
```

```R
peaks_x_kzfp_znfs_tronolab = peaks_x_kzfp_znfs_tronolab %>% dplyr::select(-V5, -V12)
```

```R
colnames(peaks_x_kzfp_znfs_tronolab) = c("znf_chr", "znf_start", "znf_end", "znf_symbol", "znf_strand", "peak_chr", "peak_start", "peak_end", "peak_ID", "peak_score", "peak_overlap")
```

```R
peaks_x_kzfp_znfs_tronolab %>% head()
```

```R
any(!grepl("sampled_peak", peaks_x_kzfp_znfs_tronolab$peak_ID))
```

```R
# adding the peak information
peaks_x_kzfp_znfs_tronolab = peaks_x_kzfp_znfs_tronolab %>% dplyr::left_join(., peaks %>% dplyr::select(peak_ID, TF_symbol, peak_lowconf, TF_peaks_total))
```

```R
# recovering the filename (i.e name of the chip itself)
peaks_x_kzfp_znfs_tronolab$filename = sapply(peaks_x_kzfp_znfs_tronolab$peak_ID, function(x) str_split(x, "_sampled_peak_")[[1]][[1]])
```

```R
peaks_x_kzfp_znfs_tronolab %>% head()
```

### Loading the Schmitges x ZNF bedtools intersect outputs

```R
i = 0
df_list = list()
for (bed in list.files('../out/schmitges_x_kzfp_znfs/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    p = paste0('../out/schmitges_x_kzfp_znfs/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        temp_df = read.table(paste0('../out/schmitges_x_kzfp_znfs/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        #temp_df$TF_symbol = schmitges_metadata %>% can be done afterwards
        temp_df$peaks_total = R.utils::countLines(paste0('../data/schmitges_kzfp_chips/hg19_peaks_chipatlas/', bed))
        i = i+1

        df_list[[i]] = temp_df}
    }
```

```R
peaks_x_kzfp_znfs_schmitges = do.call(rbind, df_list)
```

```R
peaks_x_kzfp_znfs_schmitges %>% head() %>% print()
```

```R
min(133683167 , 133683087) - max(133682343, 133682947)
#V17
```

```R
peaks_x_kzfp_znfs_schmitges$cell_line = "HEK293T"
peaks_x_kzfp_znfs_schmitges$publication = "schmitges"
peaks_x_kzfp_znfs_schmitges = peaks_x_kzfp_znfs_schmitges %>% dplyr::left_join(., schmitges_metadata %>% dplyr::select(filename, TF_symbol))
peaks_x_kzfp_znfs_schmitges %>% head() %>% print()
```

```R
peaks_x_kzfp_znfs_schmitges = peaks_x_kzfp_znfs_schmitges %>% dplyr::select(-V5, -V10, -V12, -V13, -V14, -V15, -V16)
peaks_x_kzfp_znfs_schmitges %>% head()
```

```R
colnames(peaks_x_kzfp_znfs_schmitges) = c("znf_chr", "znf_start", "znf_end", "znf_symbol", "znf_strand", 
                                      "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "filename", "peaks_total", "cell_line", "publication", "TF_symbol")
peaks_x_kzfp_znfs_schmitges %>% head()
peaks_x_kzfp_znfs_schmitges %>% dim()
```

### Loading the Najafabadi x ZNF bedtools intersect outputs

```R
# adding the najafabadi peaks at znfs
i = 0
df_list = list()
for (bed in list.files('../out/najafabadi_all_znfs_x_all_znfs/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    p = paste0('../out/najafabadi_all_znfs_x_all_znfs/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        temp_df = read.table(paste0('../out/najafabadi_all_znfs_x_all_znfs/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        temp_df$sample = str_split(filename, "_")[[1]][[1]]
        temp_df$peaks_total = R.utils::countLines(paste0('../data/najafabadi_all_znfs_chips/hg19_peaks_80/', bed))
        i = i+1

        df_list[[i]] = temp_df}
    }
```

```R
peaks_x_kzfp_znfs_najafabadi = do.call(rbind, df_list)
```

```R
peaks_x_kzfp_znfs_najafabadi %>% head()
```

```R
peaks_x_kzfp_znfs_najafabadi$cell_line = "HEK293T"
peaks_x_kzfp_znfs_najafabadi$publication = "najafabadi"
peaks_x_kzfp_znfs_najafabadi %>% head() %>% print()
```

```R
najafabadi_metadata = read.table('../data/najafabadi_all_znfs_chips/najafabadi_metadata.txt', sep = ' ', header = F) %>% 
    dplyr::rename(symbol = V1, sample = V2)
najafabadi_metadata %>% head()
```

```R
peaks_x_kzfp_znfs_najafabadi = peaks_x_kzfp_znfs_najafabadi %>% 
    dplyr::left_join(., najafabadi_metadata)
```

```R
peaks_x_kzfp_znfs_najafabadi %>% head()
```

```R
peaks_x_kzfp_znfs_najafabadi = peaks_x_kzfp_znfs_najafabadi %>% dplyr::select(-V5, -V10)
```

```R
peaks_x_kzfp_znfs_najafabadi %>% head()
```

```R
colnames(peaks_x_kzfp_znfs_najafabadi) = c("znf_chr", "znf_start", "znf_end", "znf_symbol", "znf_strand", 
                                      "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "filename", "sample", "peaks_total", "cell_line", "publication", "TF_symbol")
peaks_x_kzfp_znfs_najafabadi %>% head()
peaks_x_kzfp_znfs_najafabadi %>% dim()
```

```R
# concatenating with the encode data
peaks_x_kzfp_znfs_tronolab$cell_line = "HEK293T"
peaks_x_kzfp_znfs_tronolab$publication = "Tronolab"
```

```R
peaks_x_kzfp_znfs_encode %>% head()
```

```R
peaks_x_kzfp_znfs_tronolab %>% head()
```

```R
peaks_x_kzfp_znfs = rbind(peaks_x_kzfp_znfs_encode %>% dplyr::select("znf_chr", "znf_start", "znf_end", "znf_symbol", "znf_strand", 
                                             "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "TF_symbol", "cell_line", "publication", "peaks_total", "filename"),
      peaks_x_kzfp_znfs_tronolab %>% dplyr::rename(peaks_total = TF_peaks_total) %>% dplyr::select("znf_chr", "znf_start", "znf_end", "znf_symbol", "znf_strand", 
                                             "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "TF_symbol", "cell_line", "publication", "peaks_total", "filename"),
                         peaks_x_kzfp_znfs_najafabadi %>% dplyr::select("znf_chr", "znf_start", "znf_end", "znf_symbol", "znf_strand", 
                                             "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "TF_symbol", "cell_line", "publication", "peaks_total", "filename"),
                         peaks_x_kzfp_znfs_schmitges %>% dplyr::select("znf_chr", "znf_start", "znf_end", "znf_symbol", "znf_strand", 
                                             "peak_chr", "peak_start", "peak_end", "peak_score", "peak_overlap", "TF_symbol", "cell_line", "publication", "peaks_total", "filename"))
```

```R
peaks_x_kzfp_znfs %>% dim()
peaks_x_kzfp_znfs$TF_symbol %>% unique() %>% length()
```

```R
# computing promoters and peak length
peaks_x_kzfp_znfs = peaks_x_kzfp_znfs %>% dplyr::mutate(znf_length = znf_end-znf_start, peak_length = peak_end-peak_start)
```

```R
peaks_x_kzfp_znfs %>% head()
```

```R
# adding the gene rhythmicity information, both for TFs only since all targets are KZFPs. 
rhythm_all_genes = read.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', header = 1)
rhythm_all_genes %>% head()
```

```R
# adding the TF cycling information
peaks_x_kzfp_znfs = peaks_x_kzfp_znfs %>% 
    dplyr::left_join(x=., y=rhythm_all_genes, by = c("TF_symbol" = "symbol"))
peaks_x_kzfp_znfs %>% colnames()
```

```R
colnames(peaks_x_kzfp_znfs)[((peaks_x_kzfp_znfs %>% ncol())-(rhythm_all_genes %>% ncol())+2) : (peaks_x_kzfp_znfs %>% ncol())] = paste0("TF_", colnames(peaks_x_kzfp_znfs)[((peaks_x_kzfp_znfs %>% ncol())-(rhythm_all_genes %>% ncol())+2) : (peaks_x_kzfp_znfs %>% ncol())])
```

```R
peaks_x_kzfp_znfs %>% colnames()
peaks_x_kzfp_znfs %>% head()
```

```R
# adding the KZFP information for the ZNFs level information
peaks_x_kzfp_znfs = peaks_x_kzfp_znfs %>% dplyr::left_join(., y = kzfps_cycling_info %>% 
                                                           dplyr::rename(znf_symbol = symbol, znf_cluster = cluster, znf_phase_assigned = phase_assigned, znf_padj = padj) %>% 
                                                           dplyr::select(znf_symbol, znf_cluster, znf_phase_assigned, znf_padj))
```

```R
## end of adding the schmitges and najafabadi data for computing enrichments over ZNFs
```

```R
#rif1_high_clusters = c("chr6.1", "chr19.8", "chr19.11", "chr19.9")
rif1_high_clusters = c()
phase_centric_metrics_TFs_at_znfs = peaks_x_kzfp_znfs %>% dplyr::filter(!znf_cluster %in% rif1_high_clusters) %>% dplyr::select(znf_padj, TF_symbol, znf_symbol, znf_phase_assigned) %>% unique()
phase_centric_metrics_TFs_at_znfs %>% head()
phase_centric_metrics_TFs_at_znfs %>% dim()
```

```R
phase_centric_metrics_TFs_at_znfs$znf_symbol %>% unique() %>% length()
```

```R
phase_centric_metrics_TFs_at_znfs %>% dim()
```

```R
phase_centric_metrics %>% head()
```

```R
kzfp_centric_metrics = phase_centric_metrics %>% dplyr::filter(promoter_ensembl %in% (kzfps_cycling_info %>% dplyr::filter(!cluster %in% rif1_high_clusters))$ensembl)
kzfp_centric_metrics %>% dim()
kzfp_centric_metrics %>% dplyr::select(promoter_ensembl) %>% unique() %>% dim()
```

```R
kzfp_centric_metrics %>% head()
```

```R
kzfp_tss_centric_metrics = rbind(kzfp_centric_metrics %>%
    dplyr::select(-promoter_ensembl) %>% 
    dplyr::rename(znf_padj = promoter_padj, znf_symbol = symbol, znf_phase_assigned = phase_assigned) %>%
                             dplyr::left_join(kzfps_cycling_info %>% 
                                              dplyr::select(symbol, cluster), by = c("znf_symbol" = "symbol")))
```

```R
kzfp_tss_centric_metrics %>% head()
```

```R
kzfp_znf_centric_metrics = phase_centric_metrics_TFs_at_znfs %>%
 dplyr::left_join(kzfps_cycling_info %>% 
                  dplyr::select(symbol, cluster), by = c("znf_symbol" = "symbol"))
```

```R
kzfp_znf_centric_metrics %>% tail()
```

```R
# TODO: M-G1 vs all KZFPs, M-G1 vs all rhythmic KZPFs, 

# ... G1/S-G2 vs all
```

| | RIF1-low KZFP genes with prom. binding by TF | RIF1-low KZFP genes without prom. binding by TF | total |
|---|:---:|:---:|:---:|
| RIF1-low KZFP genes in M-lG1 | x | m-x | m |
| RIF1-low KZFP genes not in M-lG1 | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
# params: rhythmicity threshold

```

```R
any(kzfps_cycling_info$cluster %in% rif1_high_clusters)
```

```R
rif1_high_clusters
```

```R
# result list, to be merged afterwards
rhythm_threshold = 0.1
phases = list(c("M", "eG1", "lG1", "G1/S", "S1"), c("S2", "S/G2", "G2"))
i = 1
# phase_list = rhythm_all_genes$phase_assigned %>% unique()
df_list = list()
n_kzfps_in_phases = kzfps_cycling_info %>% dplyr::filter(!cluster %in% rif1_high_clusters, padj <= rhythm_threshold) %>% dplyr::filter(phase_assigned %in% phases[[i]]) %>% nrow()
n_kzfps = nrow(kzfps_cycling_info %>% dplyr::filter(!cluster %in% rif1_high_clusters, padj <= rhythm_threshold))
```

```R
n_kzfps_in_phases
n_kzfps
```

```R
kzfp_tss_centric_metrics %>% head()
```

```R
TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_tss = kzfp_tss_centric_metrics %>% dplyr::group_by(TF_symbol) %>% dplyr::summarize(rhythmic_genes_in_phases_bound = sum(((znf_padj <= rhythm_threshold) & (znf_phase_assigned %in% phases[[i]]))),
                                                                        rhythmic_genes_not_in_phases_bound = sum(((znf_padj <= rhythm_threshold) & (!znf_phase_assigned %in% phases[[i]])), na.rm = T))
```

```R
TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_tss %>% head()
```

```R
TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_tss %>% dplyr::filter(TF_symbol =="ZNF274")
```

```R
enrichment_test_phases_vs_other_rhythmic_kzfps <- function(row) {
    x_observed = as.integer(row[["rhythmic_genes_in_phases_bound"]])
    m = n_kzfps_in_phases
    n = n_kzfps - n_kzfps_in_phases
    k = as.integer(row[["rhythmic_genes_in_phases_bound"]]) + as.integer(row[["rhythmic_genes_not_in_phases_bound"]])
    
    x = 0:m # is the variable tested
    probs <- dhyper(x, m, n, k, log = FALSE)
    
    # we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
    pval_one_sided = sum(probs[x>=x_observed])

    # we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
    pval_two_sided = sum(probs[probs <= pval_one_sided])
    
    return(pval_two_sided)
}
```

```R
in_phase_vs_not_in_phase_enrichment_tss = TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_tss %>% apply(., 1, enrichment_test_phases_vs_other_rhythmic_kzfps) %>% as.data.frame()
```

```R
in_phase_vs_not_in_phase_enrichment_tss$TF_symbol = TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_tss$TF_symbol
colnames(in_phase_vs_not_in_phase_enrichment_tss)[1] = 'p_enrich'
in_phase_vs_not_in_phase_enrichment_tss %>% head()
```

```R
in_phase_vs_not_in_phase_enrichment_tss$is_kzfp = ifelse(in_phase_vs_not_in_phase_enrichment_tss$TF_symbol %in% kzfps_cycling_info$symbol, T, F)
in_phase_vs_not_in_phase_enrichment_tss$padj_enrich = p.adjust(in_phase_vs_not_in_phase_enrichment_tss$p_enrich, method = "BH")
in_phase_vs_not_in_phase_enrichment_tss = in_phase_vs_not_in_phase_enrichment_tss %>% 
    dplyr::left_join(., rhythm_all_genes %>% dplyr::select(symbol, padj, phase_assigned, stars, RepliTiming) %>% dplyr::rename(TF_symbol = symbol, padj_rhythm = padj, stars_rhythm = stars))
```

```R
in_phase_vs_not_in_phase_enrichment_tss %>% arrange(p_enrich) %>% head(10)
```

The top scorers for M-lG1 rhythmic KZFP TSS are ASH1L, ZFP69B, ZNF682, but they are not significant
The top scorers for M-lG1 rhythmic RIF1-low KZFP TSS are ZNFs: ZNF251 (L1M binder) and ZNF682. They don't even pass the non-adjusted threshold.


S-G2: SMAD1, ZBTB1, E2F5 score high

Neither pass the padj threshold.

Now, binding at the ZNF:

```R
ylim_high = 5.5 # to make the y scale similar to the other plot: 
p = in_phase_vs_not_in_phase_enrichment_tss %>% 
    dplyr::left_join(., y = TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_tss) %>% 
    dplyr::mutate(foldChange = (rhythmic_genes_in_phases_bound/(rhythmic_genes_not_in_phases_bound+1))/(n_kzfps_in_phases/n_kzfps)) %>% 
    arrange(p_enrich) %>% 
    ggplot(aes(x = log2(foldChange), y = -log10(p_enrich), col = padj_enrich < 0.05, label = TF_symbol)) + 
    geom_point() + 
    ylim(c(0, ylim_high)) +
    theme_classic() +
    scale_color_manual(values = c("grey60", "black")) + 
    theme(legend.position = "none") + 
    ggrepel::geom_label_repel(size = 3, label.size = 0, fill = alpha(c("white"),0.5)) + 
    ggtitle(paste0("TFs enriched over promoters of\nM-to-G1/S RIF1-low KZFPs"))
p
```

```R
ggsave("../out/cell_cycle_figures/TFs_enriched_at_KZFPs_promoters_M_to_G1S.svg", p, svg, height = 3.5, width = 3)
```

```R
# aLL kzfps
phases = list(c("M", "eG1", "lG1", "G1/S", "S1"), c( "S2", "S/G2", "G2"))
rhythm_threshold = 0.1
i = 1
# phase_list = rhythm_all_genes$phase_assigned %>% unique()
df_list = list()
n_kzfps_in_phases = kzfps_cycling_info %>% dplyr::filter(!cluster %in% rif1_high_clusters, padj <= rhythm_threshold, phase_assigned %in% phases[[i]]) %>% nrow()
n_kzfps_in_phases
n_kzfps = kzfps_cycling_info %>% dplyr::filter(!cluster %in% rif1_high_clusters, padj <= rhythm_threshold) %>% nrow()
n_kzfps
```

```R
TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_znfs = kzfp_znf_centric_metrics %>% dplyr::filter(!cluster %in% rif1_high_clusters, znf_padj <= rhythm_threshold) %>% dplyr::group_by(TF_symbol) %>% dplyr::summarize(rhythmic_genes_in_phases_bound = sum(((znf_phase_assigned %in% phases[[i]]))),
                                                                        rhythmic_genes_not_in_phases_bound = sum(((!znf_phase_assigned %in% phases[[i]])), na.rm = T))
```

```R
TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_znfs %>% head()
```

```R
TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_znfs %>% dplyr::filter(TF_symbol %in% c("ZNF274", "TRIM28"))
```

```R
TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_znfs %>% dplyr::filter(TF_symbol %in% c("ZNF274", "TRIM28", "ZNF334"))
```

```R
enrichment_test_phases_vs_other_rhythmic_kzfps <- function(row) {
    x_observed = as.integer(row[["rhythmic_genes_in_phases_bound"]])
    m = n_kzfps_in_phases
    n = n_kzfps - n_kzfps_in_phases
    k = as.integer(row[["rhythmic_genes_in_phases_bound"]]) + as.integer(row[["rhythmic_genes_not_in_phases_bound"]])
    
    x = 0:m # is the variable tested
    probs <- dhyper(x, m, n, k, log = FALSE)
    
    # we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
    pval_one_sided = sum(probs[x>=x_observed])

    # we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
    pval_two_sided = sum(probs[probs <= pval_one_sided])
    
    return(pval_two_sided)
}
```

```R
in_phase_vs_not_in_phase_enrichment_znfs = TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_znfs %>% apply(., 1, enrichment_test_phases_vs_other_rhythmic_kzfps) %>% as.data.frame()
```

```R
in_phase_vs_not_in_phase_enrichment_znfs$TF_symbol = TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_znfs$TF_symbol
colnames(in_phase_vs_not_in_phase_enrichment_znfs)[1] = 'p_enrich'
in_phase_vs_not_in_phase_enrichment_znfs %>% head()
```

```R
in_phase_vs_not_in_phase_enrichment_znfs$is_kzfp = ifelse(in_phase_vs_not_in_phase_enrichment_znfs$TF_symbol %in% kzfps_cycling_info$symbol, T, F)
in_phase_vs_not_in_phase_enrichment_znfs$padj_enrich = p.adjust(in_phase_vs_not_in_phase_enrichment_znfs$p_enrich, method = "BH")
in_phase_vs_not_in_phase_enrichment_znfs = in_phase_vs_not_in_phase_enrichment_znfs %>% 
    dplyr::left_join(., rhythm_all_genes %>% dplyr::select(symbol, padj, phase_assigned, stars, RepliTiming) %>% dplyr::rename(TF_symbol = symbol, padj_rhythm = padj, stars_rhythm = stars))
```

```R
peaks_x_kzfp_znfs %>% dplyr::pull(TF_symbol) %>% unique() %>% length()
```

```R
p = in_phase_vs_not_in_phase_enrichment_znfs %>% 
    dplyr::left_join(., y = TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_znfs) %>% 
    dplyr::mutate(foldChange = (rhythmic_genes_in_phases_bound/(rhythmic_genes_not_in_phases_bound+1))/(n_kzfps_in_phases/n_kzfps)) %>% 
    arrange(p_enrich) %>% 
    ggplot(aes(x = log2(foldChange), y = -log10(p_enrich), col = padj_enrich < 0.05, label = TF_symbol)) + 
    geom_point() + 
    theme_classic() +
    ylim(c(0, ylim_high)) +
    scale_color_manual(values = c("grey60", "black")) + 
    theme(legend.position = "none") + 
    ggrepel::geom_label_repel(size = 3, label.size = 0, fill = alpha(c("white"),0.5)) + 
    ggtitle(paste0("TFs enriched over ZNFs of\nM-to-G1/S RIF1-low KZFPs"))
```

```R
in_phase_vs_not_in_phase_enrichment_znfs %>% arrange(p_enrich) %>% head()
```

```R
ggsave("../out/cell_cycle_figures/TFs_enriched_at_KZFPs_ZNFs_M_to_G1S.svg", p, svg, height = 3.5, width = 3)
```

```R
n_kzfps_in_phases
```

```R
# merging the two figures
enrichment_tss_znfs_at_kzfps = rbind(in_phase_vs_not_in_phase_enrichment_tss %>%
                                         dplyr::mutate(enrichment_at = "promoters"),
                                     in_phase_vs_not_in_phase_enrichment_znfs %>%
                                         dplyr::mutate(enrichment_at = "ZNFs"))
TF_in_M_G1_KZFPs_vs_rhythmic_KZFPs = rbind(TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_tss %>%
                                           dplyr::mutate(enrichment_at = "promoters"),
                                           TF_in_M_G1_KZFPs_vs_rhytmic_KZFPs_znfs %>%
                                           dplyr::mutate(enrichment_at = "ZNFs"))
p = enrichment_tss_znfs_at_kzfps %>% 
    dplyr::left_join(., y = TF_in_M_G1_KZFPs_vs_rhythmic_KZFPs) %>% 
    dplyr::mutate(foldChange = (rhythmic_genes_in_phases_bound/(rhythmic_genes_not_in_phases_bound+1))/(n_kzfps_in_phases/n_kzfps)) %>% 
    arrange(p_enrich) %>% 
    ggplot(aes(x = log2(foldChange), y = -log10(p_enrich), col = p_enrich < 0.05, label = TF_symbol)) + 
    geom_point() + 
    theme_classic() +
    ylim(c(0, ylim_high)) +
    scale_color_manual(values = c("grey60", "black")) + 
    theme(legend.position = "none") + 
    ggrepel::geom_label_repel(size = 3, label.size = 0, fill = alpha(c("white"),0.5)) + 
    ggtitle(paste0("TFs enriched at\nM-to-G1/S RIF1-low KZFPs")) +   
    facet_wrap("enrichment_at", 1, scales = "free")
p
```

ZNF274 is the top hit for enrichment at M-S1 KZFPs. Although the enrichment is barely signif. and we had to take a padj_rhythm < 0.1 instead of 0.05, we still obtain it as the best enriched TF. 

```R
ggsave("../out/cell_cycle_figures/TFs_enriched_at_KZFPs_M_to_G1S_promoters_and_ZNFs.svg", p, svg, height = 3.5, width = 4)
```
