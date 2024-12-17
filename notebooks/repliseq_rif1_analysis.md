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

<!-- #region -->
# Assessing the enrichment of KZFP binding, in particular ZNF274, in RIF1-KO unchanged late replicating regions

## Analysis plan

1) identify the unchanged late replicating regions upon RIF1 KO in HCT cells.

2) check ZNF274 enrichment at those regions. Is ZNF274 the strongest enriching KZFP? TF? Similarly, check delta H3K9me3 upon ZNF274 KO.

3) identify the late replicating regions shifting to earlier replication upon SETDB1 and SUV39H1/2 KD in RIF1KO cells. Is ZNF274 enriched ? the most enriched KZFP? The most enriched TF?


## Decisions
- We're staying at the 50kb bin resolution for all enrichment tests, and will make sure the bin-to-feature or features-to-bin assignations are 1-to-1.
- For peaks, we'll only keep the RTbin with the greatest overlap, not all RTbins.

## observations from looking at ZNF274 binding, TRIM28 binding (Jakobsson ineurons, Flavia's CD4 resting and activated T cells), RIF1 binding and late unchanged regions
The association between ZNF274 binding and being a late unchanged region appears most clearly on chromosome 19. On other chromosomes, KZFPs are somewhat anectodal in term of gene content, and the relationship does not seem to hold.

Isolated ZNF274 peaks do not seem to associate with late unchanged replication timing.

Doesnt work for PCDH genes

On chromosome 19, by eyeballing: 7 unchanged late replicating regions and 6 of them are ZNF274 high. The relationship does not seem to hold on other chromosomes, which are much poorer in KZFP genes.

There's also the subtelomeric arm of chromosome 4, with ZNF595 etc, which has high ZNF274 and late replicating unchanged upon RIF1KO

Also the subtelomeric arm of the long arm of chromosome 5, with ZNF354B etc.

The RIF1 high KZFP cluster on chromosome 6 goes from early to a bit less early. 

The centromeric KZFP cluster on chromosome 7 with ZNF680 etc is a late unchanged region

The centromeric ZNF274 binding rich region on chromosome 7 with ZNF479 (not in a cluster) is also late unchanged

Doesnt really work for the subtelomeric region of the long arm of chr8 (ZNF517) though it makes sense since these are ZNF274 indep.

Works for the centromeric region of chrom10, although no KZNFs here.

Works for the subtelomeric region of the long arm of chr12.

Pericentromeric region of chr16, chr17, chr21, chr22

## Things to check:
Mouse KZFPs and mouse replitiming changes upon RIF1KO
<!-- #endregion -->

## Identifying control samples

High resolution RIF1KO repliseq samples at in the GEO of Klein et al. The control are not specificed, but absent, suggesting that they used Zhao et al. Genome Biology 2020 as a control. 

The standard repli-seq control for HCT cells was not found, but we speculate they used the untreated degron-fused RIF1 cell line (RIF-AID, untreated). 

```R
library(tidyverse)
library(GEOquery)
library(R.utils) # to gunzip
library(spatstat) # gaussian smoothing
library(ComplexHeatmap) # plotting the read densities
library(circlize)
library(mclust) # for gaussian kernel 1D clustering of RTdiff values as in Klein et al., Science 2021
```

## Data download

```R
dir.create('../data/repliseq_rif1/')
GEOquery::getGEOSuppFiles(GEO = "GSE160563", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE160563_HCT116_KO_log2EL_RT.bigwig")
```

```R
GEOquery::getGEOSuppFiles(GEO = "GSE160563", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE160563_HCT116_RIF1AID_Ctl_log2EL_RT.bigwig")
```

```R
GEOquery::getGEOSuppFiles(GEO = "GSE160563", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE160563_HCT116_RIF1KO_Rep1_highres_Repliseq_hg38.mat.gz")

```

```R
# H9 RIF1KO
GEOquery::getGEOSuppFiles(GEO = "GSE160563", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE160563_H9_RIF1KO_Rep1_highres_Repliseq_hg38.mat.gz")

```

```R
GEOquery::getGEOSuppFiles(GEO = "GSE160563", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE160563_HCT_RIF1KO_H3K9me3KD_KDsiRNA_log2EL_RT.bigwig")
GEOquery::getGEOSuppFiles(GEO = "GSE160563", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE160563_HCT_RIF1KO_H3K9me3KD_ctrlsiRNA_log2EL_RT.bigwig")

```

```R
# HCT WT highres repli-seq control
GEOquery::getGEOSuppFiles(GEO = "GSE137764", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE137764_HCT_GaussiansGSE137764_mooth_scaled_autosome.mat.gz")

```

```R
# H1 WT highres H9 repli-seq control
GEOquery::getGEOSuppFiles(GEO = "GSE137764", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE137764_H9_GaussiansGSE137764_mooth_scaled_autosome.mat.gz")

# H9 WT H1 highres repli-seq control
GEOquery::getGEOSuppFiles(GEO = "GSE137764", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE137764_H1_GaussiansGSE137764_mooth_scaled_autosome.mat.gz")


# mouse ESCs
GEOquery::getGEOSuppFiles(GEO = "GSE137764", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE137764_mESC_Gaussiansmooth_scaled_autosome.mat.gz")


# mouse NPCs
GEOquery::getGEOSuppFiles(GEO = "GSE137764", makeDirectory = F, baseDir = '../data/repliseq_rif1', fetch_files = T, filter_regex = "GSE137764_mNPC_Gaussiansmooth_scaled_autosome.mat.gz")
```

Missing data: the RIF1KO repli-seq from Cornacchia et al., 2012. Without this, we can't verify whether RIF1KO also leaves KZFPs relatively untouched compared to other genes.

```R
to_unzip = c("GSE137764_HCT_GaussiansGSE137764_mooth_scaled_autosome.mat.gz",
             "GSE160563_HCT116_RIF1KO_Rep1_highres_Repliseq_hg38.mat.gz",
             "GSE160563_HCT116_RIF1KO_Rep1_highres_Repliseq_hg38.mat.gz",
             "GSE137764_H9_GaussiansGSE137764_mooth_scaled_autosome.mat.gz", 
"GSE137764_H1_GaussiansGSE137764_mooth_scaled_autosome.mat.gz",
"GSE137764_mESC_Gaussiansmooth_scaled_autosome.mat.gz",
"GSE137764_mNPC_Gaussiansmooth_scaled_autosome.mat.gz",
            "GSE160563_H9_RIF1KO_Rep1_highres_Repliseq_hg38.mat.gz")

for (f in to_unzip) {
    gunzip(paste0("../data/repliseq_rif1/", f), remove=FALSE, overwrite = T)
}
```

## Reading 50-kb bigwig repli-seq data


There are two data sources in this study: the Zhao paper, and the Klein paper. The data was not processed exactly in the same way in both. We hardcode reading functions to account for that


```R
read_hires_repliseq_zhao <- function(path){
    mat = read.table(path, sep = '\t', header = F) %>% t()
    colnames(mat) = c("chrom", "start", "end", paste0("S", 1:16))
    
    # normalizing and replacing NAs by zeros
    temp_mat = mat[, paste0("S", 1:16)] %>% apply(., 2, as.numeric)
    temp_mat[which(is.na(temp_mat))] = 0
    temp_mat = temp_mat %>% apply(., 1, function(x) x*100/sum(x)) %>% t() # normalizing to percentage
    
    mat = cbind(mat[, 1:3], temp_mat %>% as.data.frame()) %>% 
                                  dplyr::mutate(across(start:end, as.integer)) %>%
                                  dplyr::mutate(across(S1:S16, as.double))
    
    }
                                  
read_hires_repliseq_klein <- function(path){ # only difference is that there is no need to normalize to 100
    mat = read.table(path, sep = '\t', header = F) %>% t()
    colnames(mat) = c("chrom", "start", "end", paste0("S", 1:16))
    
    # normalizing and replacing NAs by zeros
    temp_mat = mat[, paste0("S", 1:16)] %>% apply(., 2, as.numeric)
    temp_mat[which(is.na(temp_mat))] = 0    
    mat = cbind(mat[, 1:3], temp_mat %>% as.data.frame()) %>% 
                                  dplyr::mutate(across(start:end, as.integer)) %>%
                                  dplyr::mutate(across(S1:S16, as.double))
    
    }
                                           
draw_hires_repliseq <- function(repliseq, chr, coords) {
    region = repliseq %>% dplyr::filter(chrom == chr, start >= coords[1], end <= coords[2]) %>% arrange(start, end)
    
    mat = region %>% dplyr::select(S1:S16) %>% as.matrix() %>% t() 
    colnames(mat) = NULL

    p = mat %>% ComplexHeatmap::Heatmap(., cluster_rows = F, cluster_columns = F)
    
    return(p)
    
    }
```

Additionally, we process the bedgraph converted from bigwigs for the Hansen et al. study, see script `/scripts/download_repliseq_K562.sh`. The data are percentage-normalized reads in each cell fraction. Thus, we can reconstruct a very similar matrix as the one obtained for the hi-resolution repliseq, except that we'll have six fractions instead of twelve

Columns will be:
`chrom`, `start`, `end`, `S1`... `S6`

```R
read_repliseq_hansen <- function(dir_path) {
    df_list = list()
    sample_order = c("G1", "S1", "S2", "S3", "S4", "G2")
    for (i in 1:length(sample_order)) {
        df_list[[i]] = read.table(paste0(dir_path, "wgEncodeUwRepliSeqK562", sample_order[[i]], "PctSignalRep1.bigWig.bedgraph"), sep = '\t', header = F)
    }
    repliseq = df_list %>% reduce(full_join, by = c("V1", "V2", "V3"))
    colnames(repliseq) = c("chrom", "start", "end", paste0("S", 1:length(sample_order)))
    fill_na_cols = as.list(rep(0, times = ncol(repliseq)-3))
    names(fill_na_cols) = paste0("S", 1:length(sample_order))
    repliseq = repliseq %>% 
        tidyr::replace_na(replace = fill_na_cols) %>%
        arrange(chrom, start)
    }
```

```R
K562 = read_repliseq_hansen('../data/repliseq_k562/bedgraph/')
```

```R
K562 %>% dplyr::select(S1:S6) %>% rowSums()
```

```R
not100 = ((K562 %>% dplyr::select(S1:S6) %>% rowSums()) != 100)
```

```R
sum(not100)/nrow(K562)
```

For most rows the values do not sum to 100. How come? Will we have to renormalize it ourselves?

```R
read_repliseq_293T <- function(dir_path) {
    df_list = list()
    sample_order = c("G1_S1", "S1_S2", "S2_S3", "G2_S4")
    for (i in 1:length(sample_order)) {
        df_list[[i]] = read.table(paste0(dir_path, "WT_", sample_order[[i]], "_noWindow.bw.bedgraph"), sep = '\t', header = F)
    }
    repliseq = df_list %>% reduce(full_join, by = c("V1", "V2", "V3"))
    colnames(repliseq) = c("chrom", "start", "end", paste0("S", 1:length(sample_order)))
    fill_na_cols = as.list(rep(0, times = ncol(repliseq)-3))
    names(fill_na_cols) = paste0("S", 1:length(sample_order))
    repliseq = repliseq %>% 
        tidyr::replace_na(replace = fill_na_cols) %>%
        arrange(chrom, start)
    }


read_repliseq_274KO <- function(dir_path) {
    df_list = list()
    sample_order = c("G1_S5", "S1_S6", "S2_S7", "G2_S8")
    for (i in 1:length(sample_order)) {
        df_list[[i]] = read.table(paste0(dir_path, "274_", sample_order[[i]], "_noWindow.bw.bedgraph"), sep = '\t', header = F)
    }
    repliseq = df_list %>% reduce(full_join, by = c("V1", "V2", "V3"))
    colnames(repliseq) = c("chrom", "start", "end", paste0("S", 1:length(sample_order)))
    fill_na_cols = as.list(rep(0, times = ncol(repliseq)-3))
    names(fill_na_cols) = paste0("S", 1:length(sample_order))
    repliseq = repliseq %>% 
        tidyr::replace_na(replace = fill_na_cols) %>%
        arrange(chrom, start)
    }
```

```R
hek293t_wt = read_repliseq_293T('../data/repliseq_274ko/bedgraph/')
hek293t_274ko = read_repliseq_274KO('../data/repliseq_274ko/bedgraph/')
```

# Data visualization

## Reproducing Fig. 1B) from Klein et al.


Encapsulated as a function:

```R
format_col_names_coords <- function(start_pos, chr) {
    to_label = start_pos %>% quantile(., probs = c(0.2, 0.4, 0.6, 0.8), type = 3)
    idx_kept = match(to_label, start_pos)
    ndigits = 3
    if(range(start_pos)[2]-range(start_pos)[1] < 5e6) {
        ndigits = 4}
    coords_label = character(length = length(start_pos))
    coords_label[idx_kept] = as.character(format(to_label/1e6, digits = ndigits))
    coords_label[idx_kept[1]] = paste0(coords_label[idx_kept[1]], "Mb")
    coords_label[1] = chr
    return(coords_label)
    }
```

```R
draw_repliseq <- function(repliseq, chr, coords, perc_range = c(0, 20), use_raster = F) {    
    region = repliseq %>% dplyr::filter(chrom == chr, start >= coords[1], end <= coords[2]) %>% arrange(start, end)
        
    plot_start = (coords[1]%/%50000)*50000
    if (plot_start < 0) {plot_start = 0}
    plot_end = (coords[2]%/%50000+1)*50000
    
    all_coords = data.frame(start = seq(from = plot_start, to = plot_end - 50000, by = 50000),
                            end = seq(from = plot_start + 50000, to = plot_end, by = 50000)) %>%
        dplyr::mutate(chrom = chr)
        
    mat = all_coords %>% dplyr::left_join(., region) %>% dplyr::select("S1":tail(colnames(region), n = 1)) %>% as.matrix() %>% t() 
    
    # genomic coordinates for x axis legend
    coord_labels = format_col_names_coords(all_coords %>% pull(start), chr)
    ha = columnAnnotation(foo = anno_text(coord_labels, just = "center", rot = 0, location = 0.5))

    p = mat %>% ComplexHeatmap::Heatmap(., 
                                        cluster_rows = F,
                                        cluster_columns = F, 
                                        col = circlize::colorRamp2(seq(from = perc_range[1], to = perc_range[2], length.out = 5), viridis::cividis(5)),
                                       row_title = paste0(tail(colnames(region), n = 1)," <-- S1"),
                                       show_row_names = F,
                                        show_column_names = F,
                                        bottom_annotation = ha,
                                        column_names_rot = 0,
                                        heatmap_legend_param = list(title = "perc. repl."),
                                        use_raster = use_raster,
                                        raster_device = "png"
                                       )
    return(p)
    
}
```

```R
HCT_WT_hires = read_hires_repliseq_zhao("../data/repliseq_rif1/GSE137764_HCT_GaussiansGSE137764_mooth_scaled_autosome.mat")
```

```R
draw_repliseq(HCT_WT_hires, "chr1", c(170*1e6, 200*1e6))
```

```R
K562 %>% head()
```

Clearly this signal has not been smoothed...

Let's try gaussian smoothing using the spatstat package


Moreover, the signal is in 1kb bins, which is probably not that useful and makes computations much slower. We can average it in 50kb bins and then smoothen it.

```R
downsample_repliseq_hansen <- function(repliseq, bin_size = 50000) {
    repliseq_downsampled_chrom_list = list()
    chrom_list = unique(repliseq$chrom)
    for (chr in chrom_list) {
        repliseq_chrom = repliseq %>% dplyr::filter(chrom == chr)
        breaks = seq(from = 1, to = repliseq_chrom$start[length(repliseq_chrom$start)]+bin_size, by = bin_size)
        repliseq_chrom$bins = cut(repliseq_chrom$start, breaks = breaks)
        repliseq_downsampled_chrom_list[[chr]] = repliseq_chrom %>% 
                                                        dplyr::group_by(bins) %>% 
                                                        dplyr::reframe(chrom = unique(chrom),
                                                              start = 1+(as.numeric(bins)-1)*50000,
                                                              end = as.numeric(bins)*50000,
                                                              S1 = mean(S1),
                                                              S2 = mean(S2),
                                                              S3 = mean(S3),
                                                              S4 = mean(S4),
                                                              S5 = mean(S5),
                                                              S6 = mean(S6)) %>% 
                                                        unique() %>% 
                                                        dplyr::select(-bins)
        }
    return(do.call(rbind, repliseq_downsampled_chrom_list))
    }
```

```R
K562_ds = downsample_repliseq_hansen(K562)
```

```R
K562_ds %>% head()
```

```R
# harmonizing the start coordinate with the other repliseq experiments:
K562_ds = K562_ds %>% dplyr::mutate(start = start-1)
```

```R
draw_repliseq(K562_ds, "chr1", c(170*1e6, 200*1e6))
```

That already looks pretty good, we can smooth once and renormalize to see what we would get.

each chromosome is like one image, and we're gonna have to pad them.

The matrix to be padded is organized as follows: each fraction is a *row*, each column is a *window*. In that formulation, the first and last rows (first and last fractions) are padded, i.e. duplicated at the beginning and end of the matrix.

The code for the padding can be found here: https://github.com/oliviacamel/RIF1_KO_analysis/blob/main/high_res_repliseq/make_array.py

```R
smooth_repliseq_hansen <- function(repliseq) {
    repliseq_smoothed_dflist = list()
    for (chr in unique(repliseq$chrom)) {
        chrom_mat = repliseq %>% dplyr::filter(chrom == chr) %>% dplyr::select(S1:S6) %>% as.matrix
        chrom_mat_padded = cbind(chrom_mat[, 1], chrom_mat, chrom_mat[, ncol(chrom_mat)])
        chrom_mat_smoothed = spatstat.explore::blur(as.im(chrom_mat_padded), sigma = 1, kernel = "gaussian", normalise = F, bleed = F)
        repliseq_smoothed_dflist[[chr]] = cbind(repliseq %>% dplyr::filter(chrom == chr) %>% dplyr::select(1:3), as.data.frame(chrom_mat_smoothed[, 2:7]))
        }
    repliseq_smoothed = do.call(rbind, repliseq_smoothed_dflist)
    colnames(repliseq_smoothed) = c("chrom", "start", "end", "S1", "S2", "S3", "S4", "S5", "S6")
    rownames(repliseq_smoothed) = NULL
    repliseq_smoothed = repliseq_smoothed %>% dplyr::mutate(across(c("start", "end"), as.integer))
    return(repliseq_smoothed)
    }
```

```R
K562_smoothed = smooth_repliseq_hansen(K562_ds)
```

```R
K562_smoothed %>% head()
```

```R
draw_repliseq(K562_smoothed, "chr1", c(170*1e6, 200*1e6))
```

```R
normalize_repliseq_hansen <- function(repliseq) {
    temp_mat = repliseq %>% dplyr::select(S1:S6)
    return(cbind(repliseq %>% dplyr::select(1:3),
                 ((temp_mat)/rowSums(temp_mat)*100) %>% as.data.frame()))
    }
```

```R
K562_smoothed_norm = normalize_repliseq_hansen(K562_smoothed)
```

```R
K562_smoothed_norm %>% head()
```

```R
K562_smoothed_norm %>% dplyr::select(S1:S6) %>% rowSums() %>% head()
```

```R
draw_repliseq(K562_smoothed_norm, "chr1", c(170*1e6, 200*1e6))
```

That looks very good to me. Maybe a bit overexposed

```R
draw_repliseq(K562_smoothed_norm, "chr1", c(170*1e6, 200*1e6), perc_range = c(0, 40))
```

This looks just like the HCT116 repli-seq.

```R
HCT_RIF1KO_hires = read_hires_repliseq_klein("../data/repliseq_rif1/GSE160563_HCT116_RIF1KO_Rep1_highres_Repliseq_hg38.mat")
```

```R
HCT_RIF1KO_hires %>% dplyr::select(S1:S16) %>% rowSums() %>% head()
```

Also properly normalized

```R
p = draw_repliseq(HCT_RIF1KO_hires, "chr1", c(170*1e6, 200*1e6))
p
```

It now seems we're working on the correct data as we reobtained the plots from fig. 1B in the Klein 2021 Science paper

```R
# processing our own 293T cell repli-seq data
smooth_repliseq_293T <- function(repliseq) {
    repliseq_smoothed_dflist = list()
    for (chr in unique(repliseq$chrom)) {
        chrom_mat = repliseq %>% dplyr::filter(chrom == chr) %>% dplyr::select(S1:S4) %>% as.matrix()
        chrom_mat_padded = cbind(chrom_mat[, 1], chrom_mat, chrom_mat[, ncol(chrom_mat)])
        chrom_mat_smoothed = spatstat.explore::blur(as.im(chrom_mat_padded), sigma = 1, kernel = "gaussian", normalise = F, bleed = F)
        repliseq_smoothed_dflist[[chr]] = cbind(repliseq %>% dplyr::filter(chrom == chr) %>% dplyr::select(1:3), as.data.frame(chrom_mat_smoothed[, 2:5]))
        }
    repliseq_smoothed = do.call(rbind, repliseq_smoothed_dflist)
    colnames(repliseq_smoothed) = c("chrom", "start", "end", "S1", "S2", "S3", "S4")
    rownames(repliseq_smoothed) = NULL
    repliseq_smoothed = repliseq_smoothed %>% dplyr::mutate(across(c("start", "end"), as.integer))
    return(repliseq_smoothed)
    }
```

```R
hek293t_wt %>% head()
```

```R
p = draw_repliseq(hek293t_wt, "chr1", c(170*1e6, 200*1e6), perc_range = c(0, 3000))
p
```

```R
hek293t_wt = hek293t_wt %>% dplyr::filter(! chrom %in% c("chrM", "chrY"))
hek293t_smoothed = smooth_repliseq_293T(hek293t_wt)
```

```R
normalize_repliseq_293T <- function(repliseq) {
    temp_mat = repliseq %>% dplyr::select(S1:S4)
    return(cbind(repliseq %>% dplyr::select(1:3),
                 ((temp_mat)/rowSums(temp_mat)*100) %>% as.data.frame()))
    }
```

```R
hek293t_smoothed_norm = normalize_repliseq_293T(hek293t_smoothed)
```

```R
p = draw_repliseq(hek293t_smoothed_norm, "chr1", c(170*1e6, 200*1e6), perc_range = c(0, 40))
p
```

Same pattern as in K562 and HCT116

```R
# processing the HEK293T ZNF274 KO repli-seq data
hek293t_274ko_smoothed = hek293t_274ko %>% 
    dplyr::filter(! chrom %in% c("chrM", "chrY")) %>%
    smooth_repliseq_293T()
```

```R
hek293t_274ko_smoothed_norm = normalize_repliseq_293T(hek293t_274ko_smoothed)
hek293t_274ko_smoothed_norm %>% head()
```

```R
p = draw_repliseq(hek293t_274ko_smoothed_norm, "chr1", c(170*1e6, 200*1e6), perc_range = c(0, 40))
p
```

Looks exactly as the 293T WT


## Export repliseq heatmaps in bedgraph

```R
dir.create('../out/bedgraph_fractions', recursive = T)
```

```R
# exporting the RTindex to bedgraph file for exploration
RT_to_bedgraph <- function(RT, path, name, description = "none") {
    track_def_line = paste0("track type=bedGraph name=", name, " description=", description, " visibility=display_mode color=r,g,b altColor=r,g,b priority=priority autoScale=off")
    rbind(c(track_def_line, "", "", ""), RT) %>% write.table(., path, quote=F, sep = '\t', row.names = F, col.names = F)
    }
```

```R
fractions_to_bedgraph <- function(fractions, name, dir) {
    for (f in colnames(fractions)[4:ncol(fractions)]) {
        f_df_temp = cbind((fractions %>% dplyr::select(1:3)), fractions[f])
        RT_to_bedgraph(f_df_temp, path = paste0(dir, "/", name, "_", f, ".bedgraph"), name = name)
        }
    }

```

```R
fractions_to_bedgraph(hek293t_274ko_smoothed_norm, "hek293t_274ko_smoothed_norm", "../out/bedgraph_fractions/")
```

```R
fractions_to_bedgraph(hek293t_smoothed_norm, "hek293t_wt_smoothed_norm", "../out/bedgraph_fractions/")
```

```R
fractions_to_bedgraph(HCT_WT_hires, "hct_wt_smoothed_norm", "../out/bedgraph_fractions/")
```

```R
fractions_to_bedgraph(HCT_RIF1KO_hires, "hct_rif1ko_smoothed_norm", "../out/bedgraph_fractions/")
```

```R
fractions_to_bedgraph(K562_smoothed_norm, "k562_wt_smoothed_norm", "../out/bedgraph_fractions/")
```

## Extracting KZFP clusters and solo KZFPs from Jonas's KZFPs

```R
kzfps_jonas = read.table('../data/kzfps_jonas.csv', sep = ';', header = T, quote = "")
kzfps_jonas %>% head()
kzfps_jonas %>% dim()
```

TODO recontextualize: Extracting KZFP clusters:

```R
kzfps_jonas %>% dplyr::filter(cluster!="noCluster") %>% 
    dplyr::select(chr, start, end, cluster) %>% 
    arrange(cluster, start, end) %>%
    group_by(cluster) %>% 
    dplyr::mutate(cluster_start = min(start),
                 cluster_end = max(end)) %>%
    dplyr::select(chr, cluster_start, cluster_end, cluster) %>%
    unique() %>%
    write.table('../out/kzfps_clusters_jonas.bed', col.names = F, row.names = F,quote = F, sep = '\t')
```

```R
kzfps_jonas %>% dplyr::filter(cluster=="noCluster") %>% 
    dplyr::select(chr, start, end, assigned_gene) %>% 
    arrange(chr, start, end) %>%
    unique() %>%
    write.table('../out/kzfps_noCluster_jonas.bed', col.names = F, row.names = F,quote = F, sep = '\t')
```

<!-- #region -->
The exploration promises to be nice: it seems like even small KZFP clusters can be replicated quite late, which is not necessarily captured by the low-res repliseq.

## Identifying RIF1KO-unchanged late replicated regions

### Computing RT indices for all files
From the Klein et al. Supplementals: 

![image.png](attachment:a5510833-d706-4259-8b20-e00bfd704731.png)

For the plotting, rows full of NA were not an issue, but now they certainly are. We keep the intersection of all rows that are not full of NAs, separately for human and mouse data.

Human data in hi-resolution, 12-fraction repli-seq
- HCT_WT_hires
- HCT_RIF1KO_hires
- H1_WT_hires
- H9_WT_hires
- H9_RIF1KO_hires
Human data with 6 fractions
- K562
Human data with four fractions:
- HEK293T_WT
- HEK293T_ZNF274KO


Mouse data with 12 fractions
- mESCs_WT_hires
- mNPCs_WT_hires

## RT indices for human hi-resolution (16-fraction) data
<!-- #endregion -->

```R
H1_WT_hires = read_hires_repliseq_zhao("../data/repliseq_rif1/GSE137764_H1_GaussiansGSE137764_mooth_scaled_autosome.mat")
H9_WT_hires = read_hires_repliseq_zhao("../data/repliseq_rif1/GSE137764_H9_GaussiansGSE137764_mooth_scaled_autosome.mat")
```

```R
H9_RIF1KO_hires = read_hires_repliseq_klein("../data/repliseq_rif1/GSE160563_H9_RIF1KO_Rep1_highres_Repliseq_hg38.mat")
```

```R
# ordering them the same way to check that they have the exact same locations
HCT_WT_hires = HCT_WT_hires %>% dplyr::arrange(chrom, start, end)
HCT_RIF1KO_hires = HCT_RIF1KO_hires %>% dplyr::arrange(chrom, start, end)
H1_WT_hires = H1_WT_hires %>% dplyr::arrange(chrom, start, end)
H9_WT_hires = H9_WT_hires %>% dplyr::arrange(chrom, start, end)
H9_RIF1KO_hires = H9_RIF1KO_hires %>% dplyr::arrange(chrom, start, end)

```

```R
HCT_RIF1KO_hires %>% nrow
HCT_WT_hires %>% nrow
H1_WT_hires %>% nrow
H9_WT_hires %>% nrow
H9_RIF1KO_hires %>% nrow
```

```R
stopifnot(identical(HCT_WT_hires %>% dplyr::select(1:3),
                HCT_RIF1KO_hires %>% dplyr::select(1:3)))

stopifnot(identical(H9_WT_hires %>% dplyr::select(1:3),
                H9_RIF1KO_hires %>% dplyr::select(1:3)))

stopifnot(identical(H1_WT_hires %>% dplyr::select(1:3),
                H9_RIF1KO_hires %>% dplyr::select(1:3)))

stopifnot(identical(HCT_WT_hires %>% dplyr::select(1:3),
                H9_RIF1KO_hires %>% dplyr::select(1:3)))
```

```R
HCT_WT_hires %>% dplyr::select(1:3) %>% head()
```

```R
HCT_RIF1KO_hires %>% dplyr::select(1:3) %>% head()
```

### Filtering full NA rows in all 16-fraction repliseq matrices

```R
kept_Gilbert = rep(T, times = nrow(HCT_WT_hires))

for (m in c("HCT_WT_hires", "HCT_RIF1KO_hires", "H1_WT_hires", "H9_WT_hires", "H9_RIF1KO_hires")) {
    kept_Gilbert = kept_Gilbert & (get(m) %>% dplyr::select(S1:S16) %>% apply(., 1, function(x) sum(is.na(x))) != 16)
    }
```

```R
sum(kept_Gilbert)/(HCT_WT_hires %>% nrow)
```

```R
sum(!kept_Gilbert)
```

```R
# computing the RT index for any number of fractions
computeRTindex <- function(repliseq_bin, nfractions = 16) {
    rtidx = log2(repliseq_bin[1:(nfractions/2)]%*%((nfractions/2):1))-log2(repliseq_bin[(nfractions/2+1):nfractions]%*%(1:(nfractions/2)))
    }
```

```R
# building a bed-like RTindex repliseq DF, with the three firt columns being the window chrom, start, end and then the RTindex.
assemble_RTindex_df <- function(repliseq, kept, nfractions = 16, pseudocount = 0.01) {
    ncol_repliseq = ncol(repliseq)
    RTindex = apply(repliseq %>% dplyr::filter(kept) %>% dplyr::select(all_of((ncol_repliseq-nfractions+1):ncol_repliseq)) %>% as.matrix() + pseudocount, 1, computeRTindex, nfractions = nfractions)
    return(cbind((repliseq %>% dplyr::filter(kept))[1:3], RTindex))
}
```

```R
# computing the difference in RTindex between two samples with the same number of fractions
assemble_RTdiff_df <- function(repliseq1, repliseq2, kept, pseudocount = 0.01) {
    RTindex1 = assemble_RTindex_df(repliseq1, kept = kept)
    RTindex2 = assemble_RTindex_df(repliseq2, kept = kept)
    
    return(cbind(RTindex1[1:3], RTindex2$RTindex - RTindex1$RTindex) %>% dplyr::rename(RTdiff = 4))
}
```

```R
dir.create("../out/bedgraph", showWarnings = F)
```

```R
for (m in c("HCT_WT_hires", "HCT_RIF1KO_hires", "H1_WT_hires", "H9_WT_hires", "H9_RIF1KO_hires")) {
    RTindex = assemble_RTindex_df(get(m), kept_Gilbert, nfractions = 16)
    RT_to_bedgraph(RTindex, paste0("../out/bedgraph/RTindex_", m, ".bedgraph"), m)
    }
```

```R
RT_to_bed <- function(RT, path) {
    RT %>% dplyr::arrange(chrom, start, end) %>% write.table(path, col.names = F, row.names = F, sep = '\t', quote = F)
    }
dir.create('../out/bed/')
```

```R
for (m in c("HCT_WT_hires", "HCT_RIF1KO_hires", "H1_WT_hires", "H9_WT_hires", "H9_RIF1KO_hires")) {
        RTindex = assemble_RTindex_df(get(m), kept_Gilbert, nfractions = 16)
        RT_to_bed(RTindex, paste0("../out/bed/RTindex_", m, ".bed"))
}
```

```R
RTdiff_HCT_RIF1KO_vs_WT = assemble_RTdiff_df(HCT_WT_hires, HCT_RIF1KO_hires, kept_Gilbert)
```

```R
RT_to_bedgraph(RTdiff_HCT_RIF1KO_vs_WT, paste0("../out/bedgraph/RTdiff_HCT_RIF1KO_vs_WT.bedgraph"), "RTdiff_HCT_RIF1KO_vs_WT")
```

```R
for (m in c("RTdiff_HCT_RIF1KO_vs_WT")) {
    RT_to_bed(get(m), paste0("../out/bed/", m, ".bed"))
    }
```

### RTindex for the Hansen et al. K562 data (6 fractions)s

```R
K562_smoothed_norm %>% head()
```

```R
RTindex_K562 = assemble_RTindex_df(K562_smoothed_norm, 
                                   rep(T, times = nrow(K562_smoothed_norm)), # no need to filter out anything, there is only one sample
                                   nfractions = 6)
```

```R
RTindex_K562 %>% head()
```

```R
# no filtering
RTindex_HCT_WT = assemble_RTindex_df(HCT_WT_hires, rep(T, times = nrow(HCT_WT_hires)))
```

```R
RTindex_HCT_WT %>% head()
```

### RTindex for the Trono HEK293T repliseq (4 fractions)

```R
hek293t_274ko_smoothed_norm %>% dim()
```

```R
hek293t_smoothed_norm %>% dim()
```

Careful: the rows are not the same, meaning that the underlying coordinates are not the same between WT and ZNF274KO

```R
# any rows with only zeros?
any(hek293t_smoothed_norm %>% dplyr::select(S1:S4) %>% apply(., 1, function(x) sum(x == 0) == 4))
any(hek293t_274ko_smoothed_norm %>% dplyr::select(S1:S4) %>% apply(., 1, function(x) sum(x == 0) == 4))
```

```R
hek293t_274ko_smoothed_norm %>% 
    dplyr::select(1:3) %>% head()
```

```R
kept_293T_WT = hek293t_smoothed_norm %>% 
    dplyr::select(1:3) %>%
    dplyr::mutate(idx = 1:nrow(hek293t_smoothed_norm)) %>%
    dplyr::inner_join(., y = hek293t_274ko_smoothed_norm %>%
                     dplyr::select(1:3)) %>%
    dplyr::pull(idx)
kept_293T_WT %>% head()
```

```R
hek293t_smoothed_norm %>% dim()
hek293t_smoothed_norm = hek293t_smoothed_norm[kept_293T_WT, ]
hek293t_smoothed_norm %>% dim()
```

```R
kept_293T_ZNF274KO = hek293t_274ko_smoothed_norm %>% 
    dplyr::select(1:3) %>%
    dplyr::mutate(idx = 1:nrow(hek293t_274ko_smoothed_norm)) %>%
    dplyr::inner_join(., y = hek293t_smoothed_norm %>%
                     dplyr::select(1:3)) %>%
    dplyr::pull(idx)
kept_293T_ZNF274KO %>% head()
```

```R
hek293t_274ko_smoothed_norm %>% dim()
hek293t_274ko_smoothed_norm = hek293t_274ko_smoothed_norm[kept_293T_ZNF274KO, ]
hek293t_274ko_smoothed_norm %>% dim()
```

```R
hek293t_274ko_smoothed_norm = hek293t_274ko_smoothed_norm %>% arrange(chrom, start, end)
hek293t_smoothed_norm = hek293t_smoothed_norm %>% arrange(chrom, start, end)
```

```R
stopifnot(identical(hek293t_274ko_smoothed_norm %>% dplyr::select(1:3),
                hek293t_smoothed_norm %>% dplyr::select(1:3)))
```

```R
RTindex_HEK293T = assemble_RTindex_df(hek293t_smoothed_norm, 
                                   rep(T, times = nrow(hek293t_smoothed_norm)), # no need to filter out anything, there is only one sample
                                   nfractions = 4)
RTindex_HEK293T_274KO = assemble_RTindex_df(hek293t_274ko_smoothed_norm, 
                                   rep(T, times = nrow(hek293t_274ko_smoothed_norm)), # no need to filter out anything, there is only one sample
                                   nfractions = 4)
```

```R
RTindex_HEK293T %>% head()
```

```R
RTindex_HEK293T_274KO %>% head()
```

```R
for (RT in c("RTindex_HEK293T", "RTindex_HEK293T_274KO")) {
    RT_to_bedgraph(get(RT), paste0("../out/bedgraph/", RT, ".bedgraph"), RT)
    }
```

```R
RTindex_HEK293T %>% dim()
```

```R
RTindex_HEK293T_274KO %>% dim()
```

```R
RT_diff_293T_ZNF274KO_vs_WT = cbind(RTindex_HEK293T[1:3], data.frame(RTdiff = RTindex_HEK293T_274KO$RTindex-RTindex_HEK293T$RTindex))
```

```R
RT_diff_293T_ZNF274KO_vs_WT %>% head()
```

```R
RT_to_bedgraph(RT_diff_293T_ZNF274KO_vs_WT, paste0("../out/bedgraph/", "RTdiff_ZNF274KO_vs_WT_293T", ".bedgraph"), name = "RTdiff_ZNF274KO_vs_WT_293T")

```

```R
for (m in c("RT_diff_293T_ZNF274KO_vs_WT")) {
    RT_to_bed(get(m), paste0("../out/bed/RTdiff", m, ".bed"))
    }
```

These bedfiles can be explored in a genome browser, on hg19


Regions of the genome with a negative RTDiff upon ZNF274 KO (removal of 274 -> earlier RT) and their approx. value: 
chr1: 
- DAB1:  -0.5
- ESRRG - RP11: -0.5
- ZNF678: -1
- ZNF695/ZNF670/ZNF124: -0.8

chr2: 
- SLC8A2: -1.1
- AC007092.1: -1.06
- CTNNA2: -0.78
- ZNF2/ZNF514: -1.35
- NCKAP5: -1.5

chr3: 
- ZNF445/ZKSCAN7/ZNF501/ZNF502 etc: -1
- FOXP1: -0.48
- ZNF717: -0.6
- LSAMP: -1.20
- EPHB1: -0.81
- CLSTN2: -0.48
- SOX2-OT: -0.55

chr4: 
- ZNF595/718/141/721 etc: -1.3
- chr4:22,836,110-25,396,211: -0.8
- CXCL6-CXCL2: -0.53
- PARM1: -0.55
- ANTXR2: -0.52
- FGF5: -0.58
- RASGEF1B: -0.67
- RP11-255I10.2: -0.7
- chr4:181,148,638-181,899,600: -0.58

chr5: 
- PCDHB15 / PCDH cluster: -1.7
- ZFP2/ZNF354A/ZNF354B/ZNF897etc : -2

chr6:
- chr6:8,294,151-9,312,251: - 0.8
- KRT18P38/ID4: -0.73
- ZNF204P/ZNF391/ZNF184: -0.5
- KCNK17/KIF6: -0.57
- FAM83B: -0.5
- TPBG: -0.6
- HS3ST5: - 0.73
- GRM1: -0.56
- PARK2: -0.86

chr7:
- ZNF316/ZNF12
- POU6F2: -0.51
- chr7:48,570,486-50,504,659: -0.52
- chr7:56,284,855-56,768,398 pseudogene-rich region: -1.24
- ZNF479: -1.7
- ZNF90P3/ZNF733P/ZNF734P -> ZNF736/ZNF107/ZNF138/etc: -2.7
- WBSCR17: -0.9
- U1: -0.51
- PLXNA4: -0.63

chr8: 
- ZNF596: -1
- ERICH1: -1
- RUNX1T1: -0.51

chr9: 
- FAM27C: -0.58
- FAM27B / FAM27E3: -0.5
- RORB: -0.53
- chr9:92,652,052-93,072,219: -0.65
- TMEM246/RNF20: -0.5
- ZNF483: -0.65
- ZNF883/ZFP37: -1.3

chr10:
- ADARB2: -0.49
- KLF6: -0.53
- chr10:25,218,314-25,610,093: -0.53
- ZNF37BP/ZNF33B: -0.89
- ZNF487/ZNF239/ZNF32 etc: -0.87
- CXCL12: -0.71
- PCDH15: -0.5

chr11:
- chr11:29,078,498-30,573,212: -0.72
- NOX4: -0.53
- chr11:94,826,618-95,670,408 (y-RNA, there are several with negative scores): -0.66
- ETS1: -0.54

chr12:
- KCNA1/KCNA5/KCNA6: -0.51
- RERG: -0.76
- TRHDE: -0.84
- chr12:77,294,175-78,865,738: -0.53
- chr12:96,867,848-97,115,625: -0.552
- ZNF605/ZNF26/ZNF84/ etc: -2.7

chr14:
- SEMA6D

chr13: 
- SERP2: -0.53
- chr13:78,330,732-80,022,288: -0.75

chr14:
- STXBP6: -0.811

chr15: 
- NRF2: -1.1
- SLC5A2/COX6A2/ARMC5: -0.68
- ZNF267: -0.63
- chr16:54,203,556-54,968,246: -0.596
- ZFP90: -0.615
- ADAMTS18: -1.14
- MAF: -0.98 (several things related to eye/cornea development on chr15?)
- chr17:67,850,074-71,797,465/SLC39A1: -0.61

chr18:
- chr18:14,229,973-14,457,968/ZNF519 (borderline): -0.61
- CELF4: -0.72
- RIT2: -1.2
- RAB27B: -0.82
- CCBE1: -0.52
- CDH20: -0.84
- chr18:61,225,156-63,123,060: -0.76
- CD226: -0.59
- CYB5A: -0.54
- ZNF236/MBP: -0.53

chr19: 
- ZNF317->ZNF562:-1.53
- ZNF833->ZNF443: -2.26
- ZNF486->ZNF254: -3.2
- TSHZ3: -0.53
- ZNF302/ZNF181/ZNF30/SCN1B: -0.86
- ZNF565->ZNF573: -3.4
- ZNF404->ZNF229: -2.16
- ZNF577->ZNF331: -3.8
- ZNF787->ZNF470: -2.3
- ZNF264->ZNF274: -2.975 ZNF274 at the limit, end of the cluster doesnt change -> mirrors RIF1 high, znf274 insensitivity

chr20: 
- chr20:4,177,368-4,695,076: -0.67
- FLRT3/RNF11B/MACROD2: -0.51
- chr20:50,229,409-53,293,461: -0.91

chr21:
- chr21:31,564,316-32,026,276: -0.56

chr22:
- chr22:48,049,630-50,340,838: -0.56

chrX:
- chrX:20,613,922-20,974,298: -0.67
- DMD: -0.6
- ZNF81/ZNF182/SSX6/SSX5: -0.59
- XIAP: -0.85
- SLITRK2: -0.5
- PASD1: -0.63

Regions of the genome with a positive RTDiff (removal of 274 -> later RT): 
chr1: 
- RP4-601K24.1: +0.6

chr2: 
- IGKV2D-40: +0.67

chr3: 
- COL6A5: 0.75
- RP11-208P4.1: 0.5
- DLG1: 0.6

chr4: 
- CWH43: 0.66
- GK2: 0.7
- PRDM5: 0.46
- NUDT6/SPATA5: 0.6
- PGRMC2: 0.53
- NPY2R, MAP9: 0.65
- chr4:160,313,512-161,115,676: 0.9
- SORBS2: 0.75

chr5:
- chr5:22,403,967-23,367,793: 0.88
- chr5:30,397,661-31,132,628: 0.55
- NNT: 0.55

chr6:
- chr6:62,896,946-63,515,709: 0.53
- KCNQ5: 0.52
- chr6:120,414,313-121,644,202: 1.16
- chr6:135,918,193-136,063,257: 0.6

chr7: 
- RBMX2P4
- ISPD: 0.5
- chr7:112,921,451-113,092,151: 0.67

chr8: 
- NRG1: 0.5
- HGSNAT: 1.06
- SNTG1: 0.5

chr9: 
- SLC25A5P8: 0.59
- MAMDC2: 0.54
- TLE1: 0.57
- OR13C9: 0.64

chr10:
- NSUN6: 0.51
- chr10:38,967,875-39,179,647: centromeric: 0.68
- HELLS: 0.55
- R3HCC1L: 0.51
- ATRNL1: 0.57

chr11:
- chr11:18,669,360-19,079,577: 0.56
- OR9Q1 olfactory receptor cluster: 0.47. There seems to be a general trend of weak positive score  for the olfactory receptors.
- chr11:79,571,017-80,643,835: 0.51
- chr11:96,895,117-97,979,990: 0.78
- chr11:106,080,438-106,622,874: 0.9

chr12:
- chr12:38,640,524-39,518,926: 0.53

chr13: 
- CENPJ/SLC25A15P3: 0.67

chr14:
- SLC25A21/MIPOL1: 0.67
- SYNE2: 0.5
- chr14:83,858,820-86,255,014: 0.84
- chr14:98,913,011-99,325,156: 0.54
- UNC13C: 1.15

chr15: 
- DET1/AEN: 0.56

chr16: 
- chr16:33,830,063-34,152,136: 0.61

chr17: 
- KSR1: 0.55

chr18:
- chr18:35,645,807-35,938,596: 0.94

chr19:
- ZNF541: 0.78

chr20: 
- EBF4/CPXM1: 0.81

chr21: 
- chr21:10,558,221-11,043,911: 0.54

chrX:
- chrX:20,941,868-21,302,244: 0.55
- SPIN4: 0.59
- FAM46D: 0.97
- HDX: 0.54
- chrX:98,572,982-99,224,563: 1.10
- SLC6A14: 0.93



```R
RTindex_HCT_WT = assemble_RTindex_df(HCT_WT_hires, kept_Gilbert)

```

```R
RT_to_bedgraph(RTindex_K562, paste0("../out/bedgraph/RTindex_K562.bedgraph"), "K562")

```

```R
HCT_RIF1KO_RTdiff = assemble_RTdiff_df(HCT_WT_hires, HCT_RIF1KO_hires, kept_Gilbert)
H9_RIF1KO_RTdiff = assemble_RTdiff_df(H9_WT_hires, H9_RIF1KO_hires, kept_Gilbert)
H9_vs_HCT_RTdiff = assemble_RTdiff_df(H9_WT_hires, HCT_WT_hires, kept_Gilbert)
```

```R
for (m in c("HCT_RIF1KO_RTdiff", "H9_RIF1KO_RTdiff", "H9_vs_HCT_RTdiff")) {
    RT_to_bedgraph(RTindex, paste0("../out/bedgraph/RTindex_", m, ".bedgraph"), m)
    }
```

At this stage, we could convert the files to bed and compute the intersect with gene bodies to get average replitiming for each gene, and also the change in RT

```R
RT_to_bed <- function(RT, path) {
    RT %>% dplyr::arrange(chrom, start, end) %>% write.table(path, col.names = F, row.names = F, sep = '\t', quote = F)
    }
dir.create('../out/bed/')
```

```R
for (m in c("HCT_WT_hires", "HCT_RIF1KO_hires", "H1_WT_hires", "H9_WT_hires", "H9_RIF1KO_hires")) {
        RTindex = assemble_RTindex_df(get(m), kept_Gilbert)
        RT_to_bed(RTindex, paste0("../out/bed/RTindex_", m, ".bed"))
}
for (m in c("HCT_RIF1KO_RTdiff", "H9_RIF1KO_RTdiff", "H9_vs_HCT_RTdiff")) {
    RT_to_bed(get(m), paste0("../out/bed/RTindex_", m, ".bed"))
    }
```

```R
RT_to_bed(RTindex_K562, "../out/bed/RTindex_K562.bed")
```

## Gaussian mixture model clustering for the RTdiff values

As described in Klein et al., Science 2021

```R
# plotting the distributions of RT index values
for (m in c("HCT_RIF1KO_RTdiff", "H9_RIF1KO_RTdiff", "H9_vs_HCT_RTdiff")) {
    pdf(paste0("../out/cell_cycle_figures/", m, ".pdf"), width = 5, height = 5)
    hist(na.omit(get(m)$RTdiff), breaks = 1000, xlab = m, main = NULL)
    dev.off()
}
```

The distribution of changes in values motivates the classification. 

```R
set.seed(5)
fit = mclust::Mclust(na.omit(HCT_RIF1KO_RTdiff$RTdiff), G=3, model="V")
summary(fit)
```

```R
fit$parameters
```

```R
pdf("../out/cell_cycle_figures/RTdiff_KO_WT_Gaussian_mixture.pdf", width = 5, height = 5)
plot(fit, what="density", main="", xlab="RTdiff")
dev.off()
```

```R
HCT_RIF1KO_RTdiff$cluster = fit$classification
```

```R
HCT_RIF1KO_RTdiff %>% head()
```

```R
# adding the WT values to discriminate between late unchanged and early unchanged
HCT_RIF1KO_RTdiff$RTindex_wt = RTindex_HCT_WT$RTindex
```

```R
HCT_RIF1KO_RTdiff %>% head()
```

```R
clusters_to_cat = c("EtL", "U", "LtE")
HCT_RIF1KO_RTdiff$category = clusters_to_cat[HCT_RIF1KO_RTdiff$cluster]
HCT_RIF1KO_RTdiff %>% head()
```

```R
LU_rows = (HCT_RIF1KO_RTdiff$category == "U") & (HCT_RIF1KO_RTdiff$RTindex_wt <0)
LU_rows %>% head()
```

```R
EU_rows = (HCT_RIF1KO_RTdiff$category == "U") & (HCT_RIF1KO_RTdiff$RTindex_wt >= 0)

```

```R
HCT_RIF1KO_RTdiff[LU_rows, "category"] = "LU"
HCT_RIF1KO_RTdiff[EU_rows, "category"] = "EU"
```

```R
HCT_RIF1KO_RTdiff %>% head()
```

```R
HCT_RIF1KO_RTdiff$category %>% head()
HCT_RIF1KO_RTdiff$category %>% table()
```

### Exporting bedgraphs and bedfiles for the four types of bins

```R
for (cat in (HCT_RIF1KO_RTdiff$category %>% unique())) {
    track_def_line = paste0("track type=bedGraph name=RIF1KO_", cat, "_HCT description=Klein_et_al_processed_by_Pulver visibility=display_mode color=r,g,b altColor=r,g,b priority=priority autoScale=off")
    rbind(c(track_def_line, rep("", times = 6)), (HCT_RIF1KO_RTdiff %>% 
                                                    dplyr::filter(category == cat) %>%
                                                 dplyr::select(chrom, start, end, RTdiff, cluster, RTindex_wt, category)))  %>% write.table(., paste0('../out/bedgraph/RIF1KO_', cat, '_rif1kovswt.bedgraph'), quote=F, sep = '\t', row.names = F, col.names = F)
    }

```

```R
HCT_RIF1KO_RTdiff %>% dplyr::arrange(chrom, start, end) %>% dplyr::select(1:7) %>% write.table('../out/rif1ko_hct_RTbins.bed', col.names = F, row.names = F, sep = '\t', quote = F)
```

### Among LU regions, which are H3K9me3-dependent?

```R
# importing the log2EL of the H3K9me3 knockdown
h3k9me3_kd_ctrl_log2EL = read.table('../data/repliseq_rif1/GSE160563_HCT_RIF1KO_H3K9me3KD_ctrlsiRNA_log2EL_RT.bedgraph', header = F, sep = '\t')
colnames(h3k9me3_kd_ctrl_log2EL) = c("chrom", "start", "end", "RIF1KO_ctrl_log2EL")
h3k9me3_kd_ctrl_log2EL %>% head()
h3k9me3_kd_ctrl_log2EL %>% dim()
```

```R
h3k9me3_kd_siRNA_log2EL = read.table('../data/repliseq_rif1/GSE160563_HCT_RIF1KO_H3K9me3KD_KDsiRNA_log2EL_RT.bedgraph', header = F, sep = '\t')
colnames(h3k9me3_kd_siRNA_log2EL) = c("chrom", "start", "end", "RIF1KO_H3K9me3KD_log2EL")
h3k9me3_kd_siRNA_log2EL %>% head()
h3k9me3_kd_siRNA_log2EL %>% dim()
```

```R
# merging into the RTdiff bin table
HCT_RIF1KO_RTdiff = HCT_RIF1KO_RTdiff %>% dplyr::left_join(., h3k9me3_kd_ctrl_log2EL)
HCT_RIF1KO_RTdiff = HCT_RIF1KO_RTdiff %>% dplyr::left_join(., h3k9me3_kd_siRNA_log2EL)

```

```R
HCT_RIF1KO_RTdiff %>% head()
```

```R
HCT_RIF1KO_RTdiff$category %>% unique()
```

```R
HCT_RIF1KO_RTdiff$RIF1KO_h3k9me3kd_vs_ctrl_log2EL = HCT_RIF1KO_RTdiff$RIF1KO_H3K9me3KD_log2EL-HCT_RIF1KO_RTdiff$RIF1KO_ctrl_log2EL
```

```R
HCT_RIF1KO_RTdiff$RIF1KO_h3k9me3kd_vs_ctrl_log2EL %>% summary()
```

```R
HCT_RIF1KO_RTdiff %>% dplyr::filter(!is.na(RIF1KO_ctrl_log2EL)) %>% ggplot(., aes(x = category, y = RIF1KO_h3k9me3kd_vs_ctrl_log2EL)) + geom_violin()
```

Late unchanged regions are shifted to earlier upon H3K9me3 KD, i.e to a more positive log2EL. Thus, the deltalog2EL is expected to be positive for LU regions compared to other regions. This is what we observe. We've also seen that LtE regions tend to behave a bit similarly, particularly if they're close to the LU threshold. Thus, does it mean that LtE regions that become even more early indicate a cumulative effect of RIF1 + H3K9me3? Are these bound by ZNF274, or other KZFPs? What kind of factors enrich there? Alternatively, these may just be Late Unchanged regions that were wrongly assigned to Late to Early. We could look at the distribution of LtE with a large positive change between H3K9me3 kd and ctrl in terms of delta RT index to get a sense of what's going on, I suspect they may be borderline Late Unchanged upon RIF1KO.

As for the early to late and early unchanged, they become a bit more negative, meaning a bit later. Given the competition model, this is simply explained by the buffering of replication origin factors by the newly liberated LU and LtE regions that have lost H3K9me3.

There seems to be a lot of L1 enriched at LtE regions in general. I think that corroborates an observation made by Gilbert in a Genome Research a long time ago, the one with isochores.


### Defining a logE/L threshold to call a region H3K9me3-dependent

```R
# there is a KZFP cluster from ~21 to ~25mb, and it is H3K9me3-dependent.
HCT_RIF1KO_RTdiff %>% dplyr::filter(chrom == "chr19", start > 17710000, end < 26000000) %>% 
    dplyr::select(chrom, start, end, category, RIF1KO_h3k9me3kd_vs_ctrl_log2EL) %>% 
    dplyr::filter(!is.na(RIF1KO_h3k9me3kd_vs_ctrl_log2EL)) %>% 
    ggplot(aes(x = start, y = RIF1KO_h3k9me3kd_vs_ctrl_log2EL)) + geom_point()
```

```R
# there is another KZFP cluster from 36.7Mb to 38Mb
HCT_RIF1KO_RTdiff %>% dplyr::filter(chrom == "chr19", start > 35000000, end < 39000000) %>% 
    dplyr::select(chrom, start, end, category, RIF1KO_h3k9me3kd_vs_ctrl_log2EL) %>% 
    dplyr::filter(!is.na(RIF1KO_h3k9me3kd_vs_ctrl_log2EL)) %>% 
    ggplot(aes(x = start, y = RIF1KO_h3k9me3kd_vs_ctrl_log2EL)) + geom_point()
```

Just taking delta log2EL > 0.5 or > 1 seems like a reasonable way to assign LU regions into affected and non-affected.


Most KZFPs seem to stand in unchanged regions, predominantly late replicating. But there are also a lot in early replicating, unchanged. The major exceptions are clusters 7.2, 7.3 and 8.1, which are in EtL regions. Interestingly, these clusters are ZNF274 low. 

Maybe we should group early unchanged and late unchanged in "unchanged" for some of the comparisons involving KZFP clusters. Would be interesting to compare non-KRAB ZNFs to these as well, together with other TFs. 


```R
HCT_RIF1KO_RTdiff %>% head()
```

```R
# detecting category changes
HCT_RIF1KO_RTdiff = HCT_RIF1KO_RTdiff %>% dplyr::mutate(category_change = c(1, diff(as.integer(as.factor(category)))))
```

```R
# detecting chromosome changes
HCT_RIF1KO_RTdiff = HCT_RIF1KO_RTdiff %>% dplyr::mutate(chr_change = c(1, diff(as.integer(as.factor(chrom)))))
```

```R
# detecting discontinuities within chromosomes: if the distance to the next window is more than 50kb
HCT_RIF1KO_RTdiff = HCT_RIF1KO_RTdiff %>% dplyr::mutate(bin_break = c(1, diff(start)>50000))
```

```R
HCT_RIF1KO_RTdiff %>% head(10)
```

```R
breaks_regions = which(HCT_RIF1KO_RTdiff$category_change!=0)
breaks_regions %>% length()
```

```R
breaks_regions = which(HCT_RIF1KO_RTdiff$category_change!=0 | HCT_RIF1KO_RTdiff$chr_change!=0)
breaks_regions %>% length()
```

```R
breaks_regions = which(HCT_RIF1KO_RTdiff$category_change!=0 | HCT_RIF1KO_RTdiff$chr_change!=0 | HCT_RIF1KO_RTdiff$bin_break)
breaks_regions %>% length()
```

Indeed breaking at chromosomes too is required, and at discontinuities too

```R
HCT_RIF1KO_RTdiff[45:55, ]
```

We may need to add a break for the last line.

```R
regions = list()
for(i in 1:(length(breaks_regions))) {
    region_long = NULL
    if(i == 1) {
    region_long = HCT_RIF1KO_RTdiff[breaks_regions[i]:breaks_regions[i+1]-1, ]
        } else if(i>1 & i < length(breaks_regions)) {
        region_long = HCT_RIF1KO_RTdiff[(breaks_regions[i]+1):breaks_regions[i+1]-1, ]

        } else if (i == length(breaks_regions)) {
        region_long = HCT_RIF1KO_RTdiff[(breaks_regions[i]+1):nrow(HCT_RIF1KO_RTdiff), ]
        }
    region_row = list(chrom=region_long$chrom %>% unique(),
                  start = (region_long %>% head(1))$start,
                  end = (region_long %>% tail(1))$end,
                  category = region_long$category %>% unique(),
                  RTindex_wt_avg = mean(region_long$RTindex_wt),
                  HCT_RIF1KO_RTdiff = mean(region_long$RTdiff))
    regions[[i]] = region_row %>% as.data.frame()
    }
```

```R
regions_df = do.call(rbind, regions)
```

```R
regions_df %>% head()
regions_df %>% dim()
```

```R
regions_df %>% tail()
```

```R
HCT_RIF1KO_RTdiff %>% tail(10)
```

```R
# exporting bedfiles for the regions
dir.create("../out/rif1ko_hct_RTregions/hg38/", showWarnings = F)
for (c in regions_df$category %>% unique()) {
    regions_df %>% dplyr::filter(category == c) %>% 
        dplyr::arrange(chrom, start, end) %>% 
        write.table(paste0('../out/rif1ko_hct_RTregions/hg38/rif1ko_hct_', c, '_regions.bed'), row.names = F, col.names = F, sep = '\t', quote = F)
    }

```

```R
# as a single bed file
regions_df %>% dplyr::arrange(chrom, start, end) %>% write.table(paste0('../out/rif1ko_hct_RTregions/rif1ko_hct_regions_hg38.bed'), row.names = F, col.names = F, sep = '\t', quote = F)
```

```R
regions_df$category %>% table()
```

Our annotation of the regions is much more fragmented that that from Klein's et al. paper. It may be because:
1) they're not using the same control, which remains unspecified in their paper
2) they smoothed out the control sample even more using a gaussian kernel
3) There are some unspecified parameters in their score, in particular for the gaussian mixture model clustering.


## Enrichment of ENCODE TF and Trono lab KZFP binding sites (peaks) in late unchanged, H3K9me3-dependent regions

The script `intersect_RTregions_with_TFs_and_gene_bodies.sh` must be ran before.

```R
encode_tfs_metadata = read.table('../data/ENCODE_k562_tf_peaks/experiment_report_2024_2_19_10h_10m.tsv', skip = 1, header = T, sep = '\t')
```

```R
i = 0
df_list = list()
for (bed in list.files('../out/encode_x_rif1ko_hct_RTregions/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    p = paste0('../out/encode_x_rif1ko_hct_RTregions/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        temp_df = read.table(paste0('../out/encode_x_rif1ko_hct_RTregions/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        temp_df$TF_symbol = encode_tfs_metadata[grepl(filename, encode_tfs_metadata$Files) %>% which(), "Target.gene.symbol"]
        temp_df$peaks_total = R.utils::countLines(paste0('../data/ENCODE_k562_tf_peaks/bed_filtered_hg19/', bed))
        i = i+1

        df_list[[i]] = temp_df}
    }
```

```R
peaks_x_RTbins_encode = do.call(rbind, df_list)
```

```R
peaks_x_RTbins_encode$cell_line = "K562"
peaks_x_RTbins_encode$publication = "ENCODE"
```

```R
peaks_x_RTbins_encode %>% head() %>% print()
```

```R
# identifying the overlap column: 
min(897329, 896245) - max(895467, 895935)
```

```R
peaks_x_RTbins_encode$V5 = NULL
peaks_x_RTbins_encode$V11 = NULL
peaks_x_RTbins_encode$V12 = NULL
peaks_x_RTbins_encode$V13 = NULL
peaks_x_RTbins_encode$V14 = NULL
peaks_x_RTbins_encode$V15 = NULL
peaks_x_RTbins_encode$V16 = NULL
peaks_x_RTbins_encode$V17 = NULL # pointsource peak (where to call the peak if a single base is allowed, from the start coordinate)
peaks_x_RTbins_encode %>% dim()
```

```R
colnames(peaks_x_RTbins_encode) = c("RTregion_chr", "RTregion_start", "RTregion_end", "RTregion_RTdiff_rif1kovswt_array_avg", "RTregion_RTindex_wt_avg", "RTregion_category",
                                      "peak_chr", "peak_start", "peak_end", "peak_overlap", "filename", "TF_symbol", "peaks_total", "cell_line", "publication")
peaks_x_RTbins_encode %>% head()
peaks_x_RTbins_encode %>% dim()
```

```R
peaks_tronolab = read.table('../data/hg19_peaks_un_filt_macs80.bed', sep = '\t', header = F, )
peaks_tronolab$V6 = NULL
```

```R
colnames(peaks_tronolab) = c("peak_chr", "peak_start", "peak_end", "peak_ID", "peak_score")
```

```R
# extracting KZFP names
peaks_tronolab$TF_symbol = sapply(peaks_tronolab$peak_ID, function(x) return((str_split(x, "_", simplify = T))[1])) %>% as.vector()
```

```R
peaks_tronolab %>% head()
```

```R
peaks_tronolab$filename = sapply((peaks_tronolab$peak_ID), function(x) return(paste0((str_split(x, "_", simplify = T))[1:5], collapse = '_'))) %>% as.vector()
```

```R
# extracting file names (to be able to compare with ENCODE file names and count peaks per file and thus per TF at the end)
peaks_tronolab$TF_symbol = sapply(peaks_tronolab$peak_ID, function(x) return((str_split(x, "_", simplify = T))[1])) %>% as.vector()
```

```R
peaks_tronolab$TF_symbol %>% unique() %>% length()
```

```R
# counting all peaks per KZFP
peaks_tronolab = peaks_tronolab %>% dplyr::group_by(TF_symbol) %>% dplyr::mutate(TF_peaks_total=n()) %>% ungroup()
```

```R
peaks_tronolab %>% head()
```

```R
peaks_tronolab %>% dim()
```

```R
# annotating low confidence ChIP-seq
peaks_tronolab$peak_lowconf = grepl("lowconf", peaks_tronolab$peak_ID)
```

```R
peaks_x_RTbins_tronolab = read.table('../out/rif1ko_hct_RTbins_vs_kzfp_peaks.bed', sep = '\t', header = F)
```

```R
peaks_x_RTbins_tronolab %>% head()
```

```R
peaks_x_RTbins_tronolab %>% dim()
```

```R
peaks_x_RTbins_tronolab$V5 = NULL
peaks_x_RTbins_tronolab$V12 = NULL
peaks_x_RTbins_tronolab$V13 = NULL
```

```R
colnames(peaks_x_RTbins_tronolab) = c("RTregion_chr", "RTregion_start", "RTregion_end", "RTregion_RTdiff_rif1kovswt_array_avg", "RTregion_RTindex_wt_avg", "RTregion_category", 
                                         "peak_chr", "peak_start", "peak_end", "peak_ID", "peak_overlap")
```

```R
peaks_x_RTbins_tronolab %>% head()
```

```R
# adding the peak information
peaks_x_RTbins_tronolab = peaks_x_RTbins_tronolab %>% dplyr::left_join(., peaks_tronolab %>% dplyr::select(peak_ID, TF_symbol, peak_lowconf, TF_peaks_total, filename))
```

```R
peaks_x_RTbins_tronolab %>% head()
```

```R
peaks_x_RTbins_tronolab$cell_line = "HEK293T"
peaks_x_RTbins_tronolab$publication = "Tronolab"
```

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
for (bed in list.files('../out/schmitges_x_rif1ko_hct_RTregions/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    #p = paste0('../out/schmitges_x_promoters/', bed)
    p = paste0('../out/schmitges_x_rif1ko_hct_RTregions/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        #temp_df = read.table(paste0('../out/schmitges_x_promoters/', bed), sep = '\t', header = F)
        temp_df = read.table(paste0('../out/schmitges_x_rif1ko_hct_RTregions/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        #temp_df$TF_symbol = schmitges_metadata %>% can be done afterwards
        #temp_df$peaks_total = R.utils::countLines(paste0('../data/schmitges_kzfp_chips/hg19_peaks_chipatlas/', bed))
        temp_df$peaks_total = R.utils::countLines(paste0('../data/schmitges_all_znfs_chips/hg19_peaks_chipatlas/', bed))
        i = i+1

        df_list[[i]] = temp_df}
    }
```

```R
peaks_x_RTbins_schmitges = do.call(rbind, df_list)
```

```R
peaks_x_RTbins_schmitges %>% head() %>% print()
```

```R
peaks_x_RTbins_schmitges$cell_line = "HEK293T"
peaks_x_RTbins_schmitges$publication = "schmitges"
peaks_x_RTbins_schmitges = peaks_x_RTbins_schmitges %>% dplyr::left_join(., schmitges_metadata %>% dplyr::select(filename, TF_symbol))
peaks_x_RTbins_schmitges %>% head() %>% print()
```

```R
peaks_x_RTbins_schmitges$V5 = NULL
peaks_x_RTbins_schmitges$V12 = NULL
peaks_x_RTbins_schmitges$V13 = NULL
peaks_x_RTbins_schmitges$V14 = NULL
peaks_x_RTbins_schmitges$V15 = NULL
peaks_x_RTbins_schmitges$V16 = NULL
peaks_x_RTbins_schmitges$V17 = NULL # pointsource peak (where to call the peak if a single base is allowed, from the start coordinate)
peaks_x_RTbins_schmitges %>% dim()
```

```R
peaks_x_RTbins_schmitges %>% head() %>% print()
```

```R
colnames(peaks_x_RTbins_schmitges) = c("RTregion_chr", "RTregion_start", "RTregion_end",
                                       "RTregion_RTdiff_rif1kovswt_array_avg", "RTregion_RTindex_wt_avg", "RTregion_category", 
                                         "peak_chr", "peak_start", "peak_end", "peak_ID", "peak_overlap", 
                                       "filename", "peaks_total", "cell_line", "publication", "TF_symbol")
```

```R
peaks_x_RTbins_schmitges %>% head()
peaks_x_RTbins_schmitges %>% dim()
```

```R
# adding the najafabadi peaks at promoters
i = 0
df_list = list()
#for (bed in list.files('../out/najafabadi_x_promoters/', pattern = "*.bed")) {
for (bed in list.files('../out/najafabadi_x_rif1ko_hct_RTregions/', pattern = "*.bed")) {
    filename = (bed %>% str_split(., "\\."))[[1]][[1]]
    
    #p = paste0('../out/najafabadi_x_promoters/', bed)
    p = paste0('../out/najafabadi_x_rif1ko_hct_RTregions/', bed)
    file_info <- file.info(p)
    
    if (file_info$size != 0) {
        #temp_df = read.table(paste0('../out/najafabadi_x_promoters/', bed), sep = '\t', header = F)
        temp_df = read.table(paste0('../out/najafabadi_x_rif1ko_hct_RTregions/', bed), sep = '\t', header = F)
        temp_df$filename = filename
        temp_df$sample = str_split(filename, "_")[[1]][[1]]
        #temp_df$peaks_total = R.utils::countLines(paste0('../data/najafabadi_kzfp_chips/', bed))
        temp_df$peaks_total = R.utils::countLines(paste0('../data/najafabadi_all_znfs_chips/hg19_peaks_80/', bed))
        i = i+1

        df_list[[i]] = temp_df}
    }
```

```R
peaks_x_RTbins_najafabadi = do.call(rbind, df_list)
```

```R
peaks_x_RTbins_najafabadi %>% head()
```

```R
peaks_x_RTbins_najafabadi$cell_line = "HEK293T"
peaks_x_RTbins_najafabadi$publication = "najafabadi"
peaks_x_RTbins_najafabadi %>% head() %>% print()
```

```R
peaks_x_RTbins_najafabadi$V5 = NULL
peaks_x_RTbins_najafabadi$V12 = NULL
peaks_x_RTbins_najafabadi %>% head()
```

```R
najafabadi_metadata = read.table('../data/najafabadi_all_znfs_chips/najafabadi_metadata.txt', sep = ' ', header = F) %>% 
    dplyr::rename(symbol = V1, sample = V2)
najafabadi_metadata %>% head()
```

```R
peaks_x_RTbins_najafabadi = peaks_x_RTbins_najafabadi %>% 
    dplyr::left_join(., najafabadi_metadata)
```

```R
peaks_x_RTbins_najafabadi %>% head()
```

```R
colnames(peaks_x_RTbins_najafabadi) = c("RTregion_chr", "RTregion_start", "RTregion_end",
                                       "RTregion_RTdiff_rif1kovswt_array_avg", "RTregion_RTindex_wt_avg", "RTregion_category", 
                                         "peak_chr", "peak_start", "peak_end", "peak_ID", "peak_overlap", 
                                       "filename", "sample", "peaks_total", "cell_line", "publication", "TF_symbol")
peaks_x_RTbins_najafabadi %>% head()
peaks_x_RTbins_najafabadi %>% dim()
```

```R
# concatenating with the encode data
peaks_x_RTbins_tronolab$cell_line = "HEK293T"
peaks_x_RTbins = rbind(peaks_x_RTbins_encode %>% dplyr::select("RTregion_chr", "RTregion_start", "RTregion_end", "RTregion_category",  "RTregion_RTindex_wt_avg", "RTregion_RTdiff_rif1kovswt_array_avg", 
                                             "peak_chr", "peak_start", "peak_end", "peak_overlap", "filename", "TF_symbol", "cell_line", "peaks_total"),
      peaks_x_RTbins_tronolab %>% dplyr::rename(peaks_total = TF_peaks_total) %>% dplyr::select("RTregion_chr", "RTregion_start", "RTregion_end", "RTregion_category",  "RTregion_RTindex_wt_avg", "RTregion_RTdiff_rif1kovswt_array_avg", 
                                             "peak_chr", "peak_start", "peak_end", "peak_overlap", "filename", "TF_symbol", "cell_line", "peaks_total"),
                      peaks_x_RTbins_schmitges %>% dplyr::select("RTregion_chr", "RTregion_start", "RTregion_end", "RTregion_category",  "RTregion_RTindex_wt_avg", "RTregion_RTdiff_rif1kovswt_array_avg", 
                                             "peak_chr", "peak_start", "peak_end", "peak_overlap", "filename", "TF_symbol", "cell_line", "peaks_total"),
                      peaks_x_RTbins_najafabadi %>% dplyr::select("RTregion_chr", "RTregion_start", "RTregion_end", "RTregion_category",  "RTregion_RTindex_wt_avg", "RTregion_RTdiff_rif1kovswt_array_avg", 
                                             "peak_chr", "peak_start", "peak_end", "peak_overlap", "filename", "TF_symbol", "cell_line", "peaks_total"))
```

```R
# what are TF symbols without associated gene?
# RBM14,RBM14-RBM4 -> RBM14
# U2AF1L5,U2AF1  -> U2AF1
# COMMD3-BMI1,BMI1 -> BMI1
# ZNF724P -> ZNF724

peaks_x_RTbins$TF_symbol[which(peaks_x_RTbins$TF_symbol == "RBM14,RBM14-RBM4")] = "RBM4"
peaks_x_RTbins$TF_symbol[which(peaks_x_RTbins$TF_symbol == "U2AF1L5,U2AF1")] = "U2AF1"
peaks_x_RTbins$TF_symbol[which(peaks_x_RTbins$TF_symbol == "COMMD3-BMI1,BMI1")] = "BMI1"
peaks_x_RTbins$TF_symbol[which(peaks_x_RTbins$TF_symbol == "ZNF724P")] = "ZNF724"
```

```R
peaks_x_RTbins %>% dim()
peaks_x_RTbins$TF_symbol %>% unique() %>% length()
```

```R
# computing RTregion and peak length
peaks_x_RTbins = peaks_x_RTbins %>% dplyr::mutate(RTregion_length = RTregion_end-RTregion_start, peak_length = peak_end-peak_start)
```

```R
# splitting LU in three, with LU_K9KdLtE, LU_K9KdLtLr (late to later upon K9 knockdown) and LU (RT unchanged by K9 knockdown)
# Doing the same for LtE, i.e. LtE_K9KdLtE, LtE_K9KdLtLr and LtE.
```

```R
peaks_x_RTbins %>% head()
```

```R
HCT_RIF1KO_RTdiff %>% head()
```

```R
# adding the delta EtL
peaks_x_RTbins = peaks_x_RTbins %>% dplyr::left_join(., (HCT_RIF1KO_RTdiff %>% dplyr::select(chrom, start, end, category, RIF1KO_h3k9me3kd_vs_ctrl_log2EL)), by = c("RTregion_chr" = "chrom",
"RTregion_start" = "start", "RTregion_end" = "end", "RTregion_category" = "category"))

```

```R
peaks_x_RTbins$RTregion_category_k9 = peaks_x_RTbins$RTregion_category
idx_LU_K9KdLtE = (peaks_x_RTbins$RTregion_category == "LU") & (!is.na(peaks_x_RTbins$RIF1KO_h3k9me3kd_vs_ctrl_log2EL) & (peaks_x_RTbins$RIF1KO_h3k9me3kd_vs_ctrl_log2EL > 1))
idx_LU_K9KdLtLr = (peaks_x_RTbins$RTregion_category == "LU") & (!is.na(peaks_x_RTbins$RIF1KO_h3k9me3kd_vs_ctrl_log2EL) & (peaks_x_RTbins$RIF1KO_h3k9me3kd_vs_ctrl_log2EL < -1))

idx_LtE_K9KdLtE = (peaks_x_RTbins$RTregion_category == "LtE") & (!is.na(peaks_x_RTbins$RIF1KO_h3k9me3kd_vs_ctrl_log2EL) & (peaks_x_RTbins$RIF1KO_h3k9me3kd_vs_ctrl_log2EL > 1))
idx_LtE_K9KdLtLr = (peaks_x_RTbins$RTregion_category == "LtE") & (!is.na(peaks_x_RTbins$RIF1KO_h3k9me3kd_vs_ctrl_log2EL) & (peaks_x_RTbins$RIF1KO_h3k9me3kd_vs_ctrl_log2EL < -1))
```

```R
idx_LU_K9KdLtE %>% sum()
idx_LU_K9KdLtLr %>% sum()

idx_LtE_K9KdLtE %>% sum()
idx_LtE_K9KdLtLr %>% sum()
```

```R
peaks_x_RTbins[idx_LU_K9KdLtE, "RTregion_category_k9"] = "LU_K9KdLtE"
peaks_x_RTbins[idx_LU_K9KdLtLr, "RTregion_category_k9"] = "LU_K9KdLtLr"
peaks_x_RTbins[idx_LtE_K9KdLtE, "RTregion_category_k9"] = "LtE_K9KdLtE"
peaks_x_RTbins[idx_LtE_K9KdLtLr, "RTregion_category_k9"] = "LtE_K9KdLtLr"

```

```R
# adding the gene rhythmicity information for TFs only
rhythm_all_genes = read.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', header = 1)
rhythm_all_genes %>% head()
```

```R
# adding the TF cycling information
peaks_x_RTbins = peaks_x_RTbins %>% 
    dplyr::left_join(x=., y=rhythm_all_genes, by = c("TF_symbol" = "symbol"))
peaks_x_RTbins %>% colnames()
```

```R
colnames(peaks_x_RTbins)[((peaks_x_RTbins %>% ncol())-(rhythm_all_genes %>% ncol())+2) : (peaks_x_RTbins %>% ncol())] = paste0("TF_", colnames(peaks_x_RTbins)[((peaks_x_RTbins %>% ncol())-(rhythm_all_genes %>% ncol())+2) : (peaks_x_RTbins %>% ncol())])
```

```R
peaks_x_RTbins %>% colnames()
peaks_x_RTbins %>% head()
```

```R
# loading KZFP information
kzfps_cycling_info = read.csv('../out/tables/kzfps_cycling_info.tsv', header = T, sep = '\t')
kzfps_cycling_info %>% head()
```

## Binding profile of TFs w.r.t. RT regions

For each TF, we compute:
- Total number of peaks
- Total number of peaks in any RTregion
- Total number of peaks in late unchanged region
    - May be done for each of the four regions: `EtL`, `EU`, `LU`, `LtE`

```R
peaks_x_RTbins$RTregion_category_k9 %>% unique()
```

```R
(peaks_x_RTbins$RTregion_category_k9 == "LtE_K9KdLtLr") %>% sum()
```

```R
ChIP_centric_metrics = peaks_x_RTbins %>% dplyr::group_by(TF_symbol) %>% dplyr::summarize(#TF_peaks_total = sum(unique(peaks_total)), #TODO make this less subject to error if two numbers are the same
                                                                      TF_peaks_in_RTregions = n(), 
                                                                    TF_RTregions_bound = length(unique(RTregion_category)),
                                                                     TF_peaks_in_EtL = sum(sum(RTregion_category_k9 == "EtL"), na.rm = T),
                                                                     TF_peaks_in_EU = sum(sum(RTregion_category_k9 == "EU"), na.rm = T),
                                                                      TF_peaks_in_LU = sum(sum(RTregion_category_k9 == "LU"), na.rm = T),
                                                                    TF_peaks_in_LU = sum(sum(RTregion_category_k9 == "LU"), na.rm = T),
                                                                        TF_peaks_in_LU_K9KdLtE = sum(sum(RTregion_category_k9 == "LU_K9KdLtE"), na.rm = T),
                                                                            TF_peaks_in_LU_K9KdLtLr = sum(sum(RTregion_category_k9 == "LU_K9KdLtLr"), na.rm = T),
                                                                            TF_peaks_in_LtE_K9KdLtE = sum(sum(RTregion_category_k9 == "LtE_K9KdLtE"), na.rm = T),
                                                                                TF_peaks_in_LtE_K9KdLtLr = sum(sum(RTregion_category_k9 == "LtE_K9KdLtLr"), na.rm = T),
    TF_peaks_in_LtE = sum(sum(RTregion_category_k9 == "LtE"), na.rm = T)
)
```

```R
# adding the sum of all peaks per TF
peaks_total_df = peaks_x_RTbins %>% dplyr::select(TF_symbol, filename, peaks_total) %>% unique() %>% dplyr::group_by(TF_symbol) %>% dplyr::summarize(TF_peaks_total = sum(peaks_total))
peaks_total_df %>% head()
```

```R
ChIP_centric_metrics = ChIP_centric_metrics %>% dplyr::left_join(., peaks_total_df)
```

```R
ChIP_centric_metrics %>% head()
```

Enrichment tests: are any TFs enriched in binding to LU RTregions? 

For each TF:


| | TF peaks overlapping LU RTregions | TF peaks overlapping non-LU RTregions| total |
|---|:---:|:---:|:---:|
| peaks belonging to TF | x | m-x | m |
| peaks belonging to other TFs | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

Or the reverse:

| | RTregions with overlapping TF peak | RTregions without overlapping TF peak| total |
|---|:---:|:---:|:---:|
| RTregions of category (e.g. LU)| x | m-x | m |
| other RTregions | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
ChIP_centric_metrics$TF_peaks_in_LU_K9kdLtE %>% sum()
ChIP_centric_metrics$TF_peaks_total %>% sum()
```

Trying to invert it: 
| | Peaks belonging to TF | peaks belonging to other TFs| total |
|---|:---:|:---:|:---:|
| TF peaks overlapping LU RTregions | x | m-x | m |
| TF peaks not overlapping LU RTregions | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |


```R
enrichment_test <- function(row, RTregion_category) {
    x_observed = as.integer(row[[paste0("TF_peaks_in_", RTregion_category)]])
    m = n_peaks_in_cat
    n = n_peaks_total - m
    k = as.integer(row[["TF_peaks_total"]])
    
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
peaks_x_RTbins$RTregion_category_k9 %>% unique()
```

```R
c
```

```R
ChIP_centric_metrics %>% colnames()
```

```R
# computing constants 
for (c in (peaks_x_RTbins$RTregion_category_k9 %>% unique())[5:length((peaks_x_RTbins$RTregion_category_k9 %>% unique()))]) {
    print(c)
    
    n_peaks_in_cat = ChIP_centric_metrics[, paste0("TF_peaks_in_", c)] %>% sum()
    n_peaks_total = ChIP_centric_metrics$TF_peaks_total %>% sum()
    ratio_expected = n_peaks_in_cat/n_peaks_total
    
    RegionCat_vs_otherRegions_enrichment = ChIP_centric_metrics %>% apply(., 1, enrichment_test, RTregion_category = c) 

    names(RegionCat_vs_otherRegions_enrichment) = ChIP_centric_metrics$TF_symbol
    RegionCat_vs_otherRegions_enrichment = as.data.frame(RegionCat_vs_otherRegions_enrichment) %>% tibble::rownames_to_column("TF_symbol")
    colnames(RegionCat_vs_otherRegions_enrichment)[2] = 'p_enrich'
    RegionCat_vs_otherRegions_enrichment$is_kzfp = ifelse(RegionCat_vs_otherRegions_enrichment$TF_symbol %in% kzfps_cycling_info$symbol, T, F)

    RegionCat_vs_otherRegions_enrichment$padj_enrich = p.adjust(RegionCat_vs_otherRegions_enrichment$p_enrich, method = "BH")

    RegionCat_vs_otherRegions_enrichment$ratio_observed = (ChIP_centric_metrics[, paste0("TF_peaks_in_", c)]/ChIP_centric_metrics$TF_peaks_total)[, 1]
    RegionCat_vs_otherRegions_enrichment$ratio_expected = ratio_expected
    RegionCat_vs_otherRegions_enrichment$ratio_foldChange = RegionCat_vs_otherRegions_enrichment$ratio_observed/RegionCat_vs_otherRegions_enrichment$ratio_expected

    RegionCat_vs_otherRegions_enrichment = RegionCat_vs_otherRegions_enrichment %>% dplyr::left_join(, y = ChIP_centric_metrics[, c("TF_symbol", "TF_peaks_total", paste0("TF_peaks_in_", c, collapse = ""))])

    RegionCat_vs_otherRegions_enrichment %>% write.table(paste0('../out/tables/', c, '_vs_otherBins_enrichment_ENCODE_tronolabKZFPs.tsv'), col.names = T, row.names = F, sep = '\t')
    }

```

### TODO: what fraction of late RT genes that are RIF1 resistant but H3K9me3 dependent are explained by ZNF274 KO? 

Let's start with the proportion of late RT regions

```R
peaks_x_RTbins %>% head()
```

```R
peaks_x_RTbins %>% head() %>% colnames()
```

```R
peaks_x_RTbins %>% dplyr::select(RTregion_chr, RTregion_start, RTregion_end, RTregion_category, RT_diff_293T_ZNF274KO_vs_WT)  %>% head()
```

```R
RTbins = peaks_x_RTbins %>% dplyr::select(RTregion_chr, RTregion_start, RTregion_end, RTregion_category, RT_diff_293T_ZNF274KO_vs_WT) %>% unique()
RTbins %>% dim()
```
