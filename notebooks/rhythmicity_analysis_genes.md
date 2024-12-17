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

# Estimating the rhythmicity of gene expression throughout the cell cycle, using the 8 gate-sorted bulk RNA-seq in K562
## loading functions

```R
# libraries
library(tidyverse)
library(cosinor)
library(matrixTests)
library(ComplexHeatmap)
library(GEOquery)
library(ggpubr)
library(viridis)
```

```R
# objects

phase = list(g1 = c('G1_1', 'G1_2', 'G1_3'),
g1s = c('G1.S_1', 'G1.S_2', 'G1.S_3'),
s1 = c('S1_1', 'S1_2', 'S1_3'),
s2 = c('S2_1', 'S2_2', 'S2_3'),
s3 = c('S3_1', 'S3_2', 'S3_3'), # s3_2 was filtered out (poor quality)
s3g2 = c('S3.G2_1', 'S3.G2_2', 'S3.G2_3'),
g2 = c('G2_1', 'G2_2', 'G2_3'),
m = c('M_1', 'M_2', 'M_3')
)

phase_to_hours = list(g1 = 5,
g1s = 10,
s1 = 12,
s2 = 14,
s3 = 16, # s3_2 was filtered out (poor quality)
s3g2 = 18,
g2 = 22,
m = 23
)

to_exclude = c('S3_2', 'M_2')

new_phases = c('eG1', 'lG1', 'G1/S', 'S1', 'S2', 'G2', 'G2/M', 'M')

```

```R
# functions

# TPMs
source('../imports/functions/get_tpm_local.R')
source('../imports/functions/get_exonic_lengths.R')
source('../imports/functions/get_pseudocounts.R')

assign_phase <- function(x) {
    for(i in 1:length(phase)) {
        if (x %in% phase[[i]]) {
            return(names(phase)[[i]])
        }
    }
}

assign_gate <- function(x) {
    for(i in 1:length(phase)) {
        if (x %in% phase[[i]]) {
        return(i)
        }
    }
}

assign_stars <- function(pval) {
    if (pval <= 0.05 & pval > 0.01) {
        return('*')
    }
    if(pval <= 0.01 & pval > 0.001) {
        return('**')
    }
    if (pval <= 0.001) {
        return('***')
    }
    return('')
}

# fitting with the batch correction as a covariate
source('../imports//functions/cosinor_covariate.R')
source('../imports//functions/do_regression_cov.R')
source('../imports/functions/utils.asserts.R')
source('../imports/functions/utils.general.R')
source('../imports/functions/utils.tests.R')

# rhythmicity wrapper
estimate_rhythmicity <- function(expr, t, period = 8, cov = NULL) {
    fit = row_cosinor_cov(expr, t, period = period, cov = cov)
    
    fit$padj = p.adjust(p = fit$pvalue, method = 'BH')
    fit$stars = sapply(fit$padj, assign_stars)

    fit$phase_rounded = fit$acrophase %>% round() %% 8
    fit$phase_gating = names(phase)[fit$phase_rounded + 1]
    fit$phase_assigned = factor(new_phases[fit$phase_rounded %>% as.integer() + 1], levels = new_phases)
    
    return(fit)
}


```

### Choosing a normalization and batch effect correction methodology


As a preliminary step, we'll try plotting a few genes from the cell cycle, in log2(tpms)

```R
expr_corrected = read.csv('../data/cell_cycle/norm_vals_adj.xls', row.names = NULL, sep = '\t')
expr_corrected %>% head()
expr_corrected %>% dim()
```

```R
gene_bodies_hg19 = read.table('../data/1812_hg19_ens_all_genes_body.bed', sep = '\t', header = F)
colnames(gene_bodies_hg19) = c("chr", "start", "end", "ensembl", "score", "strand")
gene_bodies_hg19 %>% head()
gene_bodies_hg19 %>% dim()
```

```R
expr_corrected %>% dplyr::filter(!ensembl %in% gene_bodies_hg19$ensembl) %>% head()
```

All the genes are found in the bedfile of gene bodies including non coding genes.

```R
gene_metadata = expr_corrected %>% dplyr::select(ensembl, entrez, symbol, genename)
gene_metadata %>% head()
```

```R
# adding the gene body information
gene_metadata = gene_metadata %>% dplyr::left_join(., gene_bodies_hg19 %>% dplyr::select(chr, start, end, ensembl, strand))
gene_metadata %>% head()
```

Log2(TPM) transformation

```R
expr_corrected$entrez = NULL
expr_corrected$symbol = NULL
expr_corrected$genename = NULL
rownames(expr_corrected) = expr_corrected$ensembl
expr_corrected$ensembl = NULL
expr_corrected = as.matrix(expr_corrected)
```

```R
tpms_corrected = get_tpm_local(expr_corrected)
```

```R
tpms_corrected %>% head()
```

```R
# reloading
expr_corrected = read.csv('../data/cell_cycle/norm_vals_adj.xls', row.names = NULL, sep = '\t')
expr_corrected$entrez = NULL
expr_corrected %>% head()
```

```R
expr_corrected_long = expr_corrected %>% tidyr::pivot_longer(, cols = 4:27)
expr_corrected_long %>% head()
```

```R
# adding the gate information
expr_corrected_long$gate = factor(sapply(X = expr_corrected_long$name, FUN = assign_gate))
```

```R
expr_corrected_long$gate %>% levels()
```

Assigning time points in hours, to see whether it improves the fits

From Romain: 

On a period of about 24 hours, G1 lasts 10-11 hours, S lasts 8 hours, G2 lasts 4 hours and M lasts 1 hour.



```R
expr_corrected_long$phase_hours = phase_to_hours[expr_corrected_long$gate %>% as.integer()] %>% as.double()
```

```R
expr_corrected_long %>% head()
```

```R
expr_corrected_long$phase_corrected = factor(new_phases[expr_corrected_long$gate %>% as.integer()], levels = new_phases)
```

```R
colnames(expr_corrected_long)[4] = 'sample_id'
colnames(expr_corrected_long)[5] = 'norm_adj_counts'
```

```R
expr_corrected_long %>% head()
```

```R
tpms_corrected_long = tpms_corrected %>% as.data.frame() %>% rownames_to_column(var = 'ensembl') %>% tidyr::pivot_longer(cols = 2:25)
colnames(tpms_corrected_long)[2] = 'sample_id'
colnames(tpms_corrected_long)[3] = 'tpms_corrected'
tpms_corrected_long %>% head()

```

```R
expr_corrected_long = dplyr::left_join(expr_corrected_long, tpms_corrected_long)
expr_corrected_long %>% head()
```

```R
expr_corrected_long$pca_outlier = ifelse(expr_corrected_long$sample_id %in% c('S3_2', 'M_2'), yes = T, no = F)
```

```R
cell_cycle_genes = c('CCND3', 'CDK4', 'CCNE2', 'CDK2', 'CCNA2', 'CDK1', 'CCNB1', 'CDC20', 'ORC1', 'MCM6', 'WEE1', 'NUSAP1', 'AURKA', 'CCNB2')
cell_cycle_genes_romain = read.table('../data/cc_romain_genes_16_09_22.csv', skip = 1)$V1
other_genes_of_interest = c('E2F1', 'NFYC', 'E2F4', 'NFYA', 'ZNF202', 'ZNF93', 'ZNF845', 'APOBEC3B')
```

```R
cell_cycle_genes_romain
```

```R
cell_cycle_genes %in% cell_cycle_genes_romain
```

```R
# for testing purposes
cell_cycle_genes = cell_cycle_genes_romain
```

```R
expr_corrected_long %>% head()
```

Exporting for plotting only:

```R
expr_corrected_long %>% colnames()
```

```R
expr_corrected_long %>% dplyr::select(ensembl, symbol, genename, sample_id, norm_adj_counts, gate, phase_corrected, pca_outlier) %>%
write.table('../out/tables/expr_corrected_long.tsv', row.names = F, col.names = T, sep = '\t')
```

```R
options(warn = 1)
dir.create('../out/figures/cell_cycle_marker_genes/')
dir.create('../out/figures/cell_cycle_marker_genes/lognormadjcounts/')
dir.create('../out/figures/cell_cycle_marker_genes/lognormadjcounts_hours/')
dir.create('../out/figures/cell_cycle_marker_genes/logtpms_corrected/')

for(g in c(cell_cycle_genes, other_genes_of_interest)) {
    p = expr_corrected_long %>% dplyr::filter(symbol == g) %>% dplyr::arrange() %>% ggplot(., aes(x = phase_corrected, y = log2(norm_adj_counts + 1), col = pca_outlier)) + 
    geom_point(size = 2) + theme_classic() + 
    labs(x = "Phase", y = "expression [log2]") +
    ggtitle(label = g)
    ggsave(filename = paste0('../out/figures/cell_cycle_marker_genes/lognormadjcounts/', g, '_lognormadjcounts.pdf'), width = 4.5, height = 3)
}

for(g in c(cell_cycle_genes, other_genes_of_interest)) {
    p = expr_corrected_long %>% dplyr::filter(symbol == g) %>% ggplot(., aes(x = phase_corrected, y = log2(norm_adj_counts + 1), col = pca_outlier)) + geom_point(size = 2)
    ggsave(filename = paste0('../out/figures/cell_cycle_marker_genes/lognormadjcounts_hours/', g, '_lognormadjcounts_hours.pdf'), width = 4.5, height = 3)
}

for(g in c(cell_cycle_genes, other_genes_of_interest)) {
    p = expr_corrected_long %>% dplyr::filter(symbol == g) %>% ggplot(., aes(x = phase_corrected, y = log2(tpms_corrected + 1), col = pca_outlier)) + geom_point(size = 2)
    ggsave(filename = paste0('../out/figures/cell_cycle_marker_genes/logtpms_corrected/', g, '_logtpms_corrected.pdf'), width = 4.5, height = 3)
}

```

Sample S2_3 is an outlier, and this can be seen both in TPM-normalized and raw counts. We need to exclude it

Sample M_2 may be a bit of an issue, we may have to remove it. 

Overall, there does not seem to be a benefit in normalizing the data to tpms, norm adj counts seem sufficient for our needs. It would be beneficial to check the raw counts as well: to check whether the total RNA content of the cells becomes an issue in the future.

```R
expr_corrected_long %>% dim()
expr_corrected_long = expr_corrected_long %>% dplyr::filter(!sample_id %in% to_exclude)
expr_corrected_long %>% dim()
```

```R
expr_corrected %>% head()
```

```R
expr_corrected %>% ncol()
expr_corrected = expr_corrected %>% dplyr::select(-one_of(to_exclude))
expr_corrected %>% ncol()
```

```R
tpms_corrected %>% ncol()
tpms_corrected = as.data.frame(tpms_corrected) %>% dplyr::select(-one_of(to_exclude))
tpms_corrected %>% ncol()
```

## Fitting data with the cosinor from matrixStats

Already has the f-test for rhythmicity implemented, no need to implement it on top of cosinor::cosinor.lm. 

The data should be formatted as two matrices, with samples as rows and timepoints as columns. The first matrix contains the gene expression data, and the second one the timepoints in the period.

```R
load('../data/cell_cycle/countTable_gene.RData')
colnames(norm.counts) = gsub('-', '\\.', colnames(norm.counts))
norm.counts %>% head()
```

```R
# limiting to cell cycle genes
cell_cycle_genes %>% sort()
```

```R
cell_cycle_genes %>% sort() %>% length()
```

```R
(gene_metadata %>% dplyr::filter(symbol %in% cell_cycle_genes))$ensembl %>% head()
```

```R
phase_qualitative = factor(as.vector(sapply(colnames(norm.counts), assign_phase)), levels = names(phase))
phase_qualitative
```

```R
phase_hours = phase_to_hours[phase_qualitative] %>% as.double()
phase_hours
```

```R
# todouble for matrixTests
phase_qualitative = phase_qualitative %>% as.double() -1
phase_qualitative
```

```R
phase_qualitative_no_outliers = factor(as.vector(sapply(colnames(norm.counts %>% as.data.frame() %>% dplyr::select(-one_of(to_exclude))), assign_phase)), levels = names(phase))
phase_qualitative_no_outliers
```

```R
phase_qualitative_no_outliers = phase_qualitative_no_outliers %>% as.double() -1
phase_qualitative_no_outliers
```

```R
norm.counts_no_outliers = norm.counts %>% as.data.frame() %>% dplyr::select(-one_of(to_exclude)) %>% as.matrix()
```

```R
fit_arbitrary = estimate_rhythmicity(norm.counts[(gene_metadata %>% dplyr::filter(symbol %in% cell_cycle_genes))$ensembl, ], phase_qualitative, period = 8)
fit_arbitrary %>% arrange(pvalue) %>% head()
fit_arbitrary %>% arrange(pvalue) %>% tail()
```

```R
fit_arbitrary_no_outliers = estimate_rhythmicity(norm.counts_no_outliers[(gene_metadata %>% dplyr::filter(symbol %in% cell_cycle_genes))$ensembl, ], phase_qualitative_no_outliers, period = 8)
fit_arbitrary_no_outliers %>% arrange(pvalue) %>% head()
fit_arbitrary_no_outliers %>% arrange(pvalue) %>% tail()
```

Remarkably, genes known to belong to the cell cycle get adjusted p values close to the detection limit (padj < 0.05).

```R
fit_hours = estimate_rhythmicity(norm.counts[(gene_metadata %>% dplyr::filter(symbol %in% cell_cycle_genes))$ensembl, ], phase_hours, period = 24)
fit_hours %>% arrange(pvalue) %>% head()
fit_hours %>% arrange(pvalue) %>% tail()
```

```R
# fitting with the batch correction as a covariate
batch = Batch[!colnames(norm.counts) %in% to_exclude] %>% as.double() - 1

batch %>% saveRDS('../out/tables/bulk_RNAseq_chronogram_batches.RDS')
```

```R
fit_qualitative_batch = estimate_rhythmicity(norm.counts_no_outliers[(gene_metadata %>% dplyr::filter(symbol %in% cell_cycle_genes))$ensembl, ], 
                                        phase_qualitative_no_outliers, 
                                        period = 8,
                                       cov = batch)
fit_qualitative_batch %>% arrange(pvalue) %>% head()
fit_qualitative_batch %>% arrange(pvalue) %>% tail()
```

Who fits better between the two ways of writing the period (qualitative vs quantitative) ?



```R
dir.create('../out/cell_cycle_figures')
```

```R
stopifnot(all(rownames(fit_arbitrary) == rownames(fit_hours)))
```

```R
p = data.frame(qualitative = fit_arbitrary$rsquared, quantitative = fit_hours$rsquared, ensembl = rownames(fit_arbitrary)) %>% 
    tidyr::pivot_longer(!ensembl, names_to = "fit_type", values_to = "r2") %>% 
    ggplot(aes(x = fit_type, y = r2)) + 
    geom_boxplot() + 
    xlab("") +
    ggpubr::stat_compare_means(paired = T) +
    theme_classic()
p
```

```R
ggsave("../out/cell_cycle_figures/boxplot_r2_qualitative_vs_quantitative_time.svg", p, svg, width = 3.5, height = 4)
```

stat signif to keep it qualitative

```R
boxplot(data.frame(qualitative = fit_arbitrary$rsquared, quantitative = fit_hours$rsquared), ylab = 'R^2')
dev.copy(device = pdf, '../out/cell_cycle_figures/boxplot_r2_qualitative_vs_quantitative_time.pdf', width = 3.5, height = 4)
dev.off()
```

Adding quantitative information (time) does not seem to alter the fit -> we stay with an assumption less and treat samples as on a scale from 1 to 8.

```R
p = data.frame(outliers = fit_arbitrary$rsquared, no_outliers = fit_arbitrary_no_outliers$rsquared, ensembl = rownames(fit_arbitrary)) %>% 
    tidyr::pivot_longer(!ensembl, names_to = "fit_type", values_to = "r2") %>% 
    ggplot(aes(x = fit_type, y = r2)) + 
    geom_boxplot() + 
    xlab("") +
    ggpubr::stat_compare_means(paired = T) +
    theme_classic()
p
```

```R
ggsave("../out/cell_cycle_figures/boxplot_r2_w_wo_outliers.svg", p, svg, width = 2.5, height = 3)
```

Stat signif better to remove outliers -> we drop the outliers

```R
boxplot(data.frame(outliers = fit_arbitrary$rsquared, no_outliers = fit_arbitrary_no_outliers$rsquared), ylab = 'R^2')
dev.copy(device = pdf, '../out/cell_cycle_figures/boxplot_r2_w_wo_outliers.pdf', width = 3.5, height = 4)
dev.off()
```

```R
p = data.frame(no_correction = fit_arbitrary_no_outliers$rsquared, correction = fit_qualitative_batch$rsquared, ensembl = rownames(fit_arbitrary)) %>% 
    tidyr::pivot_longer(!ensembl, names_to = "fit_type", values_to = "r2") %>% 
    ggplot(aes(x = fit_type, y = r2)) + 
    geom_boxplot() + 
    xlab("") +
    ggpubr::stat_compare_means(paired = T) +
    theme_classic()
p
```

```R
ggsave("../out/cell_cycle_figures/boxplot_significance_fit_qualitative_batch.svg", p, svg, width = 2.5, height = 3)
```

Clearly better to correct for batch effect

```R
boxplot(data.frame(no_correction = fit_arbitrary_no_outliers$rsquared, correction = fit_qualitative_batch$rsquared), ylab = 'R^2')
dev.copy(device = pdf, '../out/cell_cycle_figures/boxplot_significance_fit_qualitative_batch.pdf', width = 3.5, height = 4)
dev.off()
```

```R
boxplot(data.frame(no_correction = -log10(fit_arbitrary_no_outliers$pvalue), correction = -log10(fit_qualitative_batch$pvalue)), ylab = '-log10(pvalue)')
dev.copy(device = pdf, '../out/cell_cycle_figures/boxplot_significance_fit_qualitative_batch_significances.pdf', width = 3.5, height = 4)
dev.off()
```

```R
plot(data.frame(no_correction = -log10(fit_arbitrary_no_outliers$pvalue), correction = -log10(fit_qualitative_batch$pvalue)),
     xlab = '-log10(pvalue), without batch covariate',
     ylab = '-log10(pvalue), with batch covariate')
abline(0, 1, col = 'red')
dev.copy(device = pdf, '../out/cell_cycle_figures/scatterplot_significance_fit_qualitative_batch.pdf', width = 3.5, height = 4)
dev.off()
```

Including the batch covariate inreases the additional amount of variance explained by the cosinor.

The correct test is an F-test for overall explained variance:


```R
fit_arbitrary_no_outliers %>% head()
```

```R
fit_qualitative_batch %>% head()
```

```R

```

```R

```

TPMs vs norm adj counts ?

we need to recompute TPMs from the adjusted counts, by going back to "raw count" scale

```R
mean_lib_size = countTable %>% colSums() %>% mean()
mean_lib_size
```

```R
re_scaled_norm.counts = 2^(norm.counts)*mean_lib_size/1e6
re_scaled_norm.counts %>% head()
```

```R
tpms = get_tpm_local(re_scaled_norm.counts)
```

```R
tpms_no_outliers = tpms[, !colnames(tpms) %in% to_exclude]
```

```R
fit_tpms_no_outliers = estimate_rhythmicity(log2(tpms_no_outliers[(gene_metadata %>% dplyr::filter(symbol %in% cell_cycle_genes))$ensembl, ] + 1), phase_qualitative_no_outliers, period = 8)
fit_tpms_no_outliers %>% arrange(pvalue) %>% head()
```

```R
p = data.frame(log2_norm_counts = fit_arbitrary_no_outliers$rsquared, log2_tpms = fit_tpms_no_outliers$rsquared, ensembl = rownames(fit_arbitrary)) %>% 
    tidyr::pivot_longer(!ensembl, names_to = "fit_type", values_to = "r2") %>% 
    ggplot(aes(x = fit_type, y = r2)) + 
    geom_boxplot() + 
    xlab("") +
    ggpubr::stat_compare_means(paired = T) +
    theme_classic()
p
```

```R
ggsave("../out/cell_cycle_figures/boxplot_no_outliers_tpms_vs_counts.svg", p, svg, width = 2.5, height = 3)
```

The fit is vastly better with the norm counts than the TPMs. Normalizing for total library size is probably too brutal.

```R
boxplot(data.frame(log2_norm_counts = fit_arbitrary_no_outliers$rsquared, log2_tpms = fit_tpms_no_outliers$rsquared), ylab = 'R^2')
dev.copy(device = pdf, '../out/cell_cycle_figures/boxplot_no_outliers_tpms_vs_counts.pdf', width = 4.5, height = 4)
dev.off()
```

What about non-normalized, raw counts? Does the sum of normalized counts increase along the cell cycle?

```R
phase %>% unlist() %>% as.character()
```

```R
(norm.counts %>% colSums())[phase %>% unlist() %>% as.character()] %>% plot()
```

```R
raw.counts = countTable
colnames(raw.counts) = gsub("-", "\\.", colnames(raw.counts))
raw.counts %>% head()
```

```R
raw.counts_no_outliers = raw.counts[, !colnames(raw.counts) %in% to_exclude]

```

```R
# fitting on log(raw counts + 1)
fit_qualitative_batch_rawCounts = estimate_rhythmicity(log2(raw.counts_no_outliers[(gene_metadata %>% dplyr::filter(symbol %in% cell_cycle_genes))$ensembl, ] + 1), 
                                        phase_qualitative_no_outliers, 
                                        period = 8,
                                       cov = batch)
```

```R
p = data.frame(log2_raw_counts = fit_qualitative_batch_rawCounts$rsquared, log2_norm_counts = fit_arbitrary_no_outliers$rsquared, log2_tpms = fit_tpms_no_outliers$rsquared, ensembl = rownames(fit_arbitrary)) %>% 
    tidyr::pivot_longer(!ensembl, names_to = "fit_type", values_to = "r2") %>% 
    ggplot(aes(x = fit_type, y = r2)) + 
    geom_boxplot() + 
    xlab("") +
    ggpubr::stat_compare_means(paired = T, comparisons = list(c(1, 2), c(1, 3))) +
    theme_classic()
p
```

```R
ggsave("../out/cell_cycle_figures/boxplot_no_outliers_tpms_vs_counts.svg", p, svg, width = 4.5, height = 4)
```

```R
boxplot(data.frame(log2_raw_counts = fit_qualitative_batch_rawCounts$rsquared, log2_norm_counts = fit_arbitrary_no_outliers$rsquared, log2_tpms = fit_tpms_no_outliers$rsquared), ylab = 'R^2')
dev.copy(device = pdf, '../out/cell_cycle_figures/boxplot_no_outliers_tpms_vs_counts.pdf', width = 6, height = 6)
dev.off()
```

## Fitting on all genes

The genes in norm.counts are not only protein coding genes. We could trade off on non-coding genes and increase our statistical power by removing them before estimating anything.

```R
gene_metadata %>% dim()
```

```R
coding_genes = read.table('../data/1807_hg19_ens_coding_genes_body_symbol.bed', header = F, sep = '\t')
coding_genes %>% dim()
coding_genes %>% head()
```

```R
gene_metadata$coding = ifelse((gene_metadata$ensembl %in% coding_genes$V4) & (gene_metadata$symbol != "---"), yes = T, no = F)
```

```R
table(gene_metadata$coding)
```

```R
gene_metadata %>% write.table('../out/tables/gene_metadata_for_rhythmicity.tsv', sep = '\t', quote = F, row.names = F, col.names = T)
```

```R
norm.counts_no_outliers %>% head()
```

```R
norm.counts_no_outliers_coding = norm.counts_no_outliers[(rownames(norm.counts_no_outliers) %in% (gene_metadata %>% dplyr::filter(coding) %>% pull(ensembl))), ]
```

```R
norm.counts_no_outliers_coding %>% dim()
```

```R
norm.counts_no_outliers_coding %>% head()
```

```R
norm.counts_no_outliers_coding %>% write.table('../out/tables/norm.counts_no_outliers_coding.tsv', sep = '\t', row.names = T, col.names = T, quote = F)
```

```R
phase_qualitative_no_outliers %>% saveRDS('../out/tables/phase_qualitative_no_outliers.RDS')
```

```R
estimate_rhythmicity
```

```R
fit_qualitative = estimate_rhythmicity(norm.counts_no_outliers_coding, phase_qualitative_no_outliers, period = 8, batch)

# adding the gene metadata
fit_qualitative = gene_metadata %>% dplyr::filter(coding) %>% left_join(fit_qualitative %>% tibble::rownames_to_column('ensembl'))


fit_qualitative %>% arrange(padj) %>% head()
```

```R
fit_qualitative %>% dim()
```

```R
(fit_qualitative$padj < 0.05) %>% sum()
```

```R
fit_qualitative %>% dplyr::filter(symbol == "ZNF695")
```

```R
# adding the replication timing
```

```R
# adding the replication time point
replitime = read.table('../data/cell_cycle/BaseFile_ClusteredHeatmap.csv', sep = ';', header = F)
```

```R
colnames(replitime) = str_split("chrom;start;end;name;score;strand;thickStart;thickEnd;itemRGB;blockCount;blockSizes;blockStart;RepliTiming;deepTools_group", ";", simplify = T) %>% as.vector()
```

```R
replitime %>% head()
```

```R
replitime %>% dplyr::select(RepliTiming, deepTools_group) %>% unique()
```

```R
(fit_qualitative$ensembl %in% replitime$name) %>% sum()
```

```R
fit_qualitative = dplyr::left_join(fit_qualitative, replitime %>% dplyr::select(name, RepliTiming), by = c("ensembl"="name"))
```

```R
fit_qualitative$RepliTiming = ifelse(test = is.na(fit_qualitative$RepliTiming), yes = "NA", no = fit_qualitative$RepliTiming)
```

```R
fit_qualitative$RepliTiming %>% unique() %>% sort()
```

There are no NAs. That's a really good news.

```R
# renaming the replication times to G1, S1-4, G2, as in the original PNAS publication: https://pubmed.ncbi.nlm.nih.gov/19966280/
repliTimingTranslation = c("G1", "S1", "S2", "S3", "S4", "G2")
names(repliTimingTranslation) = fit_qualitative$RepliTiming %>% unique() %>% sort()
fit_qualitative$RepliTiming = sapply(fit_qualitative$RepliTiming, function(x) return(repliTimingTranslation[x])) %>% as.vector()
```

```R
# replitiming to factor, to have them in the correct order
fit_qualitative$RepliTiming = factor(fit_qualitative$RepliTiming, levels = repliTimingTranslation)
```

```R
# continuous replication timing
replitiming_k562_cont = read.table('../out/RTregions_vs_gene_bodies/RTindex_K562.bed', sep = '\t', header = F)
replitiming_k562_cont %>% head()
```

```R
RTindexK562_gene_bodies = replitiming_k562_cont %>% 
    dplyr::select(V4, V11, V12) %>% 
    dplyr::rename(ensembl = V4, RTindex = V11, overlap = V12) %>% 
    unique() %>% dplyr::group_by(ensembl) %>% 
    dplyr::mutate(RTindex_avg = (as.double((RTindex %*% overlap)))/sum(overlap)) %>%
    dplyr::select(ensembl, RTindex_avg) %>% unique()
```

Sanity check:

```R
replitiming_k562_cont %>% head(10)
```

```R
RTindexK562_gene_bodies %>% dplyr::filter(ensembl == "ENSG00000187961")
```

```R
(3.684064*4033 + 4.059874*1094)/(4033+1094)
```

```R
fit_qualitative = fit_qualitative %>% dplyr::left_join(RTindexK562_gene_bodies)
```

```R
# exporting the gene bodies sorted by acrophase and start coordinate for all genes.
fit_qualitative %>% 
    dplyr::arrange(acrophase, chr, start) %>% 
    select(chr, start, end, ensembl, phase_assigned, strand) %>% 
    write.table('../out/gene_bodies_all_genes_sorted_acrophase_tssstart.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# adding in the RIF1 enrichment over the gene body
rif1_gene_bodies = read.table('../out/temp/gene_bodies_all_genes_sorted_acrophase_tssstart_hg38_rif1.tab', sep = '\t', skip=3)
rif1_gene_bodies %>% head()
```

```R
fit_qualitative %>% dim()
```

```R
rif1_gene_bodies %>% dim()
```

```R
# solution is to load the hg38 bed, to get the exact sequence names
gene_bodies_hg38 = read.table('../out/gene_bodies_all_genes_sorted_acrophase_tssstart_hg38.bed', sep = '\t', header = F)
hg19_genes_liftovered_on_hg38 = gene_bodies_hg38$V4
hg19_genes_liftovered_on_hg38 %>% length()
```

```R
gene_bodies_hg38 %>% duplicated() %>% length()
```

```R
gene_bodies_hg38 %>% duplicated() %>% sum()
```

```R
rif1_score = rif1_gene_bodies %>% apply(., 1, mean, na.rm = T)
rif1_score = as.data.frame(rif1_score)
rif1_score$ensembl = hg19_genes_liftovered_on_hg38
```

```R
rif1_score = rif1_score %>% dplyr::rename(rif1 = rif1_score)
```

```R
fit_qualitative = fit_qualitative %>% dplyr::left_join(., rif1_score)
```

```R
# which KZFPs are cycling ?
kzfps = read.csv('../data/human_KZFPTable.csv', sep = '\t')
kzfps %>% head(1)
```

```R
kzfps %>% dim()
```

```R
# filtering on KZFPs
kzfps = kzfps %>% dplyr::filter(KRABid != '')
kzfps %>% dim()
```

```R
colnames(kzfps)
```

```R
kzfps$chromosome %>% unique()
```

This is the hg19 assembly, but we need to convert chromosome names, as indicated on this page: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&chromInfoPage=

```R
chr_conversion = read.table("../data/hg19.chromAlias.txt", sep = '\t', header = T)
chr_conversion %>% head()
```

```R
kzfps = kzfps %>% left_join(., chr_conversion %>% dplyr::select(ucsc, refseq), by = c("chromosome" = "refseq")) %>% dplyr::rename(chr = ucsc)
kzfps %>% head()
```

```R
kzfps$chr %>% unique()
```

```R
# adding ZNF812, which is NC_000019.9_KZFP_261
kzfps$GeneName[which(kzfps$X == "NC_000019.9_KZFP_261")] = "ZNF812"
```

```R
# extracting ensemblIDs
kzfp_symbols = c(kzfps$GeneName %>% unique(), "ZNF496", "ZNF79")
kzfp_symbols %>% length()
```

None of these protein-coding KZFPs (disagreements between hg38 and h19 annotations) are in the KRABOPEDIA list. We ignore them for now.

```R
fit_kzfps = fit_qualitative %>% dplyr::filter(symbol %in% kzfp_symbols)
fit_kzfps %>% dim()
```

```R
# importing KZFPs from Jonas's list: doi: 10.1101/gr.277722.123
kzfps_jonas = read.table('../data/kzfps_jonas.csv', header = T, sep = ';', quote = "")
kzfps_jonas %>% dim()

```

```R
# detected KZFPs in the RNA-seq, that are absent from the KZFP screening results.
kzfps_jonas %>% dplyr::filter(assigned_gene %in% fit_qualitative$symbol, !assigned_gene %in% fit_kzfps$symbol)
```

```R
fit_kzfps %>% dplyr::filter(symbol %in% c("ZNF496", "ZNF79"))
```

```R
# looking for duplicated KZFPs
kzfp_symbols = kzfps %>% dplyr::filter(GeneName %in% fit_qualitative$symbol, 
                        chr %in% c(paste0("chr", 1:22), "chrX", "chrY")) %>% dplyr::select(GeneName) %>% dplyr::arrange(GeneName) %>% unlist() %>% as.character()
kzfp_symbols
```

```R
kzfp_symbols[which(kzfp_symbols %>% duplicated)]
```

These duplicated KZNFs are only so due to several KRAB domains, while the 3' ZFA remain the same. We can thus extract them using a "first" window function

```R
kzfp_zfcoordinates = kzfps %>% dplyr::filter(GeneName %in% fit_qualitative$symbol, 
                        chr %in% c(paste0("chr", 1:22), "chrX", "chrY")) %>% 
    dplyr::group_by(GeneName) %>% dplyr::arrange(GeneName) %>% 
    dplyr::filter(row_number() == 1) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, KRABStart, KRABEnd, GeneName, CanonicalsZF, Strand) %>%
    dplyr::arrange(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd)
```

```R
kzfp_zfcoordinates %>% head()
kzfp_zfcoordinates %>% dim()
```

Adding ZNF718 by hand which was not found by the KZFP scanning script, probably because its KRAB domain overlaps that of ZNF595, meaning that the algorithm did not go looking further than the first array of ZFPs (the 595 one). We add the KRAB coordinates for a control of a transcribed domain that is not the ZNF, to see if cycling explains some unspecific binding. 

The KRAB domain coordinates of ZNF718 is the same as the one of ZNF595

```R
kzfps %>% dplyr::filter(GeneName=="ZNF595") %>% dplyr::select(KRABStart, KRABEnd)
```

```R
znf718_zfcoordinates = list("chr4", 154983, 155894, 49435, 49557, "ZNF718", 10, "+")
# not added anymore, ZNF718 is not a protein-coding KZFP gene
```

Adding the KRAB and ZFA coordinates for ZNF496 and ZNF79

```R
znf496_zfcoordinates = list("chr1", 247463860, 247464369, 247473000, 247473743, "ZNF496", 4, "-")
znf79_zfcoordinates = list("chr9", 130206556, 130207464, 130197369, 130198281, "ZNF79", 11, "+")

kzfp_zfcoordinates = rbind(kzfp_zfcoordinates, znf496_zfcoordinates, znf79_zfcoordinates)
kzfp_zfcoordinates %>% head()
```

```R
kzfp_zfcoordinates = kzfp_zfcoordinates %>% arrange(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd)
```

```R
fit_kzfps = fit_kzfps %>% dplyr::left_join(., kzfp_zfcoordinates, by = c("symbol" = "GeneName"))
```

KZFPs with missing TSS as per the join between the TSS table and the KZFP table:
- ZNF761: ENSG00000160336: Processed transcript
- ZNF718: ENSG00000250312: LincRNA
- ZNF316: ENSG00000205903: Pseudogene
- ZNF888: ENSG00000213793: LincRNA
- ZNF705E: ENSG00000214534: pseudogene
- ZNF761: ENSG00000160336: processed transcript -> not in our TSS list, which is only protein coding genes


We thus exclude them from the list of KZFP genes

```R
# excluding pseudogenized KZFPs, which do not code for proteins:

fit_kzfps = fit_kzfps %>% dplyr::filter(!symbol %in% c("ZNF761", "ZNF718", "ZNF316", "ZNF888", "ZNF705E", "ZNF761"))
```

```R
# KZFPs that are present in the scanning list, but not in Jonas's list
fit_kzfps %>% dplyr::filter(!symbol %in% kzfps_jonas$assigned_gene) %>% print()
```

Could the KZFP overlapping RAB31 be this one? ENSG00000266127 aka RP11-474N24.2-001

RMDN2 doesnt seem to overlap any KZFP, or any transcript ressembling that. Must be a dead one, not even annotated as a transcript.

We will have to add cluster information for ZNF496 and ZNF79 ourselves.

```R
fit_kzfps %>% dim()
fit_kzfps = fit_kzfps %>% dplyr::filter(!symbol %in% c("RAB31", "RMDN2"))
fit_kzfps %>% dim()

```

```R
# adding cluster information
kzfps_jonas %>% colnames()
fit_kzfps %>% colnames()
```

```R
stopifnot(all(fit_kzfps$chr.x == fit_kzfps$chr.y))
```

```R
fit_kzfps$chr.y = NULL
fit_kzfps = fit_kzfps %>% dplyr::rename(chr = chr.x)
```

```R
fit_kzfps = fit_kzfps %>% dplyr::left_join(., kzfps_jonas %>% dplyr::select(chr, assigned_gene, strand, age_MA, species, z_C2H2_miss, cluster), by = c("symbol"="assigned_gene", "Strand"="strand", "chr"="chr"))
```

```R
fit_kzfps %>% head()
```

ZNF79 is not in any cluster, it is isolated on chromosome 9
ZNF496 is in a cluster together with ZNF669, ZNF124 and others.

```R
znf79_idx = which(fit_kzfps$symbol == "ZNF79")
znf496_idx = which(fit_kzfps$symbol == "ZNF496")
cluster_idx = which(colnames(fit_kzfps) == "cluster")
```

```R
fit_kzfps[znf79_idx, cluster_idx] = "noCluster"
fit_kzfps[znf496_idx, cluster_idx] = "chr1.2"

```

```R
fit_kzfps$cluster = ifelse(is.na(fit_kzfps$cluster), yes = "noCluster", no = fit_kzfps$cluster)
```

```R
cluster_names = fit_kzfps$cluster %>% unique() %>% sort()
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

```R
stopifnot(all((kzfp_cluster_order %>% sort()) == cluster_names))
```

```R
fit_kzfps$cluster = factor(fit_kzfps$cluster, levels = kzfp_cluster_order)
```

```R
fit_kzfps$cluster %>% table()
```

```R
# exporting KZFP zinc finger domains in the same order as the rhythmicity data
fit_kzfps %>% 
    arrange(acrophase) %>% 
    dplyr::select(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, symbol, CanonicalsZF, Strand) %>% 
    write.table('../out/kzfp_zfcoordinates_sorted_acrophase.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# exporting KZFP zinc finger domains, sorted by MESOR:
fit_kzfps %>% 
    arrange(mesor) %>% 
    dplyr::select(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, symbol, CanonicalsZF, Strand) %>% 
    write.table('../out/kzfp_zfcoordinates_sorted_MESOR.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# exporting KZFP KRAB domains for control of a transcribed region, that should not bind ZNF274, to see if cell cycle explains the binding at transcribed loci.
fit_kzfps %>%
    arrange(acrophase) %>% 
    dplyr::select(chr, KRABStart, KRABEnd, symbol, CanonicalsZF, Strand) %>% 
    write.table('../out/kzfp_krabcoordinates_sorted_acrophase.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# exporting KZFPs zinc finger domains for signif. cycling KZFPs only, as in the rhythmicity data:
fit_kzfps %>%
    dplyr::filter(padj < 0.05) %>%
    arrange(acrophase) %>% 
    dplyr::select(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, symbol, CanonicalsZF, Strand) %>% 
    write.table('../out/kzfp_signif_cycling_zfcoordinates_sorted_acrophase.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# exporting KZFP KRAB domains for control of a transcribed region, that should not bind ZNF274, to see if cell cycle explains the binding at transcribed loci.
fit_kzfps %>%
    dplyr::filter(padj < 0.05) %>%
    arrange(acrophase) %>% 
    dplyr::select(chr, KRABStart, KRABEnd, symbol, CanonicalsZF, Strand) %>% 
    write.table('../out/kzfp_signif_cycling_krabcoordinates_sorted_acrophase.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# exporting KZFP ZNF domain sorted by chromosome and start position
fit_kzfps %>% 
    arrange(cluster, chr, OutOfFrameZFAStart)%>% 
    dplyr::select(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, symbol, CanonicalsZF, Strand) %>% 
    write.table('../out/kzfp_zfcoordinates_sorted_coordinates.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# exporting KZFP ZNF domain sorted by cluster and acrophase
fit_kzfps %>% 
    arrange(cluster, acrophase)%>% 
    dplyr::select(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, symbol, CanonicalsZF, Strand) %>% 
    write.table('../out/kzfp_zfcoordinates_sorted_cluster_acrophase.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# exporting KZFP gene bodies, sorted by cluster and acrophase, to be used with RIF1 binding.
fit_kzfps %>% 
    arrange(cluster, acrophase)%>% 
    dplyr::select(chr, start, end, symbol, CanonicalsZF, strand) %>% 
    write.table('../out/kzfp_genebody_sorted_cluster_acrophase.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# KZFP TSS to plot the histone mark in their vicinity:
tss = read.table('../data/tss_clustered.bed', sep = '\t', header = F)
colnames(tss) = c("chr", "TSSstart", "TSSend", "ensembl", "score", "strand")
tss %>% head()
```

Exporting the TSS for all rhythmic genes, to check whether acrophase matches with a specific ATAC-seq / H3K4me3 / H3K4me1 signal.

```R
tss$ensembl %>% unique() %>% length()
```

```R
tss_all_genes = fit_qualitative %>% arrange(acrophase) %>% dplyr::left_join(., tss %>% dplyr::select(chr, TSSstart, TSSend, ensembl, strand)) %>% dplyr::filter(!is.na(chr))
tss_all_genes %>% head()
tss_all_genes %>% dim()
```

```R
tss_all_genes = tss_all_genes %>% dplyr::filter(!is.na(TSSstart))
```

```R
tss_all_genes %>% dim()
```

```R
tss_all_genes$ensembl %>% unique() %>% length()
```

```R
tss_rhythmic_genes = tss_all_genes %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase)
```

```R
tss_rhythmic_genes$ensembl %>% unique() %>% length()
```

```R
tss_rhythmic_genes %>% dim()
```

```R
#Exporting the TSS, sorted by acrophase and start coordinate, for all genes
tss_all_genes %>% 
    dplyr::arrange(acrophase, TSSstart) %>%
    dplyr::select(chr, TSSstart, TSSend, symbol, phase_assigned, strand) %>% 
    write.table('../out/tss_all_genes_sorted_acrophase_tssstart.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
tss_all_genes %>% dplyr::filter(is.na(TSSstart)) %>% dim()
```

Exporting the TSS for all rhythmic genes, to check whether acrophase matches with a specific ATAC-seq / H3K4me3 / H3K4me1 signal.

```R
h3k4me3_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k4me3.tab', sep = '\t', skip=3)
tss_all_genes$h3k4me3 = h3k4me3_all_genes_tss %>% rowSums(na.rm = T)

atac_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_atac.tab', sep = '\t', skip=3)
tss_all_genes$atac = atac_all_genes_tss %>% rowSums(na.rm = T)

h3k27ac_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_acrophase_h3k27ac.tab', sep = '\t', skip=3)
tss_all_genes$h3k27ac = h3k27ac_all_genes_tss %>% rowSums(na.rm = T)

h3k4me1_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_acrophase_h3k4me1.tab', sep = '\t', skip=3)
tss_all_genes$h3k4me1 = h3k4me1_all_genes_tss %>% rowSums(na.rm = T)

h3k9me3_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_acrophase_h3k9me3.tab', sep = '\t', skip=3)
tss_all_genes$h3k9me3 = h3k9me3_all_genes_tss %>% rowSums(na.rm = T)

h3k9me3_domains_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_acrophase_h3k9me3_domains.tab', sep = '\t', skip=3)
tss_all_genes$h3k9me3_domains = h3k9me3_domains_all_genes_tss %>% rowSums(na.rm = T)

h3k27me3_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_h3k27me3.tab', sep = '\t', skip=3)
tss_all_genes$h3k27me3 = h3k27me3_all_genes_tss %>% rowSums(na.rm = T)

ctcf_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_ctcf.tab', sep = '\t', skip=3)
tss_all_genes$ctcf = ctcf_all_genes_tss %>% rowSums(na.rm = T)

znf274_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_znf274_m01.tab', sep = '\t', skip=3)
tss_all_genes$znf274 = znf274_all_genes_tss %>% rowSums(na.rm = T)

znf75a_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_znf75a.tab', sep = '\t', skip=3)
tss_all_genes$znf75a = znf75a_all_genes_tss %>% rowSums(na.rm = T)

znf75d_all_genes_tss = read.table('../out/temp/tss_all_genes_sorted_acrophase_tssstart_znf75d.tab', sep = '\t', skip=3)
tss_all_genes$znf75d = znf75d_all_genes_tss %>% rowSums(na.rm = T)
```

```R
tss_all_genes %>% head()
tss_all_genes %>% dim()
```

```R
tss_all_genes$ensembl %>% unique() %>% length()
```

```R
# adding the ENCODE counts for RNA-seq
k562_encode_counts = read.table('../data/ENCFF352KRP.tsv', header = T, sep = '\t')
k562_encode_counts %>% dim()
k562_encode_counts %>% head()
```

```R
k562_encode_counts$gene_id %in% tss_all_genes$entrez %>% sum() 
```

TODO: why don't the entrez match?

```R
k562_encode_counts$gene_id %>% unique() %>% head()
```

```R
tss_all_genes$entrez %>% head()
```

```R
tss_all_genes$entrez %in% k562_encode_counts$gene_id %>% sum()
```

There is basically no overlap between the two sets, I'm confused as to why...

```R
tss_all_genes %>% ggplot(aes(x = log(h3k4me3), y = log(h3k27ac))) + geom_point()
```

```R
tss_all_genes$mesor %>% head()
```

```R
(tss_all_genes %>% dplyr::mutate(mesor_linear = exp(mesor)) %>% dplyr::select(h3k4me3, h3k27ac, atac, h3k4me1, h3k9me3, h3k9me3_domains, h3k27me3, ctcf, mesor_linear) + 1) %>% log() %>% cor(method = "pearson") %>% heatmap()
```

At the level of all TSS (before summerazing to the most active TSS), H3k4me3, ATAC, H3k27ac and H3k4me1 form a cluster of high correlation on their own. MESOR clusters independently, and H3k27me3 and H3k9me3 cluster together, though less strongly than activatory marks.

```R
# how many TSS must we remove?
tss_all_genes %>% dim()
```

```R
tss_all_genes$ensembl %>% unique() %>% length()
```

```R
tss_top_h3k4me3_h3k27ac_all_genes = tss_all_genes %>% dplyr::group_by(ensembl) %>% dplyr::slice_max(order_by = h3k4me3, na_rm = T, n = 1) %>% ungroup()
tss_top_h3k4me3_h3k27ac_all_genes %>% dim()
```

```R
tss_top_h3k4me3_h3k27ac_all_genes %>% colnames()
```

```R
any(is.na(tss_top_h3k4me3_h3k27ac_all_genes$rif1))
```

```R
sum(is.na(tss_top_h3k4me3_h3k27ac_all_genes$rif1))
```

```R
sum(is.na(tss_top_h3k4me3_h3k27ac_all_genes$RTindex_avg))
```

Since RIF1 comes in enrichment values, it makes sense to replace these by a zero anyways.

```R
tss_top_h3k4me3_h3k27ac_all_genes$rif1 = ifelse(is.na(tss_top_h3k4me3_h3k27ac_all_genes$rif1), yes = 0, no = tss_top_h3k4me3_h3k27ac_all_genes$rif1)
tss_top_h3k4me3_h3k27ac_all_genes$RTindex_avg = ifelse(is.na(tss_top_h3k4me3_h3k27ac_all_genes$RTindex_avg), yes = 0, no = tss_top_h3k4me3_h3k27ac_all_genes$RTindex_avg)
```

```R
# deciding ties with h3k27ac levels
tss_top_h3k4me3_h3k27ac_all_genes = tss_top_h3k4me3_h3k27ac_all_genes %>% dplyr::group_by(ensembl) %>% dplyr::slice_max(order_by = h3k27ac, na_rm = T, n = 1) %>% ungroup()
tss_top_h3k4me3_h3k27ac_all_genes %>% dim()
```

```R
tss_top_h3k4me3_h3k27ac_all_genes$entrez %>% unique() %>% length()
```

What are the 400 or so genes which don't have corresponding entrez and ensembl ID?

```R
tss_top_h3k4me3_h3k27ac_all_genes %>% dplyr::filter(duplicated(entrez))
```

```R
tss_top_h3k4me3_h3k27ac_all_genes %>% ggplot(aes(x = log(h3k4me3), y = log(h3k9me3))) + geom_point()
```

```R
tss_top_h3k4me3_h3k27ac_all_genes %>% ggplot(aes(x = log(h3k4me3), y = log(h3k9me3_domains))) + geom_point()
```

```R
tss_top_h3k4me3_h3k27ac_all_genes %>% ggplot(aes(x = log(h3k9me3), y = log(h3k9me3_domains))) + geom_point()
```

```R
tss_top_h3k4me3_h3k27ac_all_genes %>% dplyr::select(h3k9me3, h3k9me3_domains) %>% cor(., method = "spearman")
```

Only 0.47 correlation between H3K9me3 at TSS and H3K9me3 domains surrounding 25kb of the tss...

```R
library(circlize)
htmp = (tss_top_h3k4me3_h3k27ac_all_genes %>% 
        dplyr::filter(padj < 0.05) %>%
     dplyr::mutate(mesor_linear = exp(mesor), 
                   amplitude_linear = exp(amplitude), 
                   signif = -1*padj, 
                   acrophase_linear = exp(acrophase),
                  rtindex_linear = exp(RTindex_avg)) %>% 
     dplyr::select(h3k4me3, h3k27ac, atac, h3k4me1, h3k9me3, h3k27me3, mesor_linear, amplitude_linear, signif, acrophase_linear, rtindex_linear) + 1) %>% 
    dplyr::rename(H3K4me3 = h3k4me3, H3K27ac = h3k27ac, ATAC = atac, H3K4me1 = h3k4me1, H3K9me3 = h3k9me3, H3K27me3 = h3k27me3, MESOR = mesor_linear, amplitude = amplitude_linear, `signif.`=signif, acrophase = acrophase_linear, RTindex = rtindex_linear) %>%
    log() %>% 
    cor(method = "spearman") %>% 
    ComplexHeatmap::Heatmap(name = "SCC", col = colorRamp2(c(-1, 0, 1), rev(c('#f1a340','#f7f7f7','#998ec3'))),
                           row_split = 7, column_split = 7)

svg('../out/cell_cycle_figures/promoter_state_vs_cosinor_params_heatmap_rhythmic_genes.svg', width = 4, height = 3.5)
draw(htmp)
dev.off()
```

```R
stopifnot((tss_top_h3k4me3_h3k27ac_all_genes$ensembl %>% unique %>% length()) == nrow(tss_top_h3k4me3_h3k27ac_all_genes))
```

```R
htmp = (tss_top_h3k4me3_h3k27ac_all_genes %>% 
     dplyr::mutate(mesor_linear = exp(mesor), 
                   amplitude_linear = exp(amplitude), 
                   signif = -1*padj, 
                   acrophase_linear = exp(acrophase),
                  rtindex_linear = exp(RTindex_avg)) %>% 
     dplyr::select(h3k4me3, h3k27ac, atac, h3k4me1, h3k9me3, h3k27me3, mesor_linear, amplitude_linear, signif, acrophase_linear, rtindex_linear) + 1) %>% 
    dplyr::rename(H3K4me3 = h3k4me3, H3K27ac = h3k27ac, ATAC = atac, H3K4me1 = h3k4me1, H3K9me3 = h3k9me3, H3K27me3 = h3k27me3, MESOR = mesor_linear, amplitude = amplitude_linear, `signif.`=signif, acrophase = acrophase_linear, RTindex = rtindex_linear) %>%
    log() %>% 
    cor(method = "spearman") %>% 
    ComplexHeatmap::Heatmap(name = "SCC", col = colorRamp2(c(-1, 0, 1), rev(c('#f1a340','#f7f7f7','#998ec3'))),
                           row_split = 8, column_split = 8)

svg('../out/cell_cycle_figures/promoter_state_vs_cosinor_params_heatmap_all_genes.svg', width = 4, height = 3.5)
draw(htmp)
dev.off()
```

H3K4me3, ATAC and H3K27Ac cluster together with MESOR. The link between acrophase and repli timing is lost when not looking at significant genes.


### Exploring KZFP expression with respect to ChIP-seq marks

```R
kzfp_tss = fit_kzfps %>% 
    dplyr::left_join(., tss_top_h3k4me3_h3k27ac_all_genes %>% 
    dplyr::select(ensembl, TSSstart, TSSend, h3k4me3, h3k27ac, atac, h3k4me1, h3k27me3, h3k9me3, ctcf, znf274, znf75a, znf75d)) %>% 
    dplyr::rename(h3k9me3_tss = h3k9me3, znf274_tss = znf274, znf75a_tss = znf75a, znf75d_tss = znf75d) %>% 
    arrange(acrophase)
```

```R
kzfp_tss %>% dim()
```

None are found in the tss table because they are not annotated as protein-coding genes. 

None of these are on KRABOPEDIA, they are not considered bona fide KZFP genes and we have no ChIP-seq data for them.

```R
krabopedia_kzfp_table = read.table('../data/krabopedia_list.tsv', sep = '\t', header = T)
krabopedia_kzfps = (krabopedia_kzfp_table %>% as.vector() %>% unlist() %>% sapply(., function(x) gsub(" ", "", x)) %>% as.character() %>% sort() %>% unique())[-1]
```

```R
krabopedia_kzfps
```

```R
krabopedia_kzfps[which(!krabopedia_kzfps %in% fit_kzfps$symbol)]
krabopedia_kzfps[which(!krabopedia_kzfps %in% fit_kzfps$symbol)] %>% length()
```

It is expected that all ZSCAN proteins, which do not contain KRAB domains, are not in our list. Furthermore, we should look at the divergence between the krabopedia list and the list obtained from the last run of my script.

```R
krabopedia_kzfps[which(!krabopedia_kzfps %in% kzfps$GeneName)] 
krabopedia_kzfps[which(!krabopedia_kzfps %in% kzfps$GeneName)] %>% length()
```

Again, all ZSCANs are expected not to be there as they do not have KRAB domains so we can remove them.

```R
missing_kzfps = krabopedia_kzfps[which(!krabopedia_kzfps %in% kzfps$GeneName)]
missing_kzfps = missing_kzfps[which(!grepl(pattern = "ZSCAN", missing_kzfps))]
```

```R
missing_kzfps
missing_kzfps %>% length()
```

- AC010642.1 is a readthrough transcript between ZNF8 and an ERVK and only have 1 zinc finger...
- AC073343.1 only has one zinc finger
- CTD-2192J16.11 is a LncRNA
- HKR1 is ZNF875, and has a variant KRAB domain. It is in the list under ZNF875 name.
- RP11-1396O13.13 has only one zinc finger.
- ZNF121 has no KRAB domain, even though it is annotated as such in Genecards, when looking at gene family
- ZNF16 has no KRAB domain
- ZNF200 has no KRAB domain
- ZNF260 has no KRAB domain
- ZNF286B has no KRAB domain
- ZNF35 has no KRAB domain
- ZNF497 has no KRAB domain
- ZNF625-ZNF20: we have both separately in our list
- ZNF660 has no KRAB domain
- ZNF724P: initially annotated as a pseudogene. Now bona fide KZFP, we have it under the name ZNF724.
- ZNF768 as no KRAB domain
- ZNF788: ZNF788P in some annotations, putative KRAB domain protein, considered a pseudogene.
- ZNF80 has no KRAB domain
- ZNF837 has no KRAB domain.

```R
stopifnot(all((kzfp_tss %>% arrange(acrophase))$ensembl == (kzfp_tss %>% arrange(acrophase, TSSstart))$ensembl))
```

```R
# exporting TSS in the same order as the acrophase:
kzfp_tss %>% 
    arrange(acrophase) %>% 
    dplyr::select(chr, TSSstart, TSSend, symbol, phase_assigned, Strand) %>% 
    write.table('../out/kzfp_tss_sorted_acrophase.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
# exporting TSS in increasing MESOR order:
kzfp_tss %>% 
    arrange(mesor) %>% 
    dplyr::select(chr, TSSstart, TSSend, symbol, phase_assigned, Strand) %>% 
    write.table('../out/kzfp_tss_sorted_MESOR.bed', row.names=F, col.names = F, sep = '\t', quote=F)
```

```R
kzfp_tss = kzfp_tss %>% arrange(acrophase)
```

```R
h3k9me3_all_kzfp_znfs = read.table('../out/temp/kzfp_znfs_strict_sorted_acrophase_h3k9me3.tab', sep = '\t', skip=3)
h3k9me3_all_kzfp_znfs %>% head()
kzfp_tss = kzfp_tss %>% dplyr::arrange(acrophase)
kzfp_tss$h3k9me3_znfs = rowSums(h3k9me3_all_kzfp_znfs, na.rm = T)
```

H3K4me3, H3K27Ac and ATAC signal intensity follows a bimodal distribution over KZFP promoters, whereas H3K4me1 display more gradual scales.

Each mode likely represents used vs silenced promoters.

H3K27Ac correlates highly with H3K4me3, and a little bit less with H3Kme1.

CTCF at the TSS is clearly a bimodal distribution, with the lower mode at low H3k4me3 and the higher mode at high H3K4me3.

No link between either H3K4me3 or H3K9me3 and 274 binding at the TSS, as expected.

No link between ZNF274 at the TSS and H3K9me3 at the TSS; even though there are KZFPs with high signal on both. TODO: which are they?

No link between ZNF75A/D and ZNF274 or H3K9me3 at the TSS.

Slight correlation between ZNF75A and ZNF75D tss binding.


How do H3K9me3 levels at ZNFs, and ZNF274 levels at ZNF correlate between each other, and with H3K4me3 or H3K4me1?

```R
plot(log(kzfp_tss$h3k4me3+1), log(kzfp_tss$h3k9me3_znfs+1))
```

No correlation between H3K4me3 at the most active TSS and H3K9me3 at the ZNFs. Or maybe an inverse correlation for some, but not for most of them?

```R
znf274_all_kzfp_znfs = read.table('../out/temp/kzfp_znfs_strict_sorted_acrophase_znf274.tab', sep = '\t', skip=3)
znf274_all_kzfp_znfs %>% head()
kzfp_tss$znf274_znfs = rowSums(znf274_all_kzfp_znfs, na.rm = T)
```

```R
hist(log(kzfp_tss$znf274_znfs), breaks = 100)
```

```R
plot(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs*(kzfp_tss$OutOfFrameZFAEnd-kzfp_tss$OutOfFrameZFAStart)))
```

When not normalizing by ZNF length (or multiplying the ComputeMatrix signal back by the length of the ZNF array at the DNA level), slight improvement in the negative correlation between H3K4me1 at the TSS and ZNF274 binding at the ZNF, but very slight...

```R
cor(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$h3k9me3_znfs*(kzfp_tss$OutOfFrameZFAEnd-kzfp_tss$OutOfFrameZFAStart)), method = "spearman")
```

completely uncorrelated with H3K9me3 at the ZNF

```R
plot(log(kzfp_tss$h3k9me3_znfs+1), log(kzfp_tss$atac+1))
```

```R
plot(log(kzfp_tss$znf274_znfs+1), log(kzfp_tss$h3k9me3_znfs+1))
```

```R
cor(log(kzfp_tss$znf274_znfs+1), log(kzfp_tss$h3k9me3_znfs+1), method = "spearman")
```

ZNF274 and H3K9me3 signal correlate at ZNFs, even after correcting for length!

```R
plot(kzfp_tss$CanonicalsZF, log(kzfp_tss$h3k9me3_znfs+1))
```

Is the signal already normalized for the total sequence length? It doesnt seem like it, as we can see a correlation already...

```R
cor(kzfp_tss$CanonicalsZF, log(kzfp_tss$h3k9me3_znfs+1))
```

```R
h3k36me3_all_kzfp_znfs = read.table('../out/temp/kzfp_znfs_strict_sorted_acrophase_h3k36me3.tab', sep = '\t', skip=3)
h3k36me3_all_kzfp_znfs %>% head()
```

```R
kzfp_tss$h3k36me3_znfs = rowSums(h3k36me3_all_kzfp_znfs, na.rm = T)
```

```R
plot(log(kzfp_tss$znf274_znfs+1), log(kzfp_tss$h3k36me3_znfs+1))
```

No correlation between H3K36me3 at the ZNFs and ZNF274 at the ZNFs

```R
plot(log(kzfp_tss$h3k9me3_znfs+1), log(kzfp_tss$h3k36me3_znfs+1))
```

```R
cor(log(kzfp_tss$h3k9me3_znfs+1), log(kzfp_tss$h3k36me3_znfs+1), method = "spearman")
```

On the other hand, quite a strong correlation between H3K9me3 at ZNFs and H3K36me3 at ZNFs, which is the opposite of what one expects? Or may these ressemble dual domains as described here ? https://www.sciencedirect.com/science/article/pii/S1097276521011424


Checking which out of ZNF274, ZNF75A and ZNF75D at ZNFs best correlate with H3K9me3.

```R
znf75a_all_kzfp_znfs = read.table('../out/temp/kzfp_znfs_strict_sorted_acrophase_znf75a.tab', sep = '\t', skip=3)
znf75a_all_kzfp_znfs %>% head()
```

```R
kzfp_tss$znf75a_znfs = rowSums(znf75a_all_kzfp_znfs, na.rm = T)
```

```R
plot(log(kzfp_tss$znf274_znfs+1), log(kzfp_tss$znf75a_znfs+1))
```

```R
plot(log(kzfp_tss$znf75a_znfs+1), log(kzfp_tss$h3k9me3_znfs+1))
```

```R
cor(log(kzfp_tss$znf75a_znfs+1), log(kzfp_tss$h3k9me3_znfs+1), method = "spearman")
```

Correlation is there, but less than between ZNF274 and H3K9me3. Could be explained by the difference in cell type though...

```R
znf75d_all_kzfp_znfs = read.table('../out/temp/kzfp_znfs_strict_sorted_acrophase_znf75d.tab', sep = '\t', skip=3)
znf75d_all_kzfp_znfs %>% head()
```

```R
kzfp_tss$znf75d_znfs = rowSums(znf75d_all_kzfp_znfs, na.rm = T)
```

```R
plot(log(kzfp_tss$znf75d_znfs+1), log(kzfp_tss$h3k9me3_znfs+1))
```

It looks like there are two states linked to ZNF75D: bound vs. not bound. Are they linked to ZNF274 binding itself??

```R
cor(log(kzfp_tss$znf75d_znfs+1), log(kzfp_tss$h3k9me3_znfs+1), method = "pearson")
```

Correlation is less than 274 vs K9, but still quite high.

```R
plot(log(kzfp_tss$znf75d_znfs+1), log(kzfp_tss$znf274_znfs+1))
```

If anything, ZNF274 binding increases when ZNF75D increases. No obvious inverse correlation.

```R
plot(log(kzfp_tss$znf75a_tss+1), log(kzfp_tss$h3k4me3+1))
```

```R
plot(log(kzfp_tss$znf274_tss+1), log(kzfp_tss$h3k4me3+1))
```

```R
cor(log(kzfp_tss$znf274_tss+1), log(kzfp_tss$h3k4me3+1), method = "pearson")
```

```R
cor(log(kzfp_tss$znf75a_tss+1), log(kzfp_tss$h3k4me3+1), method = "pearson")
```

```R
cor(log(kzfp_tss$znf75d_tss+1), log(kzfp_tss$h3k4me3+1), method = "pearson")
```

How can there be such a strong positive correlation between ZNF274 binding at the TSS and H3K4me3 at the TSS? Aren't they supposed to oppose each other?

```R
cor(log(kzfp_tss$h3k9me3_tss+1), log(kzfp_tss$h3k4me3+1), method = "pearson")
```

```R
plot(log(kzfp_tss$znf274_znfs+1), log(kzfp_tss$h3k4me3+1))
```

At the TSS, ZNF274 signal does not correlate with H3K9me3, but rather with H3K4me3 levels. That's quite weird??

Actually yes, but the signal for ZNF274 at the TSS is 2 orders of magnitude lower than at the ZNFs. Thus, it's likely to be quite aspecific. In addition, we don't see the signal in the other ZNF274 ENCODE ChIP-seqs, either in K562 or in H1.

```R
kzfp_tss %>% 
    dplyr::filter(log(znf274_tss+1) > 7, log(h3k4me3+1) > 7) %>% 
    dplyr::arrange(desc(znf274_tss)) %>%
    dplyr::select(symbol, acrophase, padj, phase_assigned, znf274_tss, h3k4me3, znf274_znfs)
```

```R
plot(log(kzfp_tss$h3k9me3_znfs+1), kzfp_tss$rif1)
```

Checked the top candidates by hand: signal is aspecific. Caution when interpreting the correlations: the orders of magnitude better be on the expected range of the signal, not just a minor fraction of the lower signal range.

```R
summary(kzfp_tss$znf75d_tss)
summary(kzfp_tss$znf75d_znfs)
```

```R
summary(kzfp_tss$znf274_tss)
summary(kzfp_tss$znf274_znfs)
```

```R
summary(kzfp_tss$znf75a_tss)
summary(kzfp_tss$znf75a_znfs)
```

The ranges at TSS for those ZNFs is orders of magnitude lower than the ranges observed at true peaks. Thus, even though there may be correlations, the evidence is not conclusive.

```R
htmp = (kzfp_tss %>% 
     dplyr::mutate(mesor_linear = exp(mesor), amplitude_linear = exp(amplitude), signif = -1*padj) %>% 
     dplyr::select(h3k4me3, h3k27ac, atac, h3k4me1, h3k9me3_tss,  h3k27me3, ctcf, mesor_linear, amplitude_linear, signif, h3k9me3_tss, h3k9me3_znfs, znf274_znfs, h3k36me3_znfs, znf75a_znfs, znf75d_znfs) + 1) %>% 
    dplyr::rename(H3K4me3 = h3k4me3, H3K27ac = h3k27ac, ATAC = atac, H3K4me1 = h3k4me1, `prom. H3K9me3` = h3k9me3_tss, H3K27me3 = h3k27me3, CTCF = ctcf, 
                  MESOR = mesor_linear, amplitude = amplitude_linear, `signif.`=signif,
                 `C2H2 ZNF274`= znf274_znfs, `C2H2 H3K9me3` = h3k9me3_znfs, `C2H2 H3K36me3`= h3k36me3_znfs, `C2H2 ZNF75A` = znf75a_znfs, `C2H2 ZNF75D` = znf75d_znfs) %>%
    log() %>% cor(method = "spearman") %>% ComplexHeatmap::Heatmap(name = "SCC", col = colorRamp2(c(-1, 0, 1), rev(c('#f1a340','#f7f7f7','#998ec3'))))

pdf('../out/cell_cycle_figures/promoter_state_vs_cosinor_params_heatmap_kzfps.pdf', width = 4.5, height = 4)
draw(htmp)
dev.off()
```

### Verifying whether the computeMatrix scale-region output normalizes for sequence length

```R
plot(kzfp_tss$CanonicalsZF, log(kzfp_tss$znf274_znfs+1))
```

```R
plot(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$h3k9me3_znfs+1))
```

```R
cor(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$h3k9me3_znfs+1), method = "spearman")
```

There may be a very slight inverse correlation between H3K4me1 at the TSS and H3K9me3 at the ZNFs, corrected for length (because computed using scale-region). But what if we did not correct for ZNF length?

Actually the correlation coefficient is slightly positive, or very close to zero. Nice biopeotry

```R
plot(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs+1))
```

```R
plot(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs*(kzfp_tss$OutOfFrameZFAEnd-kzfp_tss$OutOfFrameZFAStart)+1))
```

```R
plot(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs+1)*(kzfp_tss$OutOfFrameZFAEnd-kzfp_tss$OutOfFrameZFAStart))
```

```R
plot(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs+1)*kzfp_tss$CanonicalsZF)
```

```R
cor(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs+1), method = "spearman")
```

```R
cor(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs*(kzfp_tss$OutOfFrameZFAEnd-kzfp_tss$OutOfFrameZFAStart)+1), method = "spearman")
```

```R
cor(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs+1)*(kzfp_tss$OutOfFrameZFAEnd-kzfp_tss$OutOfFrameZFAStart), method = "spearman")
```

```R
cor(log(kzfp_tss$h3k4me1+1), log(kzfp_tss$znf274_znfs+1)*kzfp_tss$CanonicalsZF, method = "spearman")
```

Maybe a slight inverse correlation between H3K4me1 at the promoter and ZNF274 binding at the ZNF, corrected for ZNF length.

```R
znf274_all_kzfp_znfs_bigWigAverageOverBed = read.table('../out/temp/kzfp_znfs_strict_sorted_acrophase_ZNF274_bigWigAverageOverBed.tab', sep = '\t', header = F, 
                                                       col.names = c("name", "size", "covered", "sum", "mean0", "mean"))
znf274_all_kzfp_znfs_bigWigAverageOverBed %>% dim()
```

```R
znf274_all_kzfp_znfs_bigWigAverageOverBed %>% head()
```

```R
kzfp_tss$symbol %>% head()
```

```R
kzfp_tss %>% head()
```

```R
kzfp_tss = kzfp_tss %>% dplyr::left_join(., znf274_all_kzfp_znfs_bigWigAverageOverBed %>% dplyr::select(name, sum), by = c("symbol"="name"))
```

```R
plot(log(kzfp_tss$sum+1), log(kzfp_tss$znf274_znfs+1))
```

There is still a correlation between the total ZNF274 at ZNFs and the one returned by computeMatrix...

```R
plot(kzfp_tss$CanonicalsZF, log(kzfp_tss$sum+1))
```

```R
cor(kzfp_tss$CanonicalsZF, log(kzfp_tss$sum+1), method = "spearman")
```

There is a mild but noticeable correlation between the number of ZNFs per KZFP and the total ZNF274 signal. This suggests that the computeMatrix function averages things over the window length. What lends creedance to that claim is that `my_average` in the `heatmapper` class of `deeptools computeMatrix` is called to average the coverage per base.

```R
plot(kzfp_tss$sum/znf274_all_kzfp_znfs_bigWigAverageOverBed$size, kzfp_tss$znf274_znfs)
```

```R
cor(kzfp_tss$sum/znf274_all_kzfp_znfs_bigWigAverageOverBed$size, kzfp_tss$znf274_znfs, method = "spearman")
```

Indeed, dividing the total ZNF274 signal at the bigwig by the length of each ZNF yields an equivalent to the output of `computeMatrix scale-region`.

But since ZNF274 signal at ZNFs still correlates with the number of ZNFs, this suggests that ChIP-seq signal is a function of a square of the number of ZNFs??

Or is it that many KZFPs have similar numbers of ZNFs, such that the signal correlates for a majority of them?

```R
kzfp_tss$CanonicalsZF %>% hist(breaks = 100)
```

```R
kzfp_tss %>% ggplot(aes(x = CanonicalsZF, y = sum)) + geom_point()
```

```R
kzfp_tss %>% colnames()
```

```R
kzfp_tss %>% ggplot(aes(x = CanonicalsZF, y = log(sum))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = OutOfFrameZFAEnd-OutOfFrameZFAStart, y = log(sum))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = OutOfFrameZFAEnd-OutOfFrameZFAStart, y = log(sum/(OutOfFrameZFAEnd-OutOfFrameZFAStart)))) + geom_point()
```

```R
cor((kzfp_tss$OutOfFrameZFAEnd-kzfp_tss$OutOfFrameZFAStart), log(kzfp_tss$sum/(kzfp_tss$OutOfFrameZFAEnd-kzfp_tss$OutOfFrameZFAStart)), method="spearman")
```

```R
kzfp_tss %>% ggplot(aes(x = CanonicalsZF, y = log(sum))) + geom_point()
```

```R
kzfp_tss %>% group_by(CanonicalsZF) %>% summarize(log(mean(sum))) %>% cor(., method="spearman")
```

```R
kzfp_tss %>% group_by(CanonicalsZF) %>% summarize(log(mean(znf274_znfs))) %>% cor(., method="spearman")
```

```R
kzfp_tss %>% group_by(CanonicalsZF) %>% summarize(log(median(sum))) %>% cor(., method="spearman")
```

```R
kzfp_tss %>% group_by(CanonicalsZF) %>% summarize(log(median(znf274_znfs))) %>% cor(., method="spearman")
```

Still a pretty big correlation between the number of ZNFs and the median ZNF274 signal when grouping by number of znfs. Dividing by the length once does not do away the trend completely.

```R
kzfp_tss %>% group_by(CanonicalsZF) %>% ggplot(aes(x = CanonicalsZF, y = log(sum), group = CanonicalsZF)) + geom_boxplot()
```

```R
kzfp_tss %>% group_by(CanonicalsZF) %>% ggplot(aes(x = CanonicalsZF, y = log(znf274_znfs), group = CanonicalsZF)) + geom_boxplot()
```

Still, maybe it's a bias due to the mapping that is also related to the square of sequences...

```R
kzfp_tss %>% ggplot(aes(x = log(sum), y = log(znf274_znfs*(OutOfFrameZFAEnd-OutOfFrameZFAStart)))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = CanonicalsZF, y = OutOfFrameZFAEnd-OutOfFrameZFAStart)) + geom_point()
```

### Analyzing rhythmic genes and KZFPs

```R
fit_kzfps %>% dplyr::filter(padj < 0.05, grepl('S', phase_assigned)) 
```

```R
fit_kzfps %>% dplyr::filter(padj < 0.05) %>% dim()
```

```R
fit_kzfps %>% dplyr::filter(padj < 0.05) %>% arrange(pvalue) %>% tail(15)
```

Exporting the data for GO enrichment analysis (done by Juna)


Are genes distributed around the phases or are they all the same ?

```R
fit_qualitative %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase, pvalue) %>% ggplot() + geom_histogram(mapping = aes(acrophase), bins = 50)
```

As a reminder: each phase is centered on an integer from 0 to 7, starting with G1. 

Thus, there is an elevated number of genes with phase value of about 2.5, which means they belong somewhere in S phase, between S1 and S2 (timepoints defined by Romain). 

To export the genes for GO, we split significantly oscillating genes according to each phase, each one being centered on an integer (thus the bins are integers +- 0.5).

This is equivalent to rounding to the closest integer between 0 and 7 (keeping in mind that 0 and 8 are the same, which we obtain by doing 8 modulo 8)

```R
fit_qualitative %>% dplyr::select(symbol, mesor, amplitude, acrophase, rsquared, pvalue, padj, phase_assigned) %>% dplyr::filter(padj < 0.05) %>% arrange(pvalue) %>% head(20)
```

CHAF1B is required for the assembly of histone octamers onto newly-replicated DNA (GeneCards)


```R
cell_cycle_genes_romain %in% ((fit_qualitative %>% dplyr::filter(padj < 0.05)))$symbol %>% sum()
```

```R
cell_cycle_genes_romain[which(!cell_cycle_genes_romain %in% ((fit_qualitative %>% dplyr::filter(padj < 0.05)))$symbol)]
```

```R
fit_qualitative %>% dplyr::filter(symbol == 'TTK')
```

TTK is required for centrosome duplication, thus its expression is expected to peak late in the cell cycle process. It is the only of the 168 cell cycle marker genes suggested by Romain that fails (barely) the rhythmicity test.


## GO enrichment analysis of rhythmic genes

```R
library(clusterProfiler)
library(org.Hs.eg.db)
```

```R
fit_qualitative %>% colnames()
```

```R
GOcluster_per_phase <- clusterProfiler::compareCluster(symbol~phase_assigned, 
                            data = fit_qualitative %>% dplyr::filter(padj < 0.05) %>% dplyr::select(symbol, phase_assigned), 
                            ont = "BP", 
                            keyType = 'SYMBOL',
                            universe = fit_qualitative %>% dplyr::select(symbol) %>% unique() %>% pull(), 
                            OrgDb = org.Hs.eg.db, 
                            fun="enrichGO",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.1,
                            minGSSize = 10,
                            maxGSSize = 100
)
```

```R
p = clusterProfiler::dotplot(GOcluster_per_phase, showCategory = Inf) + 
scale_colour_gradientn(colours = rev(viridis(100)), limits=c(0, 0.05))
p
```

```R
ggsave('../out/cell_cycle_figures/rhythmic_genes_enrichment_per_phase_shortened.svg', p, svg, height = 10, width = 6.5)
```

```R
# exporting the entire GO enrichment table
clusterProfiler::dotplot(GOcluster_per_phase, showCategory = Inf)$data %>% write.table('../out/tables/rhythmic_genes_enrichment_GOBP.tsv', row.names = F, col.names = T, sep = '\t')
```

<!-- #raw -->
# this is extremely slow, probably because 3200 genes is way too much for a GO analysis...
GOcluster_all_rhythmic <- clusterProfiler::compareCluster(fit_qualitative %>% dplyr::filter(padj < 0.05) %>% dplyr::pull(symbol),
                            ont = "BP", 
                            keyType = 'SYMBOL',
                            universe = fit_qualitative %>% dplyr::select(symbol) %>% unique() %>% pull(), 
                            OrgDb = org.Hs.eg.db, 
                            fun="enrichGO",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.1,
                            minGSSize = 10,
                            maxGSSize = 100
)
<!-- #endraw -->

Adding genomic region binding enrichment

```R
kzfps_binding_region_enrich = read.table('../data/KZFPs_vs_genomicPos.tsv', sep = '\t', header = T)
colnames(kzfps_binding_region_enrich)[1] = "symbol"
kzfps_binding_region_enrich %>% head()
```

```R
fit_kzfps = fit_kzfps %>% dplyr::left_join(., kzfps_binding_region_enrich)
fit_kzfps %>% head()
```

```R
fit_kzfps %>% dplyr::filter(!is.na(promoter.pval)) %>% dplyr::arrange(promoter.pval) %>% dplyr::select(symbol, promoter.pval, padj, phase_assigned) %>% head(30)
```

```R
fit_kzfps %>% dplyr::filter(!is.na(promoter.pval)) %>% dplyr::arrange(promoter.pval) %>% dplyr::select(symbol, promoter.pval, padj, phase_assigned) %>% tail(30)
```

```R
fit_kzfps %>% dplyr::filter(!is.na(exons.pval)) %>% dplyr::arrange(exons.pval) %>% dplyr::select(symbol, exons.pval, padj, phase_assigned) %>% head(30)
```

```R
fit_kzfps %>% dplyr::filter(!is.na(exons.pval)) %>% dplyr::arrange(exons.pval) %>% dplyr::select(symbol, exons.pval, padj, phase_assigned) %>% tail(30)
```

There seems to be an enrichment for promoter and exon binding KZFPs among those whose expression peaks in the S phase.


Adding the TSS and ZNF histone mark scores for KZFPs

```R
kzfp_tss %>% colnames()
```

```R
kzfp_tss %>% dim()
```

```R
fit_kzfps = fit_kzfps %>% dplyr::left_join(., kzfp_tss %>% dplyr::select(ensembl, h3k4me3, h3k27ac, atac, h3k4me1, h3k27me3, h3k9me3_tss, ctcf, znf274_tss, znf75a_tss, znf75d_tss, h3k9me3_znfs, znf274_znfs, h3k36me3_znfs, znf75a_znfs, znf75d_znfs))
```

```R
fit_qualitative %>% write.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', row.names = F)
fit_kzfps %>% write.table('../out/tables/kzfps_cycling_info.tsv', sep = '\t', row.names = F)
```

Getting average corrected expression values for future plotting

```R
norm_adj_counts = expr_corrected %>% dplyr::select(-one_of(c('ensembl', 'symbol', 'genename'))) %>% as.matrix()
rownames(norm_adj_counts) = expr_corrected$ensembl
norm_adj_counts %>% head()
```

```R
norm_adj_counts_no_outliers = norm_adj_counts[, !colnames(norm_adj_counts) %in% to_exclude]
```

```R
expr_plot = as.data.frame(log2(norm_adj_counts_no_outliers + 1))
expr_plot %>% head()
```

```R
# we use the mean and not the median as we consider the replicates as normally distributed, plus there are three replicates -> that each participates with 1/3 seems justified
expr_means = list()
for (p in names(phase)) {
    expr_means = cbind(expr_means, expr_plot %>% as.data.frame() %>% dplyr::select(dplyr::one_of(phase[[p]])) %>% apply(., 1, mean) %>% as.vector())
    }
expr_means = as.data.frame(expr_means) %>% mutate(across(where(is.list), as.double))
rownames(expr_means) = rownames(expr_plot)
colnames(expr_means) = names(phase)
expr_means %>% head()
```

```R
expr_zscore = expr_means %>% apply(., MARGIN = 1, FUN = scale) %>% t()
colnames(expr_zscore) = colnames(expr_means)
expr_zscore %>% head()
```

```R
colnames(expr_zscore) = new_phases
expr_zscore %>% head()
```

```R
expr_zscore %>% write.table('../out/tables/expr_zscore_RNAseq_cellcycle_romain.tsv', col.names = T, row.names = T, sep = '\t')
```

TODO we should really stop here, and even not do the enrichment plots beforehand.


### Computing the H3K9me3 signal around the TSS of rhythmic genes

```R
# exporting the coordinates of rhythmic genes TSS sorted by acrophase as a bed file
fit_qualitative %>% 
    dplyr::filter(padj < 0.05) %>% 
    arrange(acrophase) %>% 
    dplyr::left_join(tss_top_h3k4me3_h3k27ac_all_genes %>% 
                     dplyr::select(symbol, TSSstart, TSSend)) %>%
    dplyr::select(chr, TSSstart, TSSend, symbol) %>%
    write.table('../out/temp/tss_rhythmic_genes_sorted_acrophase.bed', col.names = F, row.names = F, sep = '\t', quote = F)
```

```R
# Computing and adding the width or intensity of the H3K9me3 trough at the promoter of rhythmic genes
h3k9me3_promoters_hires = read.table('../out/temp/tss_rhythmic_genes_sorted_acrophase_h3k9me3_domains.tab', sep = '\t', skip = 3)
```

```R
h3k9me3_promoters_hires %>% head()
```

```R
# plotting the raw signal
log(h3k9me3_promoters_hires+1) %>% as.matrix() %>% ComplexHeatmap::Heatmap(cluster_rows = F, 
                                                                    cluster_columns = F,
                                                                    col= colorRamp2(c(0, 0.01, 0.5, 1), viridis::magma(4)))
```

```R
# let's attempt smoothing the signal to observe a more general trend
library(spatstat)
```

```R
library(EBImage)
median_kernel_size = 4
gaussian_kernel_size = 7
nrow_padding = max(median_kernel_size, gaussian_kernel_size)
mat = log(h3k9me3_promoters_hires+1) %>% as.matrix()
# we pad by closing the circle, i.e. adding data from the other size of the circle
mat_padded = rbind((mat %>% tail(nrow_padding)), mat, (mat %>% head(nrow_padding)))
im.med = medianFilter(Image(mat_padded), size=4) %>% as.matrix()
im.med%>% ComplexHeatmap::Heatmap(cluster_rows = F, 
                                                                    cluster_columns = F,
                                                                    col= colorRamp2(c(0, 0.01, 0.5, 1), viridis::magma(4)))
```

```R
mat_smoothed = spatstat.explore::blur(as.im(im.med@.Data), sigma = 7, kernel = "gaussian", normalise = F, bleed = F) %>% as.matrix() %>%.[(nrow_padding+1):(nrow(mat)+nrow_padding), ]
mat_smoothed %>% ComplexHeatmap::Heatmap(cluster_rows = F, 
                                        cluster_columns = F,
                                        col= colorRamp2(c(0, 0.01, 0.5, 1), viridis::magma(4)))
```

```R
cutoff = median(summary(mat_smoothed[, 85:215] %>% as.vector()))
cutoff
```

```R
h3k9me3_signal_tss = mat_smoothed[, 85:215] %>% apply(., 1, function(x) sum(x > cutoff))
```

```R
plot(1:nrow(mat_smoothed), h3k9me3_signal_tss)
```

```R
im.med@.Data %>% as.vector() %>% summary()
```

```R
plot(1:nrow(im.med@.Data), im.med@.Data[, 85:215] %>% apply(., 1, function(x) sum(x > 0.75)))
```

```R
n_rhythmic = fit_qualitative %>% dplyr::filter(padj < 0.05) %>% nrow()
```

```R
paste0(n_rhythmic, ' rhythmic genes (padj < 0.05)')
```

```R
plot_order = (fit_qualitative %>% dplyr::filter(padj < 0.05) %>% dplyr::arrange(acrophase))$ensembl
library(circlize)
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "expr z-scored log2", direction = "horizontal")
```

```R
library("viridis")  
```

```R
fit_qualitative$RepliTiming %>% unique()
```

```R
replitiming_colors = c(viridis(length(fit_qualitative$RepliTiming %>% unique())))
replitiming_colors
```

```R
fit_qualitative$RepliTiming %>% unique() %>% sort() %>% as.character()
```

```R
fit_qualitative %>% dplyr::filter(RepliTiming == "NA") %>% dim()
```

```R
names(replitiming_colors) = fit_qualitative$RepliTiming %>% unique() %>% sort() %>% as.character()
# there's no NAs 
#replitiming_colors = replitiming_colors[-1]
```

```R
fit_qualitative$rif1 %>% summary()
```

```R
fit_qualitative$RTindex %>% summary()
```

```R
cols = list(`Repl. timing` = replitiming_colors, 
            MESOR = colorRamp2(breaks = c(-2.5, 2.5, 4, 6, 10), colors = viridis::mako(5)),
             RIF1 =  colorRamp2(c(-0.5, -0.2, 0, 0.2, 0.5), hcl.colors(n = 5, palette = "Cividis", rev = F)),
           RTindex = colorRamp2(breaks = seq(from = 2.5, to = 4.5, length.out = 4), colors = viridis::viridis(4)),
            H3K9me3 = colorRamp2(breaks = seq(from = 70, to = 145, length.out = 10), colors = viridis::magma(10)))
```

```R
plot_selection = fit_qualitative %>% dplyr::filter(padj < 0.05, ensembl %in% tss$ensembl) %>% dplyr::arrange(acrophase)
```

```R
# legacy: RepliTiming was added to row_anno in the past `Repl. timing` = RepliTiming,

row_anno = plot_selection %>% dplyr::select(RTindex_avg, rif1, mesor) %>% 
                        dplyr::rename(RIF1 = rif1, MESOR = mesor, RTindex = RTindex_avg)
row_anno$H3K9me3 = h3k9me3_signal_tss
```

```R
row_anno %>% head()
```

```R
# annotating the top 10 rhythmic gene per phase merged with the top 10 rhythmic genes overall
genes_to_label = c(fit_qualitative %>% 
    dplyr::filter(padj < 0.05, !phase_assigned %in% c()) %>% 
    dplyr::group_by(phase_assigned) %>%
    dplyr::arrange(padj) %>% dplyr::slice_head(n=1) %>% dplyr::pull(symbol),
  fit_qualitative %>%
      dplyr::arrange(padj) %>% head(10) %>% pull(symbol) %>% unique())
```

```R
genes_to_label_idx = match(genes_to_label, (fit_qualitative %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase) %>% pull(symbol)))
```

```R
genes_to_label_idx
```

```R
ha = rowAnnotation(df = row_anno, col = cols)
label_annotation = ComplexHeatmap::rowAnnotation(foo = anno_mark(at = genes_to_label_idx, side = "left", padding = 2,
    labels = genes_to_label))
plot_order = plot_selection$ensembl
```

```R
expr_color = colorRamp2(breaks = c(-2, -1, 0, 1, 2), colors = hcl.colors(5, palette = "RdBu", rev = T))
```

```R
fit_qualitative %>% arrange(padj) %>% head(10)
```

```R
#ha = rowAnnotation(RepliTiming = anno_points((fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'RepliTiming']))
cn = sapply(1:length((expr_zscore %>% colnames())), FUN = function(x) paste0((expr_zscore %>% colnames())[x], '\n', as.roman(x)))

htmp = Heatmap(expr_zscore[plot_order, ],
        show_row_names = F, 
        row_title = paste0(plot_order %>% length(), ' rhythmic genes (adj. p < 0.05)'), 
        show_column_names = F, #done via an annotation due to the automatic rotation
        heatmap_legend_param = list(
            title = "expr.",
            direction = "horizontal",
        title_position = "topcenter"),
        cluster_rows = FALSE, 
               cluster_columns = FALSE,
              right_annotation = ha,
               left_annotation = label_annotation,
              col = expr_color,
              bottom_annotation = HeatmapAnnotation(
        text = anno_text(cn, rot = 0, location = unit(1, "npc"), just = "top"),
        annotation_height = max_text_width(cn)
    )) %v% NULL
pdf('../out/cell_cycle_figures/cycling_genes_heatmap_with_RepliTiming.pdf', height = 5.5, width = 6)
draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
```

```R
svg('../out/cell_cycle_figures/cycling_genes_heatmap_with_RepliTiming.svg', height = 5.5, width = 6)
draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
```

```R
tss_top_h3k4me3_h3k27ac_all_genes %>% colnames()
```

Computing stacked barplots and enrichments for phases vs. Repli-seq and RIF1

| | rhythmic genes repl. in G1 | rhythmic genes repl. in S1-G2 | total |
|---|:---:|:---:|:---:|
| rhythmic genes in S2-M | x | m-x | m |
| rhythmic genes in eG1-S1 | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
# computing enrichment tests
# S2-M vs the rest
k = plot_selection %>% dplyr::filter(RepliTiming == "G1") %>% nrow() # number of rhythmic genes repl. in G1

acr = c("S2", "G2", "G2/M", "M")

m = plot_selection %>% dplyr::filter(phase_assigned %in% acr) %>% nrow() # number of rhythmic genes in S2-M

n = plot_selection %>% dplyr::filter(!phase_assigned %in% acr) %>% nrow() # number of rhythmic genes in eG1-S1

x = 0:m # is the variable tested

x_observed = plot_selection %>% dplyr::filter(RepliTiming == "G1", phase_assigned %in% acr) %>% nrow() # number of rhythmic genes in S-M replicated in G1,

probs <- dhyper(x, m, n, k, log = FALSE)

# we compute the probability of observing a more extreme depletion, therefore using a one sided test. 

pval_one_sided = sum(probs[x>=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])
pval_two_sided
x_observed
k
```

| | rhythmic genes repl. in S4/G2 | rhythmic genes repl. in G1-S3 | total |
|---|:---:|:---:|:---:|
| rhythmic genes in S2-M | x | m-x | m |
| rhythmic genes in eG1-S1 | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
# computing enrichment tests
# S2-M vs the rest
timing = c("S4", "G2")
k = plot_selection %>% dplyr::filter(RepliTiming %in% timing) %>% nrow() # number of rhythmic genes repl. in G1

acr = c("S2", "G2", "G2/M", "M")

m = plot_selection %>% dplyr::filter(phase_assigned %in% acr) %>% nrow() # number of rhythmic genes in S2-M

n = plot_selection %>% dplyr::filter(!phase_assigned %in% acr) %>% nrow() # number of rhythmic genes in eG1-S1

x = 0:m # is the variable tested

x_observed = plot_selection %>% dplyr::filter(RepliTiming %in% timing, phase_assigned %in% acr) %>% nrow() # number of rhythmic genes in S-M replicated in G1,

probs <- dhyper(x, m, n, k, log = FALSE)

# we compute the probability of observing a more extreme depletion, therefore using a one sided test. 

pval_one_sided = sum(probs[x<=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])
pval_two_sided
x_observed
k
```

Depletion of late replicating genes among S2-M genes.

```R
pdf('../out/cell_cycle_figures/prop_genes_in_repliTiming_per_phase.pdf', height = 2.5, width = 4.5)
p = plot_selection %>% 
    dplyr::select(phase_assigned, RepliTiming) %>% 
    dplyr::group_by(phase_assigned) %>% 
    dplyr::mutate(genes_in_phase = n()) %>% 
    dplyr::group_by(phase_assigned, RepliTiming) %>% 
    reframe(perc_in_repliTiming = n()/genes_in_phase) %>% 
    unique() %>% 
    dplyr::rename(`Repl. timing` = RepliTiming) %>%
    ggplot(aes(x = phase_assigned, y = perc_in_repliTiming, fill = `Repl. timing`)) + 
        geom_col() + 
        scale_fill_manual(values= replitiming_colors) + 
        theme_classic() + 
        xlab("") + ylab("Gene prop.") + 
        theme(axis.ticks.x=element_blank(), axis.line = element_blank()) + 
        scale_y_continuous(limits = c(0,1), expand = c(0, 0))
p
dev.off()
```

```R
svg('../out/cell_cycle_figures/prop_genes_in_repliTiming_per_phase.svg', height = 2.5, width = 4.5)
p
dev.off()
```

T-test for RIF1 levels, between S2-M and eG1-S1

```R
pdf('../out/cell_cycle_figures/RIF1_per_phase.pdf', height = 2.5, width = 4.5)
p = plot_selection %>% 
    dplyr::select(phase_assigned, rif1) %>% 
    dplyr::rename(RIF1 = rif1) %>% 
    ggplot(aes(x = phase_assigned, y = RIF1)) + 
    geom_hline(yintercept = 0) +
    geom_violin() + 
    #geom_jitter(width = 0.2, height = 0, size = 0.3, alpha = 0.3) + 
    geom_boxplot(width=0.1, outlier.size = 0.5) + 
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("RIF1 [log2 FC over input]") + 
    scale_y_continuous(limits = c(-2,2), expand = c(0, 0))

p
dev.off()


```

```R
svg('../out/cell_cycle_figures/RIF1_per_phase.svg', height = 2.5, width = 4.5)
p
dev.off()


```

```R
acr = c("G2", "G2/M")
acr_excluded = c("S1", "S2", "M")

y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(rif1) %>% 
    median(na.rm = T)

p = plot_selection %>%
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::select(phase_assigned, rif1) %>% 
    dplyr::rename(RIF1 = rif1) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "RIF1", add = "boxplot") + 
    geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 2.4, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("RIF1 [log2 FC over input]") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("G1", "G2")) + 
    coord_cartesian(ylim = c(-2, 2.5),
                      clip = 'off')

svg('../out/cell_cycle_figures/RIF1_stat_test.svg', height = 2.5, width = 2)
p
dev.off()
```

```R
pdf('../out/cell_cycle_figures/MESOR_per_phase.pdf', height = 2.5, width = 4.5)
p = plot_selection %>% 
    dplyr::select(phase_assigned, mesor) %>% 
    dplyr::rename(MESOR = mesor) %>% 
    ggplot(aes(x = phase_assigned, y = MESOR)) + 
    geom_violin() + 
#    geom_jitter(width = 0.2, height = 0, size = 0.3, alpha = 0.3) + 
    geom_boxplot(width=0.1, outlier.size = 0.5) + 
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("MESOR") + 
    scale_y_continuous(limits = c(-5,10), expand = c(0, 0))
p
dev.off()

```

```R
svg('../out/cell_cycle_figures/MESOR_per_phase.svg', height = 2.5, width = 4.5)
p
dev.off()
```

```R
acr = c("G2", "G2/M")
acr_excluded = c("S1", "S2", "M")

y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(mesor) %>% 
    median(na.rm = T)

p = plot_selection %>%
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::select(phase_assigned, mesor) %>% 
    dplyr::rename(MESOR = mesor) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "MESOR", add = "boxplot") + 
    geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 14, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("MESOR") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("G1", "G2")) + 
    coord_cartesian(ylim = c(-6, 15),
                      clip = 'off')

svg('../out/cell_cycle_figures/MESOR_stat_test.svg', height = 2.5, width = 2)
p
dev.off()
```

```R
acr = c("G2", "G2/M")
acr_excluded = c("S1", "S2", "M")
y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(RTindex_avg) %>% 
    median(na.rm = T)

p = plot_selection %>%
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::select(phase_assigned, RTindex_avg) %>% 
    dplyr::rename(RTindex = RTindex_avg) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "RTindex", add = "boxplot") + 
    geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 7.5, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("RTindex") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("G1", "G2")) + 
    coord_cartesian(ylim = c(-6, 8),
                      clip = 'off')

svg('../out/cell_cycle_figures/RTindex_stat_test.svg', height = 2.5, width = 2)
p
dev.off()
```

<!-- #raw -->
y_dashed = plot_selection %>% 
    dplyr::left_join(., tss_top_h3k4me3_h3k27ac_all_genes %>% 
                     dplyr::select(symbol, h3k9me3_domains)) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(h3k9me3_domains) %>% 
    median(na.rm = T)
<!-- #endraw -->

```R

```

```R
h3k9me3_promoters_sum = rowMeans(log(h3k9me3_promoters_hires[, 85:215]+1))
```

```R
h3k9me3_promoters_sum %>% length()
```

```R
acr = c("S1", "S2", "G2")
acr_excluded = c("G1/S", "G2/M")
y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::left_join(., tss_top_h3k4me3_h3k27ac_all_genes %>% 
                         dplyr::filter(padj < 0.05) %>%
                     dplyr::arrange(acrophase) %>%
                     dplyr::mutate(h3k9me3_promoters = h3k9me3_promoters_sum)) %>%
    dplyr::pull(h3k9me3_promoters) %>% 
    median(na.rm = T)

p = plot_selection %>% 
    dplyr::left_join(., tss_top_h3k4me3_h3k27ac_all_genes %>% 
                         dplyr::filter(padj < 0.05) %>%
                     dplyr::arrange(acrophase) %>%
                     dplyr::mutate(h3k9me3_promoters = h3k9me3_promoters_sum) %>%
                     dplyr::select(symbol, h3k9me3_promoters)) %>%
    dplyr::select(phase_assigned, h3k9me3_promoters) %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(H3K9me3 = h3k9me3_promoters) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "H3K9me3", add = "boxplot") + 
    geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 3.2, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("H3K9me3") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-G1", "S1-to-G2")) + 
    coord_cartesian(ylim = c(0, 3.5),
                      clip = 'off')

svg('../out/cell_cycle_figures/H3K9me3_domains_stat_test.svg', height = 2.5, width = 2)
p
dev.off()
```

```R
svg('../out/cell_cycle_figures/RTindex_per_phase.svg', height = 2.5, width = 4.5)
p = plot_selection %>% 
    dplyr::select(phase_assigned, RTindex_avg) %>% 
    dplyr::rename(RTindex = RTindex_avg) %>% 
    ggplot(aes(x = phase_assigned, y = RTindex)) + 
    geom_violin() + 
#    geom_jitter(width = 0.2, height = 0, size = 0.3, alpha = 0.3) + 
    geom_boxplot(width=0.1, outlier.size = 0.5) + 
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("RTindex") + 
    scale_y_continuous(limits = c(-6,6), expand = c(0, 0))
p
dev.off()

```

TODO check whether that's really necessary

```R
# are genes peaking in G1 less expressed than those peaking in s ?
expr_means_per_gene = expr_means %>% apply(., MARGIN = 1, FUN = median)
expr_means_per_gene %>% head()
```

```R
expr_means_per_gene[plot_order] %>% head()
```

```R
fit_qualitative_batch %>% head()
```

```R
expr_zscore[plot_order, ] %>% dim()
```

```R
plot_order %>% length()
ha
```

Is there a bias for highly expressed genes in cycling genes ?

```R
fit_qualitative %>% dplyr::mutate(rhythmic = padj < 0.05) %>% ggpubr::ggviolin(x = "rhythmic" , y = "mesor", add = "boxplot") + stat_compare_means(label.x = 0.9, label.y = 11, step.increase = 1)


```

Genes that are cycling tend to be more highly expressed than genes that do not cycle. 

```R
plot(fit_qualitative$mesor, -log10(fit_qualitative$padj))
abline(h = -log10(0.05))
```

```R
expr_means %>% colSums()
```

```R
plot(expr_means %>% colSums())
```

```R
plot((fit_kzfps %>% dplyr::filter(padj < 0.05))$acrophase, (fit_kzfps %>% dplyr::filter(padj < 0.05))$mesor)
```

```R
plot((fit_kzfps)$acrophase, (fit_kzfps)$mesor)
```

If the level of expression ever becomes a concern, we could: correct by the library size of each time point (unlikely to help), look for a set of internal control genes which do not vary, or even better: are highly expressed and anti-correlate with the trend in MESOR.

Or devise a permutation test for the distribution of KZFPs on the phase map, taking into account the expression level.

```R
group_cols = viridis(length(fit_qualitative$RepliTiming %>% unique()))
group_cols
```

```R
fit_qualitative$RepliTiming %>% unique()
```

```R
names(group_cols) = (fit_qualitative$RepliTiming %>% unique())
```

```R
fit_kzfps$cluster %>% levels()
```

```R
cluster_cols = viridis(length(fit_kzfps$cluster %>% levels()))
names(cluster_cols) = fit_kzfps$cluster %>% levels()
```

```R
cluster_cols["noCluster"] = "grey"
```

```R
cols = list(RepliTiming = group_cols, 
            mesor = colorRamp2(breaks = c(-2.5, 0, 2.5, 5, 10), colors = viridis::magma(5)), 
            cluster = cluster_cols,
           rif1 = colorRamp2(c(-0.5, 0, 0.5), rev(c('#f1a340','#f7f7f7','#998ec3'))))
```

```R
plot_selection = fit_kzfps %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase)
```

```R
row_anno = plot_selection %>% dplyr::select(RepliTiming, rif1, mesor, cluster)
ha = rowAnnotation(df = row_anno, col = cols)
```

<!-- #raw -->
row_anno = fit_kzfps %>% dplyr::filter(padj < 0.05) %>% arrange(TEs.pval, acrophase) %>% dplyr::mutate(signif.rhythm = -log10(padj),
                                                             signif.promoter_bind = -log10(promoter.pval),
                                                              signif.exons_bind = -log10(exons.pval),
                                                              signif.TEs_bind = -log10(TEs.pval)) %>% dplyr::select(RepliTiming, mesor, signif.rhythm, cluster,
                                                                                                                  signif.promoter_bind, signif.exons_bind, signif.TEs_bind)
ha = rowAnnotation(df = row_anno, col = c(cols, 
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.promoter_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.exons_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.TEs_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5))))
<!-- #endraw -->

```R
plot_order = (fit_kzfps %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
n_rhythmic = nrow(to_plot)
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases
htmp = Heatmap(to_plot,
        name = 'expr, z-scored log2', 
        row_title = paste0(n_rhythmic, ' rhythmic KZFPs (padj < 0.05)'),
        heatmap_legend_param = list(
            title = "expr., z-scored log2",
            direction = "horizontal",
            title_position = "topcenter"),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
       row_names_gp = grid::gpar(fontsize = 7), right_annotation = ha) %v% NULL
pdf('../out/cell_cycle_figures/cycling_KZFPs_heatmap.pdf', height = 10, width = 6)
draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
```

Note: arranging by promoter.pval or TEs.pval does not suggest any enrichment in G1-restricted or S-restricted KZFPs.

<!-- #region -->
S phase KZFPs:

We should cross them with the mass spec interactors as well !

ZNF565: binds MLTs and Alus, signif enriched on TEs, slightly enriched on promoters.


ZNF441 binds Alus and  is significantly downregulated in rCC (https://cgp.iiarjournals.org/content/cgp/19/3/305.full.pdf)

ZNF627 has no known functions and binds Tigger1

ZNF445 binds Alus and is involved in imprinting with ZNF57 in mice (https://cgp.iiarjournals.org/content/cgp/19/3/305.full.pdf)

ZNF274 is known to repress physically clustered genes, among which other KZFPs (Martina). Could it explain the downregulation of all the other KZFPs in S phase ? Mass spec (Martina): TurboID and found lot of DNA damage and repair stuff, MCM6, but these things also interact with KAP1 

ZNF764 associates with TFIIIC (https://www.embopress.org/doi/full/10.15252/embj.2018101220) and is involved in chromosome conformation 

ZKSCAN2 ??

ZNF331 targets MLT1G and H and is a putative tumor suppressor: https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-017-0417-4

ZNF689 has no targets other than SVA_F, is anti-apoptotic in cancer and blocks cells in subG1: https://www.sciencedirect.com/science/article/pii/S001448271100187X

ZNF3 (G2) interacts with MCM6 and chromsome assembly machinery (PY interactome)
<!-- #endregion -->

```R
fit_kzfps %>% dplyr::filter(symbol %in% c('ZNF3', 'ZNF597', 'ZNF20'))
```

Plot for Romain: ZNF577 ZNF649 ZNF765 ZNF93 ZNF141 ZNF560 ZNF490

```R
fit_kzfps %>% dplyr::filter(symbol %in% c('ZNF577', 'ZNF649', 'ZNF765', 'ZNF93', 'ZNF141', 'ZNF560', 'ZNF490'))
```

```R
expr_corrected_long$phase_corrected = factor(new_phases[expr_corrected_long$gate %>% as.integer()], levels = new_phases)
```

```R
expr_corrected_long %>% head()
```

```R
fit_kzfps %>% head()
```

```R
fit_kzfps %>% arrange(acrophase) %>% dplyr::mutate(significance = -log10(padj)) %>% dplyr::select(RepliTiming, mesor, significance) %>% ggplot(aes(y = significance)) + geom_histogram()
```

```R
fit_kzfps %>% colnames()
```

```R
row_anno = fit_kzfps %>% arrange(TEs.pval, acrophase) %>% dplyr::mutate(signif.rhythm = -log10(padj),
                                                             signif.promoter_bind = -log10(promoter.pval),
                                                              signif.exons_bind = -log10(exons.pval),
                                                              signif.TEs_bind = -log10(TEs.pval)) %>% dplyr::select(RepliTiming, rif1, mesor, signif.rhythm, cluster,
                                                                                                                  signif.promoter_bind, signif.exons_bind, signif.TEs_bind)
ha = rowAnnotation(df = row_anno, col = c(cols, 
                                          rif1 =  colorRamp2(c(-0.5, 0, 0.5), rev(c('#f1a340','#f7f7f7','#998ec3'))),
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.promoter_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.exons_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.TEs_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5))))
```

```R
row_anno = fit_kzfps %>% arrange(acrophase) %>% dplyr::mutate(signif.rhythm = -log10(padj),) %>% dplyr::select(RepliTiming, rif1, mesor, signif.rhythm, cluster)
ha = rowAnnotation(df = row_anno, col = c(cols, 
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5))))
```

```R
plot_order = (fit_kzfps %>% arrange(acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases
htmp = Heatmap(to_plot, 
        name = 'z-scored_log2_norm.adj.counts.plus.1', 
        row_title = 'KZFPs',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        heatmap_legend_param = list(
        title = "expr., z-scored log2",
        direction = "horizontal",
        title_position = "topcenter"),
       row_names_gp = grid::gpar(fontsize = 7), right_annotation = ha) %v% NULL

pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap.pdf', height = 28, width = 7.5)

draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")

dev.off()
```

Note: ordering by promoter.pval and TEs.pval does not seem to reveal anything.

```R
fit_kzfps %>% colnames()
```

```R
binding_274_znfs = read.table('../out/temp/kzfp_znfs_strict_sorted_coordinates_znf274.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()

# sorting by clusters, then by coordinates Are there clusters with outliers?
row_anno = fit_kzfps %>% 
    arrange(cluster, chr, OutOfFrameZFAStart, acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj),
                 signif.promoter_bind = -log10(promoter.pval),
                  signif.exons_bind = -log10(exons.pval),
                  signif.TEs_bind = -log10(TEs.pval)) %>% dplyr::select(RepliTiming, rif1, mesor, amplitude, signif.rhythm,
                                                                      signif.promoter_bind, signif.exons_bind, signif.TEs_bind)

ghost_anno = fit_kzfps %>% 
    arrange(cluster, chr, OutOfFrameZFAStart, acrophase) %>% dplyr::select(cluster)

ha = rowAnnotation(df = row_anno, col = c(cols, 
                                        amplitude = colorRamp2(breaks = c(0, 0.1, 0.2, 0.5, 1), colors = viridis::rocket(5)),
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.promoter_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.exons_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.TEs_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5))))

hl = rowAnnotation(foo = anno_block(gp = gpar(), labels = levels(ghost_anno$cluster), which = "rows", labels_rot = 0))

plot_order = (fit_kzfps %>% arrange(cluster, chr, OutOfFrameZFAStart, acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases
htmp1 = Heatmap(to_plot,
        name = 'z-scored_log2_norm.adj.counts.plus.1', 
        row_title = 'KZFPs',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,                
        heatmap_legend_param = list(
        title = "expr., z-scored log2",
        direction = "horizontal",
        title_position = "topcenter"),
       row_names_gp = grid::gpar(fontsize = 7), right_annotation = ha, left_annotation = hl)

rownames(binding_274_znfs) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
htmp2 = Heatmap(binding_274_znfs,
        col = colorRamp2(c(0, 5, 10), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,                
        heatmap_legend_param = list(
        title = "ZNF274 ChIP-seq signal",
        direction = "horizontal",
        title_position = "topcenter"),
       row_names_gp = grid::gpar(fontsize = 7), width = unit(2, "cm"))

pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster.pdf', height = 30, width = 7)
draw(htmp1 + htmp2, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
```

```R
fit_kzfps %>% dplyr::select(age_MA) %>% table()
```

```R
kzfp_tss$idx = 1:nrow(kzfp_tss)
```

```R
to_keep = ((kzfp_tss %>% dplyr::inner_join(kzfp_tss)) %>% arrange(acrophase))$idx
```

```R
kzfp_tss %>% dim()
```

```R
# ZNF274 binding motif present at ZNFs
motif_274 = read.table('../data/martina_processed_files/MB_ZNF274motif_2201_kzfp_jonas_ALL_unique.bed', sep = '\t', header = F, col.names = c("chr", "start", "end", "symbol"))
motif_274 %>% head()
motif_274 %>% dim()
```

We only need the symbol column.

```R
fit_kzfps$motif_274 = fit_kzfps$symbol %in% motif_274$symbol
```

```R
# adding Martina's expression data
ZNF274KO_vs_WT_HEK293T = read.table('../data/martina_processed_files/DE_ZNF274KO_vs_WT_HEK293T.tsv', sep = '\t', header = 1, quote = "")
ZNF274KO_vs_WT_HEK293T %>% head()
ZNF274KO_vs_WT_HEK293T %>% dim()
```

<!-- #raw -->
# adding the mutant KAP1 data # for arianna
kap1_WT_vs_6KR = read.table('../data/SD05_DE_Uninfected_WTvs6KR.csv', sep = ';', header = T)
kap1_WT_vs_6KR$X = NULL
kap1_WT_vs_6KR$X.1 = NULL
colnames(kap1_WT_vs_6KR) = c("ensembl", "transcripts", "symbol", "description", "p.value", "log2FC")
kap1_WT_vs_6KR = kap1_WT_vs_6KR %>% dplyr::mutate(foldChange = sign(log2FC)*(2^abs(log2FC)), p_adj = p.adjust(p.value, method = "BH"))
kap1_WT_vs_6KR %>% head()
<!-- #endraw -->

```R
# we want to add it as a right_annotation with the color for the fold change and the size for the significance
expr_annotation_ZNF274_WT_vs_KO = fit_kzfps %>% 
    arrange(cluster, acrophase) %>%
    dplyr::select(ensembl) %>%
    dplyr::left_join(., ZNF274KO_vs_WT_HEK293T %>% dplyr::select(ensembl, foldChange, p_adj))
expr_annotation_ZNF274_WT_vs_KO = expr_annotation_ZNF274_WT_vs_KO %>% dplyr::mutate(x_center = 1, y_center = nrow(expr_annotation_ZNF274_WT_vs_KO):1)
expr_annotation_ZNF274_WT_vs_KO %>% head()
```

Adding the delta H3K9me3 data in ZNF274_WT_vs_KO:

```R
h3k9me3_delta_znfs_ZNF274_WT_vs_KO = read.table('../out/kzfp_znfs_deltaK9_274KO_293T_sorted_cluster_acrophase.bed', header = F, sep = '\t', na.strings = '.')
colnames(h3k9me3_delta_znfs_ZNF274_WT_vs_KO) = c("znf_chr", "znf_start", "znf_end", "symbol", "znf_CanonicalsZF", "znf_strand",
                                                      "dk9_chr", "dk9_start", "dk9_end", "dk9_foldChange", "dk9_padj")
h3k9me3_delta_znfs_ZNF274_WT_vs_KO %>% head()
```

```R
nrow(h3k9me3_delta_znfs_ZNF274_WT_vs_KO)
nrow(expr_annotation_ZNF274_WT_vs_KO)
```

There's a duplicated ZNF:

```R
duplicated_znf = h3k9me3_delta_znfs_ZNF274_WT_vs_KO$symbol[which(h3k9me3_delta_znfs_ZNF274_WT_vs_KO$symbol %>% duplicated())]
h3k9me3_delta_znfs_ZNF274_WT_vs_KO %>% dplyr::filter(symbol == duplicated_znf)
```

```R
# let's just keep the top deltak9 peak per znf
# careful, slice reorders by the group_by variable!!!

h3k9me3_delta_znfs_ZNF274_WT_vs_KO = h3k9me3_delta_znfs_ZNF274_WT_vs_KO %>% 
    dplyr::mutate(temp_idx = 1:nrow(h3k9me3_delta_znfs_ZNF274_WT_vs_KO)) %>%
    dplyr::group_by(symbol) %>% 
    dplyr::slice_min(dk9_padj) %>% 
    ungroup() %>%
    dplyr::arrange(temp_idx)
```

```R
h3k9me3_delta_znfs_ZNF274_WT_vs_KO %>% head()
```

```R
h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO = h3k9me3_delta_znfs_ZNF274_WT_vs_KO %>% dplyr::mutate(x_center = 1, 
                                                                                          y_center = nrow(h3k9me3_delta_znfs_ZNF274_WT_vs_KO):1,
                                                                                         foldChange = dk9_foldChange,
                                                                                         p_adj = dk9_padj) %>% 
                                                                                        dplyr::select(x_center, foldChange, p_adj) %>%
                                                                                        as.data.frame()
h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO %>% head()
```

<!-- #raw -->
# for Arianna
expr_annotation_KAP1_WT_vs_6KR = fit_kzfps %>% 
    arrange(cluster, acrophase) %>%
    dplyr::select(ensembl) %>%
    dplyr::left_join(., kap1_WT_vs_6KR %>% dplyr::select(ensembl, foldChange, p_adj))
expr_annotation_KAP1_WT_vs_6KR = expr_annotation_KAP1_WT_vs_6KR %>% dplyr::mutate(x_center = 1, y_center = nrow(expr_annotation_KAP1_WT_vs_6KR):1)
expr_annotation_KAP1_WT_vs_6KR %>% head()
<!-- #endraw -->

```R
expr_annotation_ZNF274_WT_vs_KO$foldChange %>% summary()
```

```R
-log10(expr_annotation_ZNF274_WT_vs_KO$p_adj) %>% summary()
```

```R
scale_size <- function(x, x_min, x_max, size, offset) {
    return((x/(x_max-x_min)*size + offset) %>% tidyr::replace_na(0))
    }
```

```R
# maxing out the legend 
thresh_temp = 10
fc_to_plot_274_WT_vs_KO = ifelse(abs(expr_annotation_ZNF274_WT_vs_KO$foldChange>thresh_temp), yes = sign(expr_annotation_ZNF274_WT_vs_KO$foldChange*thresh_temp), no = expr_annotation_ZNF274_WT_vs_KO$foldChange)
#fc_to_plot_KAP1_WT_vs_6KR = ifelse(abs(expr_annotation_KAP1_WT_vs_6KR$foldChange>thresh_temp), yes = sign(expr_annotation_KAP1_WT_vs_6KR$foldChange*thresh_temp), no = expr_annotation_KAP1_WT_vs_6KR$foldChange)
fc_to_plot_h3k9_znfs_274_WT_vs_KO = ifelse(abs(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$foldChange>thresh_temp), yes = sign(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$foldChange*thresh_temp), no = h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$foldChange)
```

```R
thresh_temp_k9_min = 2.5
h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$p_adj = ifelse(abs(fc_to_plot_h3k9_znfs_274_WT_vs_KO) < thresh_temp_k9_min, yes = 1, no = h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$p_adj)
#fc_label = seq(-1*thresh_temp, thresh_temp, length.out = 7)
fc_label = c(-10, -2, 0, 2, 10)
fc_col_fun = colorRamp2(seq(-1*thresh_temp, thresh_temp, length = 7), viridis::turbo(7))

signif_lims_RNAseq = c(-log10(1), -log10(1e-3))
signif_lims_H3K9me3 = c(-log10(1), -log10(1e-8))

signif_size_fun_RNAseq <- function(x) {
    return(unit(scale_size(x, signif_lims_RNAseq[1], signif_lims_RNAseq[2], 2.5, 0.5), "mm"))
    }
signif_label = c(-log10(1), -log10(0.05), -log10(1e-3))

signif_size_fun_H3K9me3 <- function(x) {
    return(unit(scale_size(x, signif_lims_H3K9me3[1], signif_lims_H3K9me3[2], 2.5, 0.5), "mm"))
    }
```

```R
lgd_1 = ComplexHeatmap::Legend(title = "FC\nZNF274KO",
                               title_position = "leftcenter-rot",
                              labels = fc_label,
                            legend_gp = gpar(fill = fc_col_fun(fc_label)),
                             type = "points",
                            pch = 21)

lgd_2 = ComplexHeatmap::Legend(title = "adj. p\nZNF274KO",
                              title_position = "leftcenter-rot",
                              labels = paste0(10^-signif_label),
                            legend_gp = gpar(col = "white",
                                             fill = "black"),
                             type = "points",
                            pch = 21,
                              size = signif_size_fun_RNAseq(signif_label))

fc_anno = HeatmapAnnotation(`diff. expr.` = anno_points(expr_annotation_ZNF274_WT_vs_KO$x_center, 
                                            ylim = c(0.5, 1.5),
                                          width = unit(0.5, "cm"),
                                          pch = 21,
                                         gp = gpar(col = "white",
                                         fill = fc_col_fun(fc_to_plot_274_WT_vs_KO)),
                                        size = signif_size_fun_RNAseq(-log10(expr_annotation_ZNF274_WT_vs_KO$p_adj)),
                                                          axis = F,
                                                          annotation_name_rot = 90),
                            `diff. H3K9me3` = anno_points(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$x_center, 
                                            ylim = c(0.5, 1.5),
                                          width = unit(0.5, "cm"),
                                          pch = 21,
                                         gp = gpar(col = "white",
                                         fill = fc_col_fun(fc_to_plot_h3k9_znfs_274_WT_vs_KO)),
                                        size = signif_size_fun_H3K9me3(-log10(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$p_adj)),
                                                          axis = F,
                                                          annotation_name_rot = 90),
                           which = "row")
```

```R
draw(fc_anno)
```

```R
-log10(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$p_adj) %>% range(na.rm = T)
-log10(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO$p_adj) %>% summary()
```

```R
fc_to_plot_h3k9_znfs_274_WT_vs_KO %>% summary()
```

```R
h3k9me3_delta_znfs_ZNF274_WT_vs_KO %>% dplyr::filter(symbol == "ZNF747")
```

<!-- #raw -->
# for arianna
KAP1_WT_vs_6KR = anno_points(expr_annotation_KAP1_WT_vs_6KR$x_center, 
                                            ylim = c(0.5, 1.5),
                                          width = unit(0.5, "cm"),
                                          pch = 21,
                                         gp = gpar(col = "white",
                                         fill = fc_col_fun(fc_to_plot_KAP1_WT_vs_6KR)),
                                        size = signif_size_fun_RNAseq(-log10(expr_annotation_KAP1_WT_vs_6KR$p_adj)),
                                                          axis = F,
                                                          annotation_name_rot = 90),
<!-- #endraw -->

```R
binding_274_znfs = read.table('../out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
binding_274_znfs_rep2 = read.table('../out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_m01.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
binding_274_znfs_imbeault = read.table('../out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_imbeault.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
binding_274_znfs_begnis = read.table('../out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf274_begnis.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
binding_h3k9me3_znfs = read.table('../out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_h3k9me3.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
# sorting by clusters, then by acrophase. Does ZNF274 binding explain what we see well?
```

## KZFP rhythmicity vs K9 and ZNF274 Heatmaps for the paper

```R
row_anno = fit_kzfps %>% 
    arrange(cluster, acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj)) %>% 
    dplyr::select(RTindex_avg, 
                  rif1, 
                  age_MA) %>%
    dplyr::rename(RTindex = RTindex_avg,
                 RIF1 = rif1,
                 age = age_MA)

# to separate clusters in the heatmap
ghost_anno = fit_kzfps %>% 
    arrange(cluster, acrophase) %>% dplyr::select(cluster)

cols = list(`Repl. timing` = replitiming_colors, 
            MESOR = colorRamp2(breaks = c(-2.5, 2.5, 4, 6, 10), colors = viridis::mako(5)),
             RIF1 =  colorRamp2(c(-0.5, -0.2, 0, 0.2, 0.5), hcl.colors(n = 5, palette = "Cividis", rev = F)),
           RTindex = colorRamp2(breaks = seq(from = 2.5, to = 4.5, length.out = 4), colors = viridis::viridis(4)),
            H3K9me3 = colorRamp2(breaks = seq(from = 70, to = 145, length.out = 10), colors = viridis::magma(10)))

ha = rowAnnotation(df = row_anno, col = c(cols, 
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         CanonicalsZF = colorRamp2(breaks = c(0, 5, 10, 20, 50), colors = viridis::rocket(5)),
                                         age = colorRamp2(breaks = c(20, 50, 90, 105, 150, 320), colors = viridis::rocket(6))),
                          annotation_legend_param = list(title_position = "leftcenter-rot")
)
                   
hl = rowAnnotation(foo = anno_block(gp = gpar(), labels = levels(ghost_anno$cluster), which = "rows", labels_rot = 0))

plot_order = (fit_kzfps %>% arrange(cluster, acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases

htmp_expr_phases = Heatmap(to_plot,
        name = 'expr_phases', 
        row_title = 'KZFPs',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,     
        column_title_side = "top",
        heatmap_legend_param = list(
        title = "expr.\nz-scored log2",
        direction = "vertical",
        title_position = "leftcenter-rot"),
        row_names_gp = grid::gpar(fontsize = 7), 
        right_annotation = ha, 
        left_annotation = hl,
        width = unit(4, "cm"))

rownames(binding_h3k9me3_znfs) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
htmp_h3k9me3_znfs = Heatmap(binding_h3k9me3_znfs, 
        column_title = "H3K9me3",
        column_title_gp = gpar(fontsize = 12, fontfacce = "plain"),
        column_title_side = "bottom",
        column_title_rot = 90,
        #bottom_annotation_height = 0,
        col = colorRamp2(c(6, 15, 40), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_h3k9me3_znfs)),
        heatmap_legend_param = list(
        title = "H3K9me3",
        legend_direction = "vertical",
        title_position = "leftcenter-rot"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(0.5, "cm"))

rownames(binding_274_znfs_rep2) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
htmp_znf274_encode_rep2_znfs = Heatmap(binding_274_znfs_rep2, 
        column_title = "ZNF274",
        column_title_gp = gpar(fontsize = 12, fontface = "plain"),
        column_title_side = "bottom",
        column_title_rot = 90,
        #bottom_annotation_height = 0,
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        row_names_gp = grid::gpar(fontsize = 7, fontface = "italic"), 
        heatmap_legend_param = list(
        title = "ZNF274",
        legend_direction = "vertical",
        title_position = "leftcenter-rot"),
        width = unit(0.5, "cm"))

pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster_acrophase_for_paper_supplementals.pdf', height = 32, width = 10)
draw((htmp_expr_phases + fc_anno + htmp_h3k9me3_znfs + htmp_znf274_encode_rep2_znfs), 
     heatmap_legend_side="right",
     merge_legends = T,
     annotation_legend_side="right", 
     heatmap_legend_list = packLegend(lgd_1, lgd_2))
dev.off()
```

Example heatmaps for the main figures:

Canonical clusters: 

- K9 274 high chr12.1
- K9 intermediate chr19.6
- K9 low chr7.3
- K9 low chr16.2

Non-canonical RIF1-high cluster:
- chr19.8


```R
row_anno = fit_kzfps %>% 
    dplyr::filter(cluster %in% c("chr12.1", "chr7.3", "chr19.6")) %>%
    arrange(cluster, acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj)) %>% 
    dplyr::select(RTindex_avg, 
                  rif1, 
                  age_MA) %>%
    dplyr::rename(RTindex = RTindex_avg,
                 RIF1 = rif1,
                 age = age_MA)

# to separate clusters in the heatmap
ghost_anno = fit_kzfps %>%     
    dplyr::filter(cluster %in% c("chr12.1", "chr7.3", "chr19.6")) %>%

    arrange(cluster, acrophase) %>% dplyr::select(cluster)
ghost_anno$cluster = factor(ghost_anno$cluster)


ha = rowAnnotation(df = row_anno, col = c(cols, 
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         CanonicalsZF = colorRamp2(breaks = c(0, 5, 10, 20, 50), colors = viridis::rocket(5)),
                                         age = colorRamp2(breaks = c(20, 50, 90, 105, 150, 320), colors = viridis::rocket(6))),
                          annotation_legend_param = list(title_position = "leftcenter-rot")
)
                   
hl = rowAnnotation(foo = anno_block(gp = gpar(), labels = levels(ghost_anno$cluster), which = "rows", labels_rot = 0))

plot_order = (fit_kzfps %>% 
            dplyr::filter(cluster %in% c("chr12.1", "chr7.3", "chr19.6")) %>%
            arrange(cluster, acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases

# subindexing the H3K9me3 and ZNF274 heatmaps:
plot_order_no_subselection = (fit_kzfps %>% 
            arrange(cluster, acrophase))$ensembl

idx_subselection = match(plot_order, plot_order_no_subselection)


htmp_expr_phases = Heatmap(to_plot,
        name = 'expr_phases', 
        row_title = 'RIF1-low KZFP clusters',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,     
        column_title_side = "top",
        heatmap_legend_param = list(
        title = "expr.\nz-scored log2",
        direction = "vertical",
        title_position = "leftcenter-rot"),
        row_names_gp = grid::gpar(fontsize = 7), 
        right_annotation = ha, 
        left_annotation = hl,
        width = unit(4, "cm"))

rownames(binding_h3k9me3_znfs) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order_no_subselection, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order_no_subselection, 'stars'])
htmp_h3k9me3_znfs = Heatmap(binding_h3k9me3_znfs[idx_subselection, ],
        column_title = "H3K9me3",
        column_title_gp = gpar(fontsize = 12, fontfacce = "plain"),
        column_title_side = "bottom",
        column_title_rot = 90,
        #bottom_annotation_height = 0,
        col = colorRamp2(c(6, 15, 40), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_h3k9me3_znfs)),
        heatmap_legend_param = list(
        title = "H3K9me3",
        legend_direction = "vertical",
        title_position = "leftcenter-rot"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(0.5, "cm"),
        use_raster = T,
        raster_device = "png")

rownames(binding_274_znfs_rep2) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order_no_subselection, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order_no_subselection, 'stars'])
htmp_znf274_encode_rep2_znfs = Heatmap(binding_274_znfs_rep2[idx_subselection, ], 
        column_title = "ZNF274",
        column_title_gp = gpar(fontsize = 12, fontface = "plain"),
        column_title_side = "bottom",
        column_title_rot = 90,
        #bottom_annotation_height = 0,
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        row_names_gp = grid::gpar(fontsize = 7, fontface = "italic"), 
        heatmap_legend_param = list(
        title = "ZNF274",
        legend_direction = "vertical",
        title_position = "leftcenter-rot"),
        width = unit(0.5, "cm"),
        use_raster = T,
        raster_device = "png")

fc_anno = HeatmapAnnotation(`diff. expr.` = anno_points(expr_annotation_ZNF274_WT_vs_KO[idx_subselection, "x_center"], 
                                            ylim = c(0.5, 1.5),
                                          width = unit(0.5, "cm"),
                                          pch = 21,
                                         gp = gpar(col = "white",
                                         fill = fc_col_fun(fc_to_plot_274_WT_vs_KO[idx_subselection])),
                                        size = signif_size_fun_RNAseq(-log10(expr_annotation_ZNF274_WT_vs_KO[idx_subselection, "p_adj"])),
                                                          axis = F,
                                                          annotation_name_rot = 90),
                            `diff. H3K9me3` = anno_points(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO[idx_subselection, "x_center"], 
                                            ylim = c(0.5, 1.5),
                                          width = unit(0.5, "cm"),
                                          pch = 21,
                                         gp = gpar(col = "white",
                                         fill = fc_col_fun(fc_to_plot_h3k9_znfs_274_WT_vs_KO[idx_subselection])),
                                        size = signif_size_fun_H3K9me3(-log10(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO[idx_subselection, "p_adj"])),
                                                          axis = F,
                                                          annotation_name_rot = 90),
                           which = "row")

pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster_acrophase_for_paper_main.pdf', height = 6, width = 10)
draw((htmp_expr_phases + fc_anno + htmp_h3k9me3_znfs + htmp_znf274_encode_rep2_znfs), 
     heatmap_legend_side="right",
     merge_legends = T,
     annotation_legend_side="right", 
     heatmap_legend_list = packLegend(lgd_1, lgd_2))
dev.off()

svg('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster_acrophase_for_paper_main.svg', height = 6, width = 10)
draw((htmp_expr_phases + fc_anno + htmp_h3k9me3_znfs + htmp_znf274_encode_rep2_znfs), 
     heatmap_legend_side="right",
     merge_legends = T,
     annotation_legend_side="right", 
     heatmap_legend_list = packLegend(lgd_1, lgd_2))
dev.off()
```

```R
row_anno = fit_kzfps %>% 
    dplyr::filter(cluster %in% c("chr19.8")) %>%
    arrange(cluster, acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj)) %>% 
    dplyr::select(RTindex_avg, 
                  rif1, 
                  age_MA) %>%
    dplyr::rename(RTindex = RTindex_avg,
                 RIF1 = rif1,
                 age = age_MA)

# to separate clusters in the heatmap
ghost_anno = fit_kzfps %>%     
    dplyr::filter(cluster %in% c("chr19.8")) %>%

    arrange(cluster, acrophase) %>% dplyr::select(cluster)
ghost_anno$cluster = factor(ghost_anno$cluster)


ha = rowAnnotation(df = row_anno, col = c(cols, 
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         CanonicalsZF = colorRamp2(breaks = c(0, 5, 10, 20, 50), colors = viridis::rocket(5)),
                                         age = colorRamp2(breaks = c(20, 50, 90, 105, 150, 320), colors = viridis::rocket(6))),
                          annotation_legend_param = list(title_position = "leftcenter-rot")
)
                   
hl = rowAnnotation(foo = anno_block(gp = gpar(), labels = levels(ghost_anno$cluster), which = "rows", labels_rot = 0))

plot_order = (fit_kzfps %>% 
            dplyr::filter(cluster %in% c("chr19.8")) %>%
            arrange(cluster, acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases

# subindexing the H3K9me3 and ZNF274 heatmaps:
plot_order_no_subselection = (fit_kzfps %>% 
            arrange(cluster, acrophase))$ensembl

idx_subselection = match(plot_order, plot_order_no_subselection)


htmp_expr_phases = Heatmap(to_plot,
        name = 'expr_phases', 
        row_title = 'RIF1-high KZFP clusters',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,     
        column_title_side = "top",
        heatmap_legend_param = list(
        title = "expr.\nz-scored log2",
        direction = "vertical",
        title_position = "leftcenter-rot"),
        row_names_gp = grid::gpar(fontsize = 7), 
        right_annotation = ha, 
        left_annotation = hl,
        width = unit(4, "cm"))

rownames(binding_h3k9me3_znfs) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order_no_subselection, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order_no_subselection, 'stars'])
htmp_h3k9me3_znfs = Heatmap(binding_h3k9me3_znfs[idx_subselection, ],
        column_title = "H3K9me3",
        column_title_gp = gpar(fontsize = 12, fontfacce = "plain"),
        column_title_side = "bottom",
        column_title_rot = 90,
        #bottom_annotation_height = 0,
        col = colorRamp2(c(6, 15, 40), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_h3k9me3_znfs)),
        heatmap_legend_param = list(
        title = "H3K9me3",
        legend_direction = "vertical",
        title_position = "leftcenter-rot"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(0.5, "cm"),
        use_raster = T,
        raster_device = "png")

rownames(binding_274_znfs_rep2) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order_no_subselection, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order_no_subselection, 'stars'])
htmp_znf274_encode_rep2_znfs = Heatmap(binding_274_znfs_rep2[idx_subselection, ], 
        column_title = "ZNF274",
        column_title_gp = gpar(fontsize = 12, fontface = "plain"),
        column_title_side = "bottom",
        column_title_rot = 90,
        #bottom_annotation_height = 0,
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        row_names_gp = grid::gpar(fontsize = 7, fontface = "italic"), 
        heatmap_legend_param = list(
        title = "ZNF274",
        legend_direction = "vertical",
        title_position = "leftcenter-rot"),
        width = unit(0.5, "cm"),
        use_raster = T,
        raster_device = "png")

fc_anno = HeatmapAnnotation(`diff. expr.` = anno_points(expr_annotation_ZNF274_WT_vs_KO[idx_subselection, "x_center"], 
                                            ylim = c(0.5, 1.5),
                                          width = unit(0.5, "cm"),
                                          pch = 21,
                                         gp = gpar(col = "white",
                                         fill = fc_col_fun(fc_to_plot_274_WT_vs_KO[idx_subselection])),
                                        size = signif_size_fun_RNAseq(-log10(expr_annotation_ZNF274_WT_vs_KO[idx_subselection, "p_adj"])),
                                                          axis = F,
                                                          annotation_name_rot = 90),
                            `diff. H3K9me3` = anno_points(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO[idx_subselection, "x_center"], 
                                            ylim = c(0.5, 1.5),
                                          width = unit(0.5, "cm"),
                                          pch = 21,
                                         gp = gpar(col = "white",
                                         fill = fc_col_fun(fc_to_plot_h3k9_znfs_274_WT_vs_KO[idx_subselection])),
                                        size = signif_size_fun_H3K9me3(-log10(h3k9me3_delta_znfs_annotation_ZNF274_WT_vs_KO[idx_subselection, "p_adj"])),
                                                          axis = F,
                                                          annotation_name_rot = 90),
                           which = "row")

pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster_acrophase_rif1_high_for_paper_main.pdf', height = 3, width = 10)
draw((htmp_expr_phases + fc_anno + htmp_h3k9me3_znfs + htmp_znf274_encode_rep2_znfs), 
     heatmap_legend_side="right",
     merge_legends = T,
     annotation_legend_side="right", 
     heatmap_legend_list = packLegend(lgd_1, lgd_2))
dev.off()

svg('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster_acrophase_rif1_high_for_paper_main.svg', height = 3.5, width = 10)
draw((htmp_expr_phases + fc_anno + htmp_h3k9me3_znfs + htmp_znf274_encode_rep2_znfs), 
     heatmap_legend_side="right",
     merge_legends = T,
     annotation_legend_side="right", 
     heatmap_legend_list = packLegend(lgd_1, lgd_2))
dev.off()
```

<!-- #raw -->
# not used: annotation titles on top of the heatmap
decorate_annotation("RTindex", { 
    grid.text("RTindex", y = unit(1, "npc") + unit(2, "mm"), just = "left", rot = 90) 
})
decorate_annotation("RIF1", { 
    grid.text("RIF1", y = unit(1, "npc") + unit(2, "mm"), just = "left", rot = 90) 
})
decorate_annotation("age", { 
    grid.text("age", y = unit(1, "npc") + unit(2, "mm"), just = "left", rot = 90) 
})
dev.off()

decorate_annotation("diff_expr.", { 
    grid.text("diff. expr", y = unit(1, "npc") + unit(2, "mm"), just = "left", rot = 90) 
})

decorate_annotation("diff_H3K9me3", { 
    grid.text("diff. H3K9me3", y = unit(1, "npc") + unit(2, "mm"), just = "left", rot = 90) 
})
dev.off()
<!-- #endraw -->

```R

x = ghost_anno$cluster %>% as.integer() 
slices_annotation = which(x[-1] != x[-length(x)])

decorate_heatmap_body(heatmap = "expr_phases", 
                      code = {
                          
    x = fit_kzfps %>% 
                          arrange(cluster, acrophase) %>%
                          dplyr::pull(acrophase)
    grid.lines(x = c(x, x), y = c(0, 1), gp = gpar(lwd = 2, lty = 2))
}, #slice = c(1, slices_annotation-1, nrow(ghost_anno))
)
```

```R
slices_annotation
```

```R

```

```R
fit_kzfps %>% 
    arrange(cluster, acrophase) %>% dplyr::pull(acrophase) %>% range()
```

<!-- #raw -->
# legacy: ZNF75D does not correlate with expression rhythmicity
binding_75a_znfs_imbeault = read.table('../out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf75a.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
binding_75d_znfs_imbeault = read.table('../out/temp/kzfp_znfs_strict_sorted_cluster_acrophase_znf75d.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
<!-- #endraw -->

<!-- #raw -->
# legacy: huge exploratory heatmap which is not necessarily useful for the paper
row_anno = fit_kzfps %>% 
    arrange(cluster, acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj),
                 signif.promoter_bind = -log10(promoter.pval),
                  signif.exons_bind = -log10(exons.pval),
                  signif.TEs_bind = -log10(TEs.pval)) %>% dplyr::select(RepliTiming, rif1, mesor, amplitude, signif.rhythm, age_MA,
                                                                      signif.promoter_bind, signif.exons_bind, signif.TEs_bind, 
                                                                        motif_274,
                                                                       CanonicalsZF)

ghost_anno = fit_kzfps %>% 
    arrange(cluster, acrophase) %>% dplyr::select(cluster)


ha = rowAnnotation(df = row_anno, col = c(cols, 
                                        amplitude = colorRamp2(breaks = c(0, 0.1, 0.2, 0.5, 1), colors = viridis::rocket(5)),
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                          age_MA = colorRamp2(breaks = c(0, 30, 45, 100, 300), colors = rev(viridis::viridis(5))),
                                         signif.promoter_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.exons_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.TEs_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         motif_274 = viridis::viridis(2),
                                         CanonicalsZF = colorRamp2(breaks = c(0, 5, 10, 20, 50), colors = viridis::rocket(5))))
                   
hl = rowAnnotation(foo = anno_block(gp = gpar(), labels = levels(ghost_anno$cluster), which = "rows", labels_rot = 0))

plot_order = (fit_kzfps %>% arrange(cluster, acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases
htmp1 = Heatmap(to_plot,
        name = 'z-scored_log2_norm.adj.counts.plus.1', 
        row_title = 'KZFPs',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,                
        heatmap_legend_param = list(
        title = "expr., z-scored log2",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        right_annotation = ha, 
        left_annotation = hl,
        width = unit(4, "cm"))

# to show rownames
htmp2 = Heatmap(binding_274_znfs, 
        column_title = "ZNF274 (K562)",
        col = colorRamp2(c(0, 5, 10), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))


htmp3 = Heatmap(binding_274_znfs_rep2, 
        column_title = "ZNF274 (K562) rep2",
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

rownames(binding_274_znfs_imbeault) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])

htmp4 = Heatmap(binding_274_znfs_imbeault, 
        column_title = "ZNF274 (293T I)",
        col = colorRamp2(c(0, 0.1, 0.2), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

rownames(binding_274_znfs_begnis) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])

htmp5 = Heatmap(binding_274_znfs_begnis, 
        column_title = "ZNF274 (293T B)",
        col = colorRamp2(c(0, 10, 20), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))
<!-- #endraw -->

<!-- #raw -->
# to correct but anyway CTCF does not associate with ZNF274 binding in K562, at least at KZFPs. 
rownames(binding_ctcf_tss) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])

htmp6 = Heatmap(binding_ctcf_tss, 
        column_title = "CTCF (K562)",
        col = colorRamp2(c(0, 40, 70), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_ctcf_tss)),
        heatmap_legend_param = list(
        title = "CTCF",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))
<!-- #endraw -->

<!-- #raw -->
# legacy: huge exploratory heatmap which is not necessarily useful for the paper

rownames(binding_75a_znfs_imbeault) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
htmp7 = Heatmap(binding_75a_znfs_imbeault, 
        column_title = "ZNF75A (293T)",
        col = colorRamp2(c(0, 0.1, 0.2), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF75A",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

rownames(binding_75d_znfs_imbeault) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
htmp8 = Heatmap(binding_75d_znfs_imbeault, 
        column_title = "ZNF75D (293T)",
        col = colorRamp2(c(0, 1.5, 2), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF75D",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

rownames(binding_h3k9me3_znfs) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
htmp9 = Heatmap(binding_h3k9me3_znfs, 
        column_title = "H3K9me3 (K562)",
        col = colorRamp2(c(6, 15, 40), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_h3k9me3_znfs)),
        heatmap_legend_param = list(
        title = "H3K9me3 (ZNFs)",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster_acrophase.pdf', height = 32, width = 10)
#draw((htmp1 + htmp2+ htmp3 + htmp4 + htmp5 + htmp6 + htmp7 + htmp8), heatmap_legend_side="top", annotation_legend_side="right")
#draw((htmp1 + fc_anno + htmp3 + htmp8 + htmp7 + htmp6), heatmap_legend_side="top", annotation_legend_side="right", heatmap_legend_list = packLegend(lgd_1, lgd_2))
draw((htmp1 + fc_anno + htmp3 + htmp8 + htmp7+ htmp9), heatmap_legend_side="top", annotation_legend_side="right", heatmap_legend_list = packLegend(lgd_1, lgd_2))
dev.off()
<!-- #endraw -->

```R
cols$rif1 = colorRamp2(c(-1, 0, 1), rev(c('#f1a340','#f7f7f7','#998ec3')))
```

<!-- #raw -->
# legacy: not used for paper
row_anno = fit_kzfps %>% 
    arrange(cluster, acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj),
                  signif.promoter_bind = -log10(promoter.pval)) %>%
    dplyr::select(RepliTiming, rif1, age_MA, motif_274)

ghost_anno = fit_kzfps %>% 
    arrange(cluster, acrophase) %>% dplyr::select(cluster)


ha = rowAnnotation(df = row_anno, col = c(cols, 
                                        #amplitude = colorRamp2(breaks = c(0, 0.1, 0.2, 0.5, 1), colors = viridis::rocket(5)),
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                          age_MA = colorRamp2(breaks = c(0, 30, 45, 100, 300), colors = rev(viridis::viridis(5))),
                                         signif.promoter_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         #signif.exons_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         #signif.TEs_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         motif_274 = viridis::viridis(2)
                                         #CanonicalsZF = colorRamp2(breaks = c(0, 5, 10, 20, 50), colors = viridis::rocket(5))))
                                          ))
                   
hl = rowAnnotation(foo = anno_block(gp = gpar(), labels = levels(ghost_anno$cluster), which = "rows", labels_rot = 0))

plot_order = (fit_kzfps %>% arrange(cluster, acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases

htmp1 = Heatmap(to_plot,
        name = 'z-scored_log2_norm.adj.counts.plus.1', 
        row_title = 'KZFPs',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,                
        heatmap_legend_param = list(
        title = "expr., z-scored log2",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        right_annotation = ha, 
        left_annotation = hl,
        width = unit(4, "cm"))

rownames(binding_274_znfs_rep2) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])

htmp3 = Heatmap(binding_274_znfs_rep2, 
        column_title = "ZNF274 (K562) rep2",
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))
pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster_acrophase_simplified.pdf', height = 32, width = 10)
#draw((htmp1 + htmp2+ htmp3 + htmp4 + htmp5 + htmp6 + htmp7 + htmp8), heatmap_legend_side="top", annotation_legend_side="right")
draw((htmp1 + fc_anno + htmp3), heatmap_legend_side="top", annotation_legend_side="right", heatmap_legend_list = packLegend(lgd_1, lgd_2))
dev.off()
<!-- #endraw -->

<!-- #raw -->
# legacy: not used for paper
row_anno = fit_kzfps %>% 
    arrange(cluster, acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj),
                  signif.promoter_bind = -log10(promoter.pval)) %>%
    dplyr::select(RepliTiming, rif1, age_MA, motif_274)

ghost_anno = fit_kzfps %>% 
    arrange(cluster, acrophase) %>% dplyr::select(cluster)


ha = rowAnnotation(df = row_anno, col = c(cols, 
                                        #amplitude = colorRamp2(breaks = c(0, 0.1, 0.2, 0.5, 1), colors = viridis::rocket(5)),
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                          age_MA = colorRamp2(breaks = c(0, 30, 45, 100, 300), colors = rev(viridis::viridis(5))),
                                         signif.promoter_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         #signif.exons_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         #signif.TEs_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         motif_274 = viridis::viridis(2)
                                         #CanonicalsZF = colorRamp2(breaks = c(0, 5, 10, 20, 50), colors = viridis::rocket(5))))
                                          ))
                   
hl = rowAnnotation(foo = anno_block(gp = gpar(), labels = levels(ghost_anno$cluster), which = "rows", labels_rot = 0))

plot_order = (fit_kzfps %>% arrange(cluster, acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases

htmp1 = Heatmap(to_plot,
        name = 'z-scored_log2_norm.adj.counts.plus.1', 
        row_title = 'KZFPs',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,                
        heatmap_legend_param = list(
        title = "expr., z-scored log2",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        right_annotation = ha, 
        left_annotation = hl,
        width = unit(4, "cm"))

htmp3 = Heatmap(binding_274_znfs_rep2, 
        column_title = "ZNF274 (K562) rep2",
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

binding_h3k9me3_znfs
htmp9 = Heatmap(binding_h3k9me3_znfs, 
        column_title = "H3K9me3 (K562)",
        col = colorRamp2(c(6, 15, 40), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_h3k9me3_znfs)),
        heatmap_legend_param = list(
        title = "H3K9me3 (ZNFs)",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))
pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_bycluster_acrophase_noSymbol.pdf', height = 12, width = 6)
#draw((htmp1 + htmp2+ htmp3 + htmp4 + htmp5 + htmp6 + htmp7 + htmp8), heatmap_legend_side="top", annotation_legend_side="right")
draw((htmp1 + fc_anno + htmp3 + htmp9), heatmap_legend_side="top", annotation_legend_side="right", heatmap_legend_list = packLegend(lgd_1, lgd_2))
dev.off()
<!-- #endraw -->

<!-- #raw -->
# legacy: not for paper
# all kzfps without p-adj, with the marks but not by cluster, rather by acrophase
temp_idx_swap = match((fit_kzfps %>% dplyr::arrange(acrophase))$ensembl,
                      (fit_kzfps %>% arrange(cluster, acrophase))$ensembl)


row_anno = fit_kzfps %>% 
#    dplyr::filter(!cluster %in% c("chr19.11, chr19.9")) %>%
    arrange(acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj),
                 signif.promoter_bind = -log10(promoter.pval),
                  signif.exons_bind = -log10(exons.pval),
                  signif.TEs_bind = -log10(TEs.pval)) %>% dplyr::select(RepliTiming, mesor, amplitude, signif.rhythm, age_MA,
                                                                      signif.promoter_bind, signif.exons_bind, signif.TEs_bind, 
                                                                        motif_274,
                                                                       CanonicalsZF)


ha = rowAnnotation(df = row_anno, col = c(cols, 
                                        amplitude = colorRamp2(breaks = c(0, 0.1, 0.2, 0.5, 1), colors = viridis::rocket(5)),
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                          age_MA = colorRamp2(breaks = c(0, 30, 45, 100, 300), colors = rev(viridis::viridis(5))),
                                         signif.promoter_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.exons_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.TEs_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         motif_274 = viridis::viridis(2),
                                         CanonicalsZF = colorRamp2(breaks = c(0, 5, 10, 20, 50), colors = viridis::rocket(5))))
                   
plot_order = (fit_kzfps %>% 
                # dplyr::filter(!cluster %in% c("chr19.11, chr19.9")) %>% 
                  arrange(acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases
htmp1 = Heatmap(to_plot,
        name = 'z-scored_log2_norm.adj.counts.plus.1', 
        row_title = 'KZFPs',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,                
        heatmap_legend_param = list(
        title = "expr., z-scored log2",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        right_annotation = ha, 
        width = unit(4, "cm"))

# to show rownames
htmp2 = Heatmap(binding_274_znfs[temp_idx_swap, ], 
        column_title = "ZNF274 (K562)",
        col = colorRamp2(c(0, 5, 10), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))


htmp3 = Heatmap(binding_274_znfs_rep2[temp_idx_swap, ], 
        column_title = "ZNF274 (K562) rep2",
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

htmp4 = Heatmap(binding_274_znfs_imbeault[temp_idx_swap, ], 
        column_title = "ZNF274 (293T I)",
        col = colorRamp2(c(0, 0.1, 0.2), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

htmp5 = Heatmap(binding_274_znfs_begnis[temp_idx_swap, ], 
        column_title = "ZNF274 (293T B)",
        col = colorRamp2(c(0, 10, 20), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))


htmp7 = Heatmap(binding_75a_znfs_imbeault[temp_idx_swap, ], 
        column_title = "ZNF75A (293T)",
        col = colorRamp2(c(0, 0.1, 0.2), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF75A",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

htmp8 = Heatmap(binding_75d_znfs_imbeault[temp_idx_swap, ], 
        column_title = "ZNF75D (293T)",
        col = colorRamp2(c(0, 1.5, 2), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_split = ghost_anno$cluster,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF75D",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

pdf('../out/cell_cycle_figures/cycling_KZFPs_unfiltered_heatmap_acrophase.pdf', height = 32, width = 14)
draw((htmp1 + htmp3 + htmp8 + htmp7), heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
<!-- #endraw -->

<!-- #raw -->
# legacy, not in the paper (simplified by the ZNF274 signal violins)
# all kzfps without p-adj, with the marks but not by cluster, rather by acrophase
temp_idx_swap = match((fit_kzfps %>% dplyr::filter(!cluster %in% c("chr19.11, chr19.9")) %>% dplyr::arrange(acrophase))$ensembl,
                     (fit_kzfps %>% arrange(cluster, acrophase))$ensembl)


row_anno = fit_kzfps %>% 
    dplyr::filter(!cluster %in% c("chr19.11, chr19.9")) %>%
    arrange(acrophase) %>% 
    dplyr::mutate(signif.rhythm = -log10(padj),
                 signif.promoter_bind = -log10(promoter.pval),
                  signif.exons_bind = -log10(exons.pval),
                  signif.TEs_bind = -log10(TEs.pval)) %>% dplyr::select(RepliTiming, mesor, amplitude, signif.rhythm, age_MA,
                                                                      signif.promoter_bind, signif.exons_bind, signif.TEs_bind, 
                                                                        motif_274,
                                                                       CanonicalsZF)


ha = rowAnnotation(df = row_anno, col = c(cols, 
                                        amplitude = colorRamp2(breaks = c(0, 0.1, 0.2, 0.5, 1), colors = viridis::rocket(5)),
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                          age_MA = colorRamp2(breaks = c(0, 30, 45, 100, 300), colors = rev(viridis::viridis(5))),
                                         signif.promoter_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.exons_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         signif.TEs_bind=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         motif_274 = viridis::viridis(2),
                                         CanonicalsZF = colorRamp2(breaks = c(0, 5, 10, 20, 50), colors = viridis::rocket(5))))
                   
plot_order = (fit_kzfps %>% dplyr::filter(!cluster %in% c("chr19.11, chr19.9")) %>% arrange(acrophase))$ensembl
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases
htmp1 = Heatmap(to_plot,
        name = 'z-scored_log2_norm.adj.counts.plus.1', 
        row_title = 'KZFPs',
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,                
        heatmap_legend_param = list(
        title = "expr., z-scored log2",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        right_annotation = ha, 
        width = unit(4, "cm"))

# to show rownames
htmp2 = Heatmap(binding_274_znfs[temp_idx_swap, ], 
        column_title = "ZNF274 (K562)",
        col = colorRamp2(c(0, 5, 10), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))


htmp3 = Heatmap(binding_274_znfs_rep2[temp_idx_swap, ], 
        column_title = "ZNF274 (K562) rep2",
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

htmp4 = Heatmap(binding_274_znfs_imbeault[temp_idx_swap, ], 
        column_title = "ZNF274 (293T I)",
        col = colorRamp2(c(0, 0.1, 0.2), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

htmp5 = Heatmap(binding_274_znfs_begnis[temp_idx_swap, ], 
        column_title = "ZNF274 (293T B)",
        col = colorRamp2(c(0, 10, 20), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))

#rownames(binding_ctcf_tss[temp_idx_swap, ]) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])

pdf('../out/cell_cycle_figures/cycling_KZFPs_noweird19clusters_heatmap_acrophase.pdf', height = 32, width = 14)
draw((htmp1 + htmp2+ htmp3 + htmp4 + htmp5), heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
<!-- #endraw -->

## Enrichment of Delta H3K9me3 and DE at M-to-S1 vs. S2-to-G2 KZFPs:

The universe will be all 337 KZFPs for which we have rhythmicity data.

```R
h3k9me3_delta_znfs_ZNF274_WT_vs_KO %>% head()
```

```R
enrich_stats = fit_kzfps %>% 
    dplyr::arrange(cluster, acrophase) %>% 
    dplyr::select(padj, phase_assigned, cluster) %>%
    cbind(., h3k9me3_delta_znfs_ZNF274_WT_vs_KO)

colnames(enrich_stats)
```

On all RIF1-low KZFPs:

| | KZFPs with K9 loss | KZFPs without K9 loss | total |
|---|:---:|:---:|:---:|
| M-to-S1 KZFPs | x | m-x | m |
| S2-to-G2 KZFPs | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
clusters_rif1high = c("chr6.1", "chr19.8", "chr19.9", "chr19.11")

```

```R
k = nrow(enrich_stats %>% dplyr::filter(dk9_padj < 0.05, dk9_foldChange < 0))
k
```

```R
acr = c('S2','G2','G2/M')
```

```R
x_observed = nrow(enrich_stats %>% dplyr::filter(dk9_padj < 0.05, dk9_foldChange < 0, !phase_assigned %in% acr))
x_observed
```

```R
m = nrow(enrich_stats %>% dplyr::filter(!phase_assigned %in% acr))
n = nrow(enrich_stats %>% dplyr::filter(phase_assigned %in% acr))

m
n
```

```R
x = 0:m # is the variable tested
probs <- dhyper(x, m, n, k, log = FALSE)


# we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
pval_one_sided = sum(probs[x>=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])

pval_two_sided

```

```R
enrich_stats %>% head()
```

```R
enrich_data_for_plot = enrich_stats %>% 
dplyr::mutate(dk9_padj = ifelse(is.na(dk9_padj), yes = 1, no = dk9_padj),
              dk9_foldChange = ifelse(is.na(dk9_foldChange), yes = 0, no = dk9_foldChange)) %>%
dplyr::mutate(H3K9me3 = c("H3K9me3 unaffected", "H3K9me3 loss")[ifelse(test = ((dk9_padj < 0.05) & (dk9_foldChange < 0)),
                                     yes = T,
                                     no = F)+1],
              acr_group = c("M-S1", "S2-G2")[phase_assigned %in% acr+1]) %>%
    dplyr::select(H3K9me3, acr_group) %>% 
    table() %>% 
    as.data.frame()
enrich_data_for_plot
```

```R
total_kzfps = enrich_data_for_plot %>% dplyr::group_by(acr_group) %>% dplyr::summarize(total_KZFPs = sum(Freq))
total_kzfps$y_label = 1.05
total_kzfps$x_label = c(1, 2)
```

```R
enrich_data_for_plot
```

```R
enrich_data_for_plot %>%
    ggplot(aes(x = acr_group, y = Freq, fill = H3K9me3)) + 
    geom_col(position = "fill") + 
    theme_classic() + 
    scale_fill_manual(values = c("black", "grey60")) +
    ylab(expression(paste("Frac. ", italic("KZFPs"), collapse = " "))) + 
    xlab("") + 
    geom_text(data = total_kzfps, aes(x = x_label, y = y_label, label = paste0("n=",as.character(total_KZFPs)), fill = NULL)) +
    theme(axis.line.x = element_blank()) + 
    ggtitle(paste0("p=", format(pval_two_sided, scientific = T, digits = 3), ", Fisher's Exact Test"))

dev.copy(svg, "../out/cell_cycle_figures/deltaK9_enrichment_KZFPs.svg", width = 4, height = 2.5)
dev.off()
```

There is a significant overlap between H3K9me3 loss upon ZNF274 KO and being an S2-to-G2 KZFP.

Note: it works whether focusing on RIF1-low only, and whether focusing on rhythmic KZFPs only.

Fold changes are not significant, but that's because there's not a fold change for every gene, but only when there are two overlapping K9 regions.


Expression: the overlap is only borderline (not) significant, 

But for fold changes we see a clear difference:


| | KZFPs upregulated | KZFPs not upregulated | total |
|---|:---:|:---:|:---:|
| M-to-S1 KZFPs | x | m-x | m |
| S2-to-G2 KZFPs | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
enrich_stats = expr_annotation_ZNF274_WT_vs_KO %>% dplyr::left_join(fit_kzfps %>% dplyr::select(ensembl, symbol, padj, phase_assigned, cluster))
```

```R
enrich_stats %>% head()
```

```R
x_observed = nrow(enrich_stats %>% dplyr::filter(!phase_assigned %in% acr, p_adj < 0.05, foldChange > 0))
x_observed
```

```R
k = nrow(enrich_stats %>% dplyr::filter(p_adj < 0.05, foldChange > 0))
k
```

```R
m = nrow(enrich_stats %>% dplyr::filter(!phase_assigned %in% acr, !is.na(foldChange)))
n = nrow(enrich_stats %>% dplyr::filter(phase_assigned %in% acr, !is.na(foldChange)))

m
n
```

```R
x = 0:m # is the variable tested
probs <- dhyper(x, m, n, k, log = FALSE)


# we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
pval_one_sided = sum(probs[x>=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])

pval_two_sided
```

<!-- #region -->
How can that not be an enrichment??? from 16% expected to 93% ??


What if we take fold changes?
<!-- #endregion -->

Checking whether the fold change is different from zero based on rhythmicity

```R

enrich_data_for_plot = enrich_stats %>% 
    dplyr::mutate(log2FC = sign(foldChange)*log2(abs(foldChange))) %>%
    group_by(c("M-S1", "S2-G2")[phase_assigned %in% acr+1]) %>%
    summarise(P = wilcox.test(log2FC, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", "ns"),
              MaxWidth = max(log2FC)) %>% 
    dplyr::rename(acr_group = 1)
enrich_data_for_plot$y_pos = 11
enrich_data_for_plot
```

Note:s the whole thing works even when only considering rhythmic KZFPs.

```R
# not separating by RIF1-high 
enrich_stats %>% 
    dplyr::mutate(log2FC = sign(foldChange)*log2(abs(foldChange))) %>%
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log2FC)) + 
    geom_violin() + 
    geom_jitter(width = 0.1, height = 0, size = 0.5) +
    geom_boxplot(outlier.shape=NA, width = 0.1) + 
    geom_hline(yintercept = 0, lty = 2) +
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
    geom_text(aes(x = acr_group, 
                   y = y_pos, 
                   label = paste0("Wilcoxon\np=", format(P, scientific = T, digits = 3))), 
               data = enrich_data_for_plot
              )
dev.copy(svg, "../out/cell_cycle_figures/znf274ko_foldchange_mtos1_vs_s2tog2.svg", height = 4, width = 2.5)
dev.off()
```

## RIF1 values per KZFP cluster:

```R
fit_kzfps$age_MA
```

```R
point_size = 0.7
p = fit_kzfps %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(`RIF1 high` = quantile(rif1, probs = 0.75)>0.5,
                    Age = median(age_MA, na.rm = T)) %>%
    ggplot(aes(x = reorder(cluster, rif1, median, decreasing = T), y = rif1, col = `RIF1 high`)) + 
    geom_hline(yintercept = 0.5, lty = 2) +
    geom_boxplot(outlier.size = point_size) + 
    geom_jitter(width = 0.1, height = 0, size = point_size) + 
    #geom_point(aes(y = Age)) + 
    theme_classic() + 
    ylab("RIF1 [log2 FC over input]") + 
    xlab("") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.line.x = element_blank(),
             axis.ticks.x = element_blank(),
         legend.position = "none") +
    scale_color_manual(values = c("grey60", "black"))
p
```

```R
ggsave('../out/cell_cycle_figures/rif1_per_cluster.svg', p, svg, height = 2.5, width = 4)

```

## Differential expression analysis ZNF274KO vs WT, G1 and G2/M

```R
load('../data/znf274ko_cellcycle_facs_sorted/countTable_gene.RData')
countTable %>% head()
sinfo
```

```R
library(edgeR)

```

```R
y = DGEList(counts=countTable, genes=rownames(countTable), group = factor(sinfo$Group))
y = calcNormFactors(y)
y = estimateDisp(y)

```

```R
factor(sinfo$Group)
```

```R
# DE in G1
et = exactTest(y, pair=c("WT_G1", "ZNF274KO_G1"))
de_genes = topTags(et, n = nrow(et$table))
de_genes$table %>% head()

```

```R
fit_kzfps = fit_kzfps %>% dplyr::left_join(de_genes$table %>% 
                                dplyr::select(-logCPM, -PValue) %>% 
                               dplyr::rename(ensembl = genes, 
                                             logFC_ZNF274KO_vs_WT_G1 = logFC, 
                                             padj_DE_ZNF274KO_vs_WT_G1 = FDR))
```

```R
# DE in S/G2:
et = exactTest(y, pair=c("ZNF274KO_S", "WT_S"))
de_genes = topTags(et, n = nrow(et$table))
de_genes$table %>% head()
```

```R
fit_kzfps = fit_kzfps %>% dplyr::left_join(de_genes$table %>% 
                                dplyr::select(-logCPM, -PValue) %>% 
                               dplyr::rename(ensembl = genes, 
                                             logFC_ZNF274KO_vs_WT_S = logFC, 
                                             padj_DE_ZNF274KO_vs_WT_S = FDR))
```

```R
fit_kzfps %>% 
    dplyr::filter(padj_DE_ZNF274KO_vs_WT_G1 < 0.05) %>%
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = -logFC_ZNF274KO_vs_WT_G1)) + 
    geom_violin() + 
    geom_jitter(width = 0.1, height = 0, size = 0.5) +
    geom_boxplot(outlier.shape=NA, width = 0.1) + 
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
    ggpubr::stat_compare_means(label.y = 10.5)
```

Can we show that the log fold change is centered on zero for S2-G2 KZFPs, but not M-S1 KZFPs?

```R
fit_kzfps %>% 
    dplyr::filter(padj_DE_ZNF274KO_vs_WT_G1 < 0.05) %>%
    group_by(c("M-S1", "S2-G2")[phase_assigned %in% acr+1]) %>%
    summarise(P = wilcox.test(-logFC_ZNF274KO_vs_WT_G1, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", "ns"),
              MaxWidth = max(logFC_ZNF274KO_vs_WT_G1))
```

```R
# not significant DE
fit_kzfps %>% 
    group_by(c("M-S1", "S2-G2")[phase_assigned %in% acr+1]) %>%
    summarise(P = wilcox.test(-logFC_ZNF274KO_vs_WT_G1, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", "ns"),
              MaxWidth = max(logFC_ZNF274KO_vs_WT_G1))
```

```R
# DE in S:
fit_kzfps %>% 
    dplyr::filter(padj_DE_ZNF274KO_vs_WT_S < 0.05) %>%
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = -logFC_ZNF274KO_vs_WT_S)) + 
    geom_violin() + 
    geom_jitter(width = 0.1, height = 0, size = 0.5) +
    geom_boxplot(outlier.shape=NA, width = 0.1) + 
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
    ggpubr::stat_compare_means(label.y = 10.5)
```

```R
fit_kzfps %>% 
    dplyr::filter(padj_DE_ZNF274KO_vs_WT_S < 0.05) %>%
    group_by(c("M-S1", "S2-G2")[phase_assigned %in% acr+1]) %>%
    summarise(P = wilcox.test(-logFC_ZNF274KO_vs_WT_S, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", "ns"),
              MaxWidth = max(logFC_ZNF274KO_vs_WT_S))
```

```R
fit_kzfps %>% 
    group_by(c("M-S1", "S2-G2")[phase_assigned %in% acr+1]) %>%
    summarise(P = wilcox.test(-logFC_ZNF274KO_vs_WT_S, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", "ns"),
              MaxWidth = max(logFC_ZNF274KO_vs_WT_S))
```

Indeed, the median fold change for S2-G2 KZFPs is zero, in both groups. Note: there is a borderline non-significant p-value for S2-G2 KZFPs when one does not filter on DE genes. 

Now, we need to show that the fold changes are greater in S than in G1 for the same genes.

```R
df1 = fit_kzfps %>% 
    dplyr::select(symbol, ensembl, phase_assigned, padj, padj_DE_ZNF274KO_vs_WT_G1) %>%
    tidyr::pivot_longer(cols = c("padj_DE_ZNF274KO_vs_WT_G1"),  names_to = "FACS_fraction_padj", values_to = c("padj_DE_ZNF274KO_vs_WT")) %>%
    dplyr::mutate(FACS_fraction = "G1")
df1$FACS_fraction_padj = NULL
df1 %>% head()

```

```R
df2 = fit_kzfps %>% 
    dplyr::select(symbol, ensembl, phase_assigned, padj, logFC_ZNF274KO_vs_WT_G1) %>%
    tidyr::pivot_longer(cols = c("logFC_ZNF274KO_vs_WT_G1"),  names_to = "FACS_fraction_padj", values_to = c("logFC_ZNF274KO_vs_WT")) %>%
    dplyr::mutate(FACS_fraction = "G1")
df2$FACS_fraction_padj = NULL
df2 %>% head()
```

```R
df3 = fit_kzfps %>% 
    dplyr::select(symbol, ensembl, phase_assigned, padj, padj_DE_ZNF274KO_vs_WT_S) %>%
    tidyr::pivot_longer(cols = c("padj_DE_ZNF274KO_vs_WT_S"),  names_to = "FACS_fraction_padj", values_to = c("padj_DE_ZNF274KO_vs_WT")) %>%
    dplyr::mutate(FACS_fraction = "S")
df3$FACS_fraction_padj = NULL
df3 %>% head()
```

```R
df4 = fit_kzfps %>% 
    dplyr::select(symbol, ensembl, phase_assigned, padj, logFC_ZNF274KO_vs_WT_S) %>%
    tidyr::pivot_longer(cols = c("logFC_ZNF274KO_vs_WT_S"),  names_to = "FACS_fraction_padj", values_to = c("logFC_ZNF274KO_vs_WT")) %>%
    dplyr::mutate(FACS_fraction = "S")
df4$FACS_fraction_padj = NULL
df4 %>% head()
```

```R
kzfps_znf274ko_facs = rbind((df1 %>% dplyr::full_join(., df2)), (df3 %>% dplyr::full_join(., df4)))
kzfps_znf274ko_facs = kzfps_znf274ko_facs %>% dplyr::left_join(fit_kzfps %>% dplyr::select(symbol, cluster))
kzfps_znf274ko_facs %>% head()
```

```R
# optional filtering step: DE genes in the S phase comparison:
et = exactTest(y, pair=c("ZNF274KO_S", "WT_S"))
de_genes = topTags(et, n = nrow(et$table))
de_genes$table %>% head()
```

<!-- #raw -->
kzfps_znf274ko_facs %>% dim()
kzfps_znf274ko_facs = kzfps_znf274ko_facs %>% 
    dplyr::filter(ensembl %in% (de_genes$table %>% 
                    dplyr::filter(FDR < 0.05, logFC < 0) %>% dplyr::pull(genes)))
kzfps_znf274ko_facs %>% dim()
<!-- #endraw -->

```R
kzfps_znf274ko_facs %>% 
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = -logFC_ZNF274KO_vs_WT, fill = FACS_fraction)) + 
    geom_violin() + 
    geom_jitter(width = 0.1, height = 0, size = 0.5) +
    geom_boxplot(outlier.shape=NA, width = 0.1) + 
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
    ggpubr::stat_compare_means(label.y = 10.5)
```

```R
kzfps_znf274ko_facs %>% 
    dplyr::filter(!phase_assigned %in% acr) %>%
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = FACS_fraction, y = -logFC_ZNF274KO_vs_WT)) + 
    geom_violin() + 
    geom_jitter(width = 0.1) +
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
        ggpubr::stat_compare_means(label.y = 10.5, paired = T)
```

No difference of fold change between G1 and S fractions for S2-G2 KZFPs -> no bias in cell phase specificity

```R
kzfps_znf274ko_facs %>% 
    dplyr::filter(!phase_assigned %in% acr, !cluster %in% clusters_rif1high) %>%
#    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = FACS_fraction, y = -logFC_ZNF274KO_vs_WT)) + 
    geom_violin() + 
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
        ggpubr::stat_compare_means(label.y = 10.5, paired = T)
```

In RIF1-low clusters, we have a nice, higher DE fold change in S than in G1. M-to-S1-peaking KZFPs react the most in terms of fold change, and particularly during S. Thus, ZNF274-mediated repression is stronger in S.

The result holds if focusing on rhythmic KZFPs

```R
kzfps_znf274ko_facs %>% 
    dplyr::filter(!phase_assigned %in% acr, cluster %in% clusters_rif1high) %>%
#    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = FACS_fraction, y = -logFC_ZNF274KO_vs_WT)) + 
    geom_violin() + 
    geom_jitter(width = 0.1) + 
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
        ggpubr::stat_compare_means(label.y = 10.5, paired = T)
```

What is happening here? Could it be due to selecting too many fold changes? Maybe we should filter on DE in at least one of the contexts, if not in the S fraction. To be tried.

Still holds when focusing on rhythmic KZFPs.

Still holds when limiting to KZFPs with significant changes in the S fraction. So there seems to be no problem in upregulating KZFPs in RIF1-high clusters in S. But is their upregulation the same as KZFPs in RIF-low clusters?


```R
kzfps_znf274ko_facs %>% 
    dplyr::filter(!phase_assigned %in% acr, FACS_fraction == "S") %>%
#    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = cluster %in% clusters_rif1high, y = -logFC_ZNF274KO_vs_WT)) + 
    geom_violin() + 
    geom_jitter(width = 0.1) + 
    geom_boxplot() +
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
        ggpubr::stat_compare_means(label.y = 10.5, paired = F, method = "t.test")
```

```R
kzfps_znf274ko_facs %>% 
    dplyr::filter(FACS_fraction == "S") %>%
#    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = phase_assigned %in% acr, y = -logFC_ZNF274KO_vs_WT)) + 
    geom_violin() + 
    geom_jitter(width = 0.1) + 
    geom_boxplot() +
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
        ggpubr::stat_compare_means(label.y = 10.5, paired = F)
```

```R
kzfps_znf274ko_facs %>% 
    dplyr::filter(phase_assigned %in% acr, cluster %in% clusters_rif1high) %>%
    ggplot(aes(x = FACS_fraction, y = -logFC_ZNF274KO_vs_WT)) + 
    geom_violin() + 
    geom_jitter(width = 0.1) +
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
        ggpubr::stat_compare_means(label.y = 10.5, paired = T)
```

```R
kzfps_znf274ko_facs %>% 
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = logFC_ZNF274KO_vs_WT, fill = FACS_fraction)) + 
    geom_violin() + 
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
    ggpubr::stat_compare_means(label.y = 10.5)
```

```R
fit_kzfps %>% 
    dplyr::filter(padj_DE_ZNF274KO_vs_WT_S < 0.05, !cluster %in% clusters_rif1high) %>%
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = -logFC_ZNF274KO_vs_WT_S)) + 
    geom_violin() + 
    geom_jitter(width = 0.1, height = 0, size = 0.5) +
    geom_boxplot(outlier.shape=NA, width = 0.1) + 
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
    ggpubr::stat_compare_means(label.y = 10.5)
```

```R
fit_kzfps %>% 
    dplyr::filter(padj_DE_ZNF274KO_vs_WT_G1 < 0.05, !cluster %in% clusters_rif1high) %>%
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = -logFC_ZNF274KO_vs_WT_G1)) + 
    geom_violin() + 
    geom_jitter(width = 0.1, height = 0, size = 0.5) +
    geom_boxplot(outlier.shape=NA, width = 0.1) + 
    theme_classic() + 
    ylab("log2(FC)") +
    xlab("") + 
    ggpubr::stat_compare_means(label.y = 10.5)
```

## Other ZNF274 targets

ZNF274 also binds to KRAB-less ZNFs, which tend to cluster with KZFPs. We assign them to KZFP clusters based on distance: 250kb from the middle of the closest KZFP (Jonas's definition of a cluster).

Other targets: check in Martina's paper.

```R
tfs = read.csv('../data/DatabaseExtract_v_1.01.csv', sep = ',', header = T, check.names = T)
tfs %>% head()
```

```R
# keeping only tfs expressed in our data
tfs %>% dim()

tfs_expressed = tfs %>% dplyr::filter(Is.TF. == 'Yes', Ensembl.ID %in% rownames(norm.counts_no_outliers))
tfs_expressed %>% dim()
```

```R
spring_pastels = c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7")
fit_tfs = fit_qualitative %>% dplyr::filter(ensembl %in% tfs_expressed$Ensembl.ID) %>% arrange(acrophase)
```

```R
fit_tfs$is_kzfp = fit_tfs$ensembl %in% fit_kzfps$ensembl
summary(fit_tfs$is_kzfp)
```

```R
tfs_expressed$DBD %>% unique()
```

```R
fit_tfs$is_c2h2 = (tfs_expressed %>% column_to_rownames('Ensembl.ID'))[fit_tfs$ensembl, 'DBD'] %>% grepl('C2H2 ZF', .)
summary(fit_tfs$is_c2h2)
```

```R
fit_tfs$is_tf = fit_tfs$ensembl %in% (tfs_expressed$Ensembl.ID)
summary(fit_tfs$is_tf)
```

```R
fit_tfs$tf_category = ifelse(fit_tfs$is_kzfp, yes = 'KZFP', 
                                    no = ifelse(fit_tfs$is_c2h2, yes = 'C2H2 ZF', 
                                                no = ifelse(fit_tfs$is_tf, yes = 'TF', no = 'not TF')))
table(fit_tfs$tf_category)
```

```R
fit_qualitative %>% colnames()
```

```R
fit_c2h2_nokrab = fit_tfs %>% dplyr::filter(tf_category == "C2H2 ZF") %>% dplyr::arrange(acrophase)
fit_c2h2_nokrab %>% dim()
```

```R
# which KZFPs are cycling ?
c2h2s = read.csv('../data/human_KZFPTable.csv', sep = '\t')
c2h2s %>% head(1)
```

```R
c2h2s %>% dim()
```

```R
# filtering on c2h2
c2h2s = c2h2s %>% dplyr::filter(KRABid == '')
c2h2s %>% dim()
```

```R
colnames(c2h2s)
```

```R
c2h2s$chromosome %>% unique()
```

This is the hg19 assembly, but we need to convert chromosome names, as indicated on this page: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&chromInfoPage=

```R
chr_conversion = read.table("../data/hg19.chromAlias.txt", sep = '\t', header = T)
chr_conversion %>% head()
```

```R
c2h2s = c2h2s %>% left_join(., chr_conversion %>% dplyr::select(ucsc, refseq), by = c("chromosome" = "refseq")) %>% dplyr::rename(chr = ucsc)
c2h2s %>% head()
```

```R
c2h2s = c2h2s %>% dplyr::mutate(entrez = gsub("GeneID:", "", GeneID))
```

```R
c2h2s %>% colnames()
```

```R
# adding the zinc finger array information to the c2h2 fit
fit_c2h2_nokrab = fit_c2h2_nokrab %>% dplyr::left_join(., c2h2s %>% dplyr::select(entrez, chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, Strand, CanonicalsZF), by = c("entrez"))
```

```R
# finding and removing duplicates
fit_c2h2_nokrab %>% dplyr::filter(duplicated(symbol))
```

Some genes, such as ZFHX4, have multiple exons coding for C2H2 ZFs, and thus multiple out of frame zinc finger arrays. We only keep the largest one, i.e. the one with the most canonical zinc fingers.

```R
fit_c2h2_nokrab = fit_c2h2_nokrab %>% dplyr::group_by(ensembl) %>% dplyr::slice_max(n = 1, order_by = CanonicalsZF, with_ties=T) %>% dplyr::ungroup()
fit_c2h2_nokrab %>% dim()
```

This already gets rid of most ZNF doublets.

```R
fit_c2h2_nokrab %>% dplyr::group_by(ensembl) %>% dplyr::filter(n()>1) %>% dplyr::arrange(desc(n()))
```

Many of these multiple out of frame C2H2 ZFs correspond to existing and coding exons. The only way to deal with those, i.e. to choose one of those, is to conservatively take the one with the highest average ZNF274 ChIP-seq signal.


There are also all KZFPs on non-standard chromosomes, which we can remove.

```R
fit_c2h2_nokrab = fit_c2h2_nokrab %>% dplyr::filter((chr.x %in% paste0("chr", c(1:22, "X", "Y")))&(chr.y %in% paste0("chr", c(1:22, "X", "Y"))))
fit_c2h2_nokrab %>% dplyr::group_by(ensembl) %>% dplyr::filter(n()>1)
```

```R
stopifnot(all(fit_c2h2_nokrab$chr.x == fit_c2h2_nokrab$chr.y))
fit_c2h2_nokrab$chr.y = NULL
fit_c2h2_nokrab = fit_c2h2_nokrab %>% dplyr::rename(chr = chr.x)
```

Removing non standard chromosomes removes ~20 duplicated KZNFs.


Checking the remaining Zinc finger array-duplicated ZNFs by hand:
- ZNF64: two exons
- IKZF2: two exons
- CTCF: two exons
- ZBTB16: three exons
- PRDM4: two exons
- BCL6: two exons
- PRDM2: single exon -> merge
- ZNF236: three exons
- ZIC5: two exons
- PRDM15: two exons
- PRDM16: two exons
- IKZF3: two exons
- GFI1: two exons
- ZNF143: two exons
- ZBTB4: two exons
- IKZF1: two exons
- MTF1: two exons
- ZNF292: single exon -> merge

```R
# merging out of frame zinc finger arrays
row_new = (fit_c2h2_nokrab %>% dplyr::filter(symbol=="PRDM2"))[1, ]
```

```R
fit_c2h2_nokrab %>% dplyr::filter(symbol=="PRDM2") %>% dplyr::select(CanonicalsZF) %>% colSums() %>% as.vector()
```

```R
row_new$CanonicalsZF = fit_c2h2_nokrab %>% dplyr::filter(symbol=="PRDM2") %>% dplyr::select(CanonicalsZF) %>% colSums() %>% as.vector()
```

```R
(fit_c2h2_nokrab %>% dplyr::filter(symbol=="PRDM2"))[2, ]$OutOfFrameZFAEnd
```

```R
row_new$OutOfFrameZFAEnd = (fit_c2h2_nokrab %>% dplyr::filter(symbol=="PRDM2"))[2, ]$OutOfFrameZFAEnd
```

```R
fit_c2h2_nokrab = fit_c2h2_nokrab %>% dplyr::filter(!symbol=="PRDM2")
fit_c2h2_nokrab = rbind(fit_c2h2_nokrab, row_new)
```

```R
row_new = (fit_c2h2_nokrab %>% dplyr::filter(symbol=="ZNF292"))[1, ]
```

```R
row_new$CanonicalsZF = fit_c2h2_nokrab %>% dplyr::filter(symbol=="ZNF292") %>% dplyr::select(CanonicalsZF) %>% colSums() %>% as.vector()
```

```R
row_new$OutOfFrameZFAEnd = (fit_c2h2_nokrab %>% dplyr::filter(symbol=="ZNF292"))[2, ]$OutOfFrameZFAEnd
```

```R
fit_c2h2_nokrab = fit_c2h2_nokrab %>% dplyr::filter(!symbol=="ZNF292")
fit_c2h2_nokrab = rbind(fit_c2h2_nokrab, row_new)
```

```R
fit_c2h2_nokrab %>% dplyr::group_by(ensembl) %>% dplyr::filter(n()>1)
```

```R
# exporting znf coordinates for each, ordered by acrophase and DNA start
fit_c2h2_nokrab = fit_c2h2_nokrab %>% arrange(acrophase, OutOfFrameZFAStart)
fit_c2h2_nokrab %>% dim()
fit_c2h2_nokrab %>% dplyr::select(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, symbol, CanonicalsZF, Strand) %>% write.table('../out/c2h2s_nokrab_zfcoordinates_sorted_acrophase_DNAStart.bed', row.names=F, col.names=F, sep='\t', quote = F)
```

```R
znf274_znfs_c2h2_nokrab = read.table('../out/temp/c2h2s_nokrab_zfcoordinates_strict_sorted_acrophase_dnastart_znf274m01.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
```

```R
znf274_znfs_c2h2_nokrab %>% rowSums(na.rm = T) %>% summary()
```

```R
fit_c2h2_nokrab$znf274_znfs = znf274_znfs_c2h2_nokrab %>% rowSums(na.rm = T)
```

```R
fit_c2h2_nokrab %>% dim()
fit_c2h2_nokrab %>% dplyr::group_by(ensembl) %>% dplyr::slice_max(order_by = znf274_znfs, n=1, na_rm=T) %>% ungroup() %>% dim()
```

```R
fit_c2h2_nokrab_top_znf274_znfs = fit_c2h2_nokrab %>% dplyr::group_by(ensembl) %>% dplyr::slice_max(order_by = znf274_znfs, n=1, na_rm=T) %>% ungroup() %>% arrange(acrophase)
```

```R
fit_c2h2_nokrab_top_znf274_znfs %>% dim()
```

```R
fit_c2h2_nokrab_top_znf274_znfs %>% 
    dplyr::select(chr, OutOfFrameZFAStart, OutOfFrameZFAEnd, symbol, CanonicalsZF, Strand) %>% 
    write.table('../out/c2h2s_nokrab_zfcoordinates_top274_sorted_acrophase.bed', row.names=F, col.names=F, sep='\t', quote = F)
```

Median much lower than the mean, therefore there are probably a few 274 binders, and most non-binders.

```R
plot_selection = fit_c2h2_nokrab_top_znf274_znfs %>% dplyr::arrange(acrophase) %>% as.data.frame()
```

```R
replitiming_cols = viridis(length(fit_qualitative$RepliTiming %>% unique()))
replitiming_cols
```

```R
fit_qualitative$RepliTiming %>% unique()
```

```R
names(replitiming_cols) = (fit_qualitative$RepliTiming %>% unique())
```

```R
cols = list(RepliTiming = replitiming_cols, mesor = colorRamp2(breaks = c(-2.5, 0, 2.5, 5, 10), colors = viridis::magma(5)))
```

```R
row_anno = plot_selection %>% dplyr::select(RepliTiming, mesor)
ha = rowAnnotation(df = row_anno, col = cols)
```

```R
group_cols = viridis(length(fit_c2h2_nokrab_top_znf274_znfs$RepliTiming %>% unique()))
group_cols
```

```R
names(group_cols) = (fit_c2h2_nokrab_top_znf274_znfs$RepliTiming %>% unique()) %>% sort()
```

```R
group_cols
```

```R
#cols = list(RepliTiming = group_cols, mesor = colorRamp2(breaks = c(-2.5, 0, 2.5, 5, 10), colors = viridis::magma(5)), cluster = cluster_cols)
cols = list(RepliTiming = group_cols, mesor = colorRamp2(breaks = c(-2.5, 0, 2.5, 5, 10), colors = viridis::magma(5)), rif1 = colorRamp2(c(-1, 0, 1), rev(c('#f1a340','#f7f7f7','#998ec3'))))
```

```R
cols %>% names()
```

```R
row_anno = plot_selection %>% dplyr::mutate(signif.rhythm = -log10(padj)) %>% dplyr::select(RepliTiming, rif1, mesor, amplitude, signif.rhythm)
ha = rowAnnotation(df = row_anno, col = c(cols,
                                          amplitude = colorRamp2(breaks = c(0, 0.1, 0.2, 0.5, 1), colors = viridis::rocket(5)),
                                          signif.rhythm=colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5))))
```

```R
plot_order = plot_selection$ensembl
to_plot = expr_zscore[plot_order, ]
```

```R
to_plot %>% head()
to_plot %>% dim()
```

```R
n_rhythmic = nrow(to_plot)
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases


htmp1 = Heatmap(to_plot,
        name = 'expr, z-scored log2', 
        row_title = paste0(n_rhythmic, ' KRAB-less ZFP genes'),
        heatmap_legend_param = list(
            title = "expr., z-scored log2",
            direction = "horizontal",
            title_position = "topcenter"),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
       row_names_gp = grid::gpar(fontsize = 7), right_annotation = ha)

binding_274_znfs_rep2_c2h2_nokrab = read.table('../out/temp/c2h2s_nokrab_zfcoordinates_top274_sorted_acrophase_znf274m01.tab', sep = '\t', skip = 3, header = F) %>% as.matrix()
rownames(binding_274_znfs_rep2_c2h2_nokrab) = rownames(to_plot)
htmp2 = Heatmap(binding_274_znfs_rep2_c2h2_nokrab, 
        column_title = "ZNF274 (K562) rep2",
        col = colorRamp2(c(0, 70, 110), viridis::magma(3)),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        cluster_row_slices = FALSE,
        column_labels = rep("", times = ncol(binding_274_znfs)),
        heatmap_legend_param = list(
        title = "ZNF274",
        direction = "horizontal",
        title_position = "topcenter"),
        row_names_gp = grid::gpar(fontsize = 7), 
        width = unit(2, "cm"))


pdf('../out/cell_cycle_figures/noKRAB_ZFPs_heatmap.pdf', height = 20, width = 5)
draw((htmp1 + htmp2), heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
```

TODO: find a way to simplify this for the paper, maybe through an enrichment test?


KRAB-less ZNFs with 274 binding are all M-G1 and match those found DE up by martina upon ZNF274 KO in HEK293Ts


We may have to "add" clusters if some ZNFs are close to "noCluster" KZFPs.

We now look for all KRAB-less ZNFs whose middle coordinate is located within 250kb of one of these clusters. All others are noCluster.

First: are there groups of 3 ZNFs or more, KRAB or not, that include noCluster KZFPs? These would be "new clusters" we'd have to create ourselves.

```R
kzfps_jonas %>% colnames()
```

```R
cluster_coords = kzfps_jonas %>% 
    dplyr::select(chr, start, end, cluster) %>% 
    dplyr::mutate(l = end-start, middle = end-as.integer(l/2)) %>% 
    dplyr::group_by(chr, cluster) %>% 
    dplyr::mutate(cluster_start = min(middle), cluster_end = max(middle)) %>%
    dplyr::select(chr, cluster, cluster_start, cluster_end) %>% unique() %>% 
    arrange(chr, cluster, cluster_start)
```

```R
cluster_coords %>% head()
cluster_coords %>% dim()
```

### Protocadherins

```R
plot_selection = fit_qualitative %>% dplyr::filter(grepl("protocadherin", genename)) %>% dplyr::arrange(acrophase)
```

```R
row_anno = plot_selection %>% dplyr::mutate(signif.rhythm = -log10(padj)) %>% dplyr::select(RepliTiming, mesor, signif.rhythm, amplitude)
ha = rowAnnotation(df = row_anno, col = c(cols, 
                                          signif.rhythm = colorRamp2(breaks = c(0, 1, 2, 5, 10), colors = viridis::rocket(5)),
                                         amplitude = colorRamp2(breaks = c(0, 0.1, 0.2, 0.5, 1), colors = viridis::rocket(5))))
```

```R
plot_order = plot_selection$ensembl
to_plot = expr_zscore[plot_order, ]
```

```R
n_rhythmic = nrow(to_plot)
rownames(to_plot) = paste((gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol'], (fit_qualitative %>% column_to_rownames('ensembl'))[plot_order, 'stars'])
colnames(to_plot) = new_phases
htmp = Heatmap(to_plot,
        name = 'expr, z-scored log2', 
        row_title = paste0(n_rhythmic, ' protocadherin genes'),
        heatmap_legend_param = list(
            title = "expr., z-scored log2",
            direction = "horizontal",
            title_position = "topcenter"),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
       row_names_gp = grid::gpar(fontsize = 7), right_annotation = ha) %v% NULL
pdf('../out/cell_cycle_figures/cycling_PCDH_heatmap.pdf', height = 4, width = 5)
draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
```

Most of them seem really M-G1 peaking, with those being exceptions as: (we check their significance in Martina's paper, the Neuronal progenitor and 293T RNA seq):
- PCDHB5: ns, ns
- PCDHGB6: ns, ns
- PCDHGA3: ns, ns
- PCDH9: ns, ns
- PCDHGA10: ns, ns
- PCDHGA12: ns, ns

In other words, PCDH genes that tend to reach max expression in S-G2 are those that do not get signif. upregulated upon ZNF274 KO in WT cells


TODO: compute this as an enrichment test for the paper


Other gene clusters upregulated upon ZNF274 KO in 293T:
- Kallikrein KLK
- pregnancy-specific glycoprotein (PSG)
- keratin (KRT) 
- NOD-like receptor protein (NLRP)

And delta K9 at SNORD116 gene cluster

```R
fit_qualitative %>% dplyr::filter(grepl("kallikrein", genename)) %>% dplyr::arrange(acrophase)
```

Doesnt work for KLK, they still peak in S.

But their RT changes upon ZNF274 KO!

```R
fit_qualitative %>% dplyr::filter(grepl("PSG", symbol)) %>% dplyr::arrange(acrophase)
```

```R
fit_qualitative %>% dplyr::filter(grepl("pregnancy", genename)) %>% dplyr::arrange(acrophase)
```

The PSGs are not found in our dataset

very ery slightly earlier RT, but probably not significant

```R
fit_qualitative %>% dplyr::filter(grepl("keratin", genename)) %>% dplyr::arrange(acrophase)
```

Doesnt work for KRT genes (lots of S-G2 genes)

very very slight earlier RT, maybe significant if doing a paired test on RTindex.

```R
fit_qualitative %>% dplyr::filter(grepl("NLRP", symbol)) %>% dplyr::arrange(acrophase)
```

Doesnt work either, but a very slight earlier RT upon ZNF274KO (NLRP11-NLRP5)

So the other clustered, ZNF274KO upregulated genes in NPCs do not appear as enriched for M-G1 cycling in K562. Note that the upregulation of PCDH was seen in 293Ts, which are not neuronal cells.

```R
fit_tfs %>% colnames()
```

```R
tfs %>% head()
```

```R
plot(fit_kzfps$mesor, fit_kzfps$statistic)
```

```R
cor(fit_kzfps$mesor, fit_kzfps$statistic, method = "spearman")
```

```R
plot(fit_kzfps$amplitude, -log10(fit_kzfps$padj))
```

```R
cor(fit_kzfps$amplitude, fit_kzfps$statistic, method = "spearman")
```

```R
plot(kzfp_tss$acrophase, log(kzfp_tss$znf274_znfs))
```

```R
kzfp_tss %>% dplyr::filter(padj < 0.05) %>% ggplot(aes(x = acrophase, y = log(znf274_znfs))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = acrophase, y = log(znf274_znfs))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = acrophase, y = log(znf274_tss))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = acrophase, y = log(h3k9me3_tss))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = mesor, y = log(h3k9me3_tss+1))) + geom_point()
```

```R
cor(kzfp_tss$mesor, log(kzfp_tss$h3k9me3_tss+1), method = "spearman")
```

```R
kzfp_tss %>% ggplot(aes(x = mesor, y = log(h3k9me3_znfs+1), col = acrophase)) + geom_point()
```

```R
cor(kzfp_tss$mesor, log(kzfp_tss$h3k9me3_znfs+1), method = "spearman")
```

```R
cor(kzfp_tss$mesor, log(kzfp_tss$atac+1), method = "spearman")
```

```R
cor(kzfp_tss$mesor, log(kzfp_tss$h3k4me3+1), method = "pearson")
```

```R
cor(kzfp_tss$mesor, log(kzfp_tss$h3k27ac+1), method = "pearson")
```

```R
kzfp_tss %>% ggplot(aes(x = mesor, y = log(atac+1))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = mesor, y = log(h3k4me1+1))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = mesor, y = log(h3k27ac+1))) + geom_point()
```

```R
kzfp_tss %>% ggplot(aes(x = mesor, y = acrophase)) + geom_point()
```

```R
kzfp_tss %>% dplyr::filter() %>% ggplot(aes(x = acrophase, y = znf274_znfs)) + geom_point() + ylim(c(0, 10))
```

```R
kzfp_tss %>% ggplot(aes(x = phase_assigned, y = log(znf274_znfs))) + geom_boxplot()
```

The trend is here, there is a lower median logZNF274 signal at the ZNFs for S1-G2, than M-G1/S, but that may not be the best cutoff value.


When plotting the signal on the first set of ZNF274 ChIP-seq, we saturated at signal = 10. Let's try to use that function.

```R
ncol_temp = binding_274_znfs %>% ncol()
```

```R
threshold = 7
```

```R
znf274_znfs_sat_mean = binding_274_znfs %>% ifelse(. < threshold, yes = ., no = threshold) %>% apply(MARGIN=1, mean, na.rm = T)
```

```R
kzfp_tss = kzfp_tss %>% arrange(cluster, acrophase)
kzfp_tss$znf274_znfs_sat_mean = znf274_znfs_sat_mean
```

```R
kzfp_tss$znf274_znfs_bins_higher_than_thresh = (binding_274_znfs > threshold) %>% apply(MARGIN = 1, sum, na.rm = T)
```

```R
kzfp_tss %>% ggplot(aes(x = phase_assigned, y = znf274_znfs_bins_higher_than_thresh)) + geom_boxplot() + geom_point()
```

```R
acr = c("G1/S", "S1", "S2", "G2", "G2/M")
kzfp_tss$s2_g2 = ifelse(kzfp_tss$phase_assigned %in% acr, yes = T, no = F)
```

```R
kzfp_tss %>% ggplot(aes(x = s2_g2, y = znf274_znfs_bins_higher_than_thresh)) + geom_boxplot() + geom_point()
```

```R
acr = c("S2", "G2", "G2/M")

```

```R
t.274_znfs = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$znf274_znfs),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$znf274_znfs),
        alternative = "two.sided",
      var.equal = F)
t.274_znfs
```

```R
paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2))
```

```R
kzfp_tss %>% dplyr::filter(phase_assigned %in% acr) %>% pull(znf274_znfs)
```

```R
acr
```

```R
kzfp_tss %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log(znf274_znfs+1))) + 
geom_violin() + 
geom_jitter(width = 0.1, height = 0, size = 0.5) +
geom_boxplot(outlier.shape=NA, width = 0.1) + 
geom_hline(yintercept = median(log(kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr) %>% pull(znf274_znfs) + 1)), lty = 2) +
theme_classic() + 
ylab("log(ZNF274 ChIP-seq)") +
xlab("") +
ggpubr::stat_compare_means(label.y = 10.5)
#ggtitle(paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2)))
dev.copy(svg, "../out/cell_cycle_figures/znf274_znfs_s2-g2_vs_m-g1s.svg", height = 3, width = 2.5)
dev.off()
```

```R
kzfp_tss %>% dplyr::filter(padj < 0.05) %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log(znf274_znfs+1))) + 
geom_violin() + 
geom_jitter(width = 0.1, height = 0, size = 0.5) +
geom_boxplot(outlier.shape=NA, width = 0.1) + 
geom_hline(yintercept = median(log(kzfp_tss %>% dplyr::filter(padj < 0.05, !phase_assigned %in% acr) %>% pull(znf274_znfs) + 1)), lty = 2) +
theme_classic() + 
ylab("log(ZNF274 ChIP-seq)") +
xlab("") +
ggpubr::stat_compare_means(label.y = 10.5)
#ggtitle(paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2)))
dev.copy(svg, "../out/cell_cycle_figures/znf274_rhythmic_KZFPs_znfs_s2-g2_vs_m-g1s.svg", height = 3, width = 2.5)
dev.off()
```

Barely works even for rhythmic KZFPs, without having to remove RIF1 high clusters

```R
# excluding the RIF1 high clusters
kzfp_tss %>% 
    dplyr::filter(padj < 0.05, !cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1")) %>% 
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log(znf274_znfs+1))) + 
geom_violin() + 
geom_jitter(width = 0.1, height = 0, size = 0.5) +
geom_boxplot(outlier.shape=NA, width = 0.1, alpha = 0) + 
geom_hline(yintercept = median(log(kzfp_tss %>% 
                                   dplyr::filter(padj < 0.05, 
                                                 !cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1"),
                                                 !phase_assigned %in% acr) %>% 
                                   pull(znf274_znfs) + 1)), lty = 2) +
theme_classic() + 
ggpubr::stat_compare_means(label.y = 10) +
ylab("log(ZNF274 ChIP-seq)") +
xlab("")
#ggtitle(paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2)))
dev.copy(svg, "../out/cell_cycle_figures/znf274_znfs_no_RIF1_high_rhythmic_MG1_vs_SG2.svg", height = 3, width = 2.5)
dev.off()
```

```R
# excluding the RIF1 high clusters, not limiting to rhythmic KZFPs
kzfp_tss %>% 
    dplyr::filter(!cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1")) %>% 
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log(znf274_znfs+1))) + 
geom_violin() + 
geom_jitter(width = 0.1, height = 0, size = 0.5) +
geom_boxplot(outlier.shape=NA, width = 0.1, alpha = 0) + 
geom_hline(yintercept = median(log(kzfp_tss %>% 
                                   dplyr::filter(
                                                 !cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1"),
                                                 !phase_assigned %in% acr) %>% 
                                   pull(znf274_znfs) + 1)), lty = 2) +
theme_classic() + 
ggpubr::stat_compare_means(label.y = 10) +
ylab("log(ZNF274 ChIP-seq)") +
xlab("")
# ggtitle(paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2)))
dev.copy(svg, "../out/cell_cycle_figures/znf274_znfs_no_RIF1_high_MG1_vs_SG2.svg", height = 3, width = 2.5)
dev.off()
```

Same but for H3K9me3 at ZNFs

```R
# excluding the RIF1 high clusters
kzfp_tss %>% 
    dplyr::filter(padj < 0.05, !cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1")) %>% 
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log(h3k9me3_znfs+1))) + 
geom_violin() + 
geom_jitter(width = 0.1, height = 0, size = 0.5) +
geom_boxplot(outlier.shape=NA, width = 0.1, alpha = 0) + 
geom_hline(yintercept = median(log(kzfp_tss %>% 
                                   dplyr::filter(padj < 0.05, 
                                                 !cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1"),
                                                 !phase_assigned %in% acr) %>% 
                                   pull(h3k9me3_znfs) + 1)), lty = 2) +
theme_classic() + 
ggpubr::stat_compare_means(label.y = 10) +
ylab("log(H3K9me3 ChIP-seq)") +
xlab("")
#ggtitle(paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2)))
dev.copy(svg, "../out/cell_cycle_figures/h3k9me3_znfs_no_RIF1_high_rhythmic_MG1_vs_SG2.svg", height = 3, width = 2.5)
dev.off()
```

```R
# excluding the RIF1 high clusters, not limiting to rhythmic KZFPs
kzfp_tss %>% 
    dplyr::filter(!cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1")) %>% 
    ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log(h3k9me3_znfs+1))) + 
geom_violin() + 
geom_jitter(width = 0.1, height = 0, size = 0.5) +
geom_boxplot(outlier.shape=NA, width = 0.1, alpha = 0) + 
geom_hline(yintercept = median(log(kzfp_tss %>% 
                                   dplyr::filter(
                                                 !cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1"),
                                                 !phase_assigned %in% acr) %>% 
                                   pull(h3k9me3_znfs) + 1)), lty = 2) +
theme_classic() + 
ggpubr::stat_compare_means(label.y = 10) +
ylab("log(H3K9me3 ChIP-seq)") +
xlab("")
# ggtitle(paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2)))
dev.copy(svg, "../out/cell_cycle_figures/h3k9me3_znfs_no_RIF1_high_MG1_vs_SG2.svg", height = 3, width = 2.5)
dev.off()
```

Tendency but not signif...

```R
kzfp_tss %>% colnames()
```

```R
t.274_znfs = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr, !cluster %in% c("chr6.1", "chr19.8", "chr19.9", "chr19.11")))$h3k9me3_tss),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr, !cluster %in% c("chr6.1", "chr19.8", "chr19.9", "chr19.11")))$h3k9me3_tss),
        alternative = "two.sided",
      var.equal = F)
t.274_znfs
```

```R
t.274_znfs = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$h3k9me3_tss),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$h3k9me3_tss),
        alternative = "two.sided",
      var.equal = F)
t.274_znfs
```

```R
t.274_znfs = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$h3k4me3),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$h3k4me3),
        alternative = "two.sided",
      var.equal = F)
t.274_znfs
```

```R
(kzfp_tss %>% dplyr::filter(phase_assigned %in% acr)) %>% nrow()
```

```R
(kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr)) %>% nrow()
```

```R
(kzfp_tss %>% dplyr::filter(acrophase > 3, acrophase < 6)) %>% nrow()
```

```R
(kzfp_tss %>% dplyr::filter(acrophase < 3 | acrophase > 6)) %>% nrow()
```

```R
t.274_znfs_canonicalsZF = t.test(x = log(
    ((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$znf274_znfs)*
    ((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$CanonicalsZF)),
       y = log(
           ((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$znf274_znfs)*
           ((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$CanonicalsZF)),
        alternative = "two.sided",
      var.equal = F)
t.274_znfs_canonicalsZF
```

```R
kzfp_tss %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log(znf274_znfs*CanonicalsZF+1))) + 
geom_boxplot(outlier.shape=NA) + 
geom_point(position = position_jitter(w = 0.1, h = 0), size = 0.5) +
theme_classic() + 
ylab("log(ZNF274 ChIP-seq)") +
xlab("")+
ggtitle(paste0("p-val=", formatC(t.274_znfs_canonicalsZF$p.value, format = "e", digits = 2)))
dev.copy(pdf, "../out/cell_cycle_figures/znf274_znfs_s2-g2_vs_m-g1s_canonicalszf.pdf", height = 2, width = 2)
dev.off()
```

```R
t.h3k9me3_znfs = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$h3k9me3_znfs),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$h3k9me3_znfs),
        alternative = "two.sided",
      var.equal = F)
t.h3k9me3_znfs
```

The significance for the one-sided t-test is p < 0.00054

```R
kzfp_tss %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], 
                                    y = log(h3k9me3_znfs+1))) + 
geom_boxplot(outlier.shape=NA) +
geom_point(position = position_jitter(w = 0.1, h = 0), size = 0.5) +
theme_classic() + 
ylab("log(H3K9me3 ChIP-seq)") +
xlab("")+
ggtitle(paste0("p-val=", formatC(t.h3k9me3_znfs$p.value, format = "e", digits = 2)))
dev.copy(pdf, "../out/cell_cycle_figures/h3k9me3_znfs_s2-g2_vs_m-g1s.pdf", height = 2, width = 2)
dev.off()
```

```R
t.h3k9me3_znfs_canonicalsZF = t.test(x = log(
    ((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$h3k9me3_znfs)*
    ((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$CanonicalsZF)),
       y = log(
           ((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$h3k9me3_znfs)*
           ((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$CanonicalsZF)),
        alternative = "two.sided",
      var.equal = F)
t.h3k9me3_znfs_canonicalsZF
```

```R
kzfp_tss %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], y = log(h3k9me3_znfs*CanonicalsZF+1))) + 
geom_boxplot(outlier.shape=NA) + 
geom_point(position = position_jitter(w = 0.1, h = 0), size = 0.5) +
theme_classic() + 
ylab("log(H3K9me3 ChIP-seq)") +
xlab("")+
ggtitle(paste0("p-val=", formatC(t.h3k9me3_znfs_canonicalsZF$p.value, format = "e", digits = 2)))
dev.copy(pdf, "../out/cell_cycle_figures/h3k9me3_znfs_s2-g2_vs_m-g1s_canonicalszf.pdf", height = 2, width = 2)
dev.off()
```

```R
t.h3k4me3_tss = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$h3k4me3+1),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$h3k4me3+1),
        alternative = "two.sided",
      var.equal = F)
t.h3k4me3_tss
```

```R
kzfp_tss %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], 
                                    y = log(h3k4me1+1))) + 
geom_boxplot(outlier.shape=NA) + 
geom_point(position = position_jitter(w = 0.1, h = 0), size = 0.5) +
theme_classic() + 
ylab("log(H3K4me3 ChIP-seq)") +
xlab("")+
ggtitle(paste0("p-val=", formatC(t.h3k4me3_tss$p.value, format = "e", digits = 2)))
dev.copy(pdf, "../out/cell_cycle_figures/h3k4me3_tss_s2-g2_vs_m-g1s.pdf", height = 2, width = 2)
dev.off()
```

```R
t.h3k27ac_tss = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$h3k27ac+1),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$h3k27ac+1),
        alternative = "two.sided",
      var.equal = F)
t.h3k27ac_tss
```

```R
kzfp_tss %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], 
                                    y = log(h3k27ac+1))) + 
geom_boxplot(outlier.shape=NA) + 
geom_point(position = position_jitter(w = 0.1, h = 0), size = 0.5) +
theme_classic() + 
ylab("log(H3K27ac ChIP-seq)") +
xlab("")+
ggtitle(paste0("p-val=", formatC(t.h3k27ac_tss$p.value, format = "e", digits = 2)))
dev.copy(pdf, "../out/cell_cycle_figures/h3k27ac_tss_s2-g2_vs_m-g1s.pdf", height = 2, width = 2)
dev.off()
```

```R
t.h3k4me1_tss = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$h3k4me1+1),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$h3k4me1+1),
        alternative = "two.sided",
      var.equal = F)
t.h3k4me1_tss
```

```R
kzfp_tss %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], 
                                    y = log(h3k4me1+1))) + 
geom_boxplot(outlier.shape=NA) + 
geom_point(position = position_jitter(w = 0.1, h = 0), size = 0.5) +
theme_classic() + 
ylab("log(H3K4me1 ChIP-seq)") +
xlab("")+
ggtitle(paste0("p-val=", formatC(t.h3k4me1_tss$p.value, format = "e", digits = 2)))
dev.copy(pdf, "../out/cell_cycle_figures/h3k4me1_tss_s2-g2_vs_m-g1s.pdf", height = 2, width = 2)
dev.off()
```

```R
t.atac_tss = t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$atac+1),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$atac+1),
        alternative = "two.sided",
      var.equal = F)
t.atac_tss
```

```R
kzfp_tss %>% ggplot(aes(x = c("M-S1", "S2-G2")[phase_assigned %in% acr+1], 
                                    y = log(atac+1))) + 
geom_boxplot(outlier.shape=NA) + 
geom_point(position = position_jitter(w = 0.1, h = 0), size = 0.5) +
theme_classic() + 
ylab("log(ATAC-seq)") +
xlab("")+
ggtitle(paste0("p-val=", formatC(t.atac_tss$p.value, format = "e", digits = 2)))
dev.copy(pdf, "../out/cell_cycle_figures/ATAC_tss_s2-g2_vs_m-g1s.pdf", height = 2, width = 2)
dev.off()
```

```R
kzfp_tss %>% ggplot(aes(x = s2_g2, y = log(atac+1))) + geom_boxplot() + geom_point()
```

```R
t.test(x = log((kzfp_tss %>% dplyr::filter(!phase_assigned %in% acr))$atac+1),
       y = log((kzfp_tss %>% dplyr::filter(phase_assigned %in% acr))$atac+1),
        alternative = "two.sided",
      var.equal = T)
```

TODO: correct for the number of statistical tests. 274, K9 and K3me1 at promoter should all hold.

Multiply by the number of ZNFs? By the length of the C2H2? 

```R
kzfp_tss$atac
```

```R
plot(kzfp_tss$acrophase, log(kzfp_tss$znf274_znfs))
```

```R
plot(kzfp_tss$acrophase, log(kzfp_tss$h3k9me3_znfs))
```

```R
# adding ZNF274 binding at the ZNF
m2 = matrix(rnorm(50*10), nrow = 50)
ha = rowAnnotation(foo = anno_density(m2, type = "heatmap", width = unit(6, "cm")))
```

```R
fit_qualitative %>% dplyr::filter(symbol == "ZNF850")
```

```R
plot_order = c('ENSG00000196109', 'ENSG00000269067', 'ENSG00000267041', 'ENSG00000001167', 'ENSG00000120837', 'ENSG00000066136', 'ENSG00000198824')
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = (gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol']
rownames(to_plot) = paste0(rownames(to_plot), (fit_qualitative %>% tibble::column_to_rownames('ensembl'))[plot_selection, 'stars'])
colnames(to_plot) = new_phases
htmp = Heatmap(to_plot, 
        heatmap_legend_param = list(
            title = "expr., z-scored log2",
            direction = "horizontal",
            title_position = "topcenter"),
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
       row_names_gp = grid::gpar(fontsize = 7), 
        row_title = "LTR12C-HERV9NC binders")
    

draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")

dev.copy(pdf, '../out/cell_cycle_figures/LTR12_regulators.pdf', height = 3, width = 3)
dev.off()
```

### KZFPs vs other TFs

```R
tfs = read.csv('../data/DatabaseExtract_v_1.01.csv', sep = ',', header = T, check.names = T)
tfs %>% head()
```

```R
# keeping only tfs expressed in our data
tfs %>% dim()

tfs_expressed = tfs %>% dplyr::filter(Is.TF. == 'Yes', Ensembl.ID %in% rownames(norm.counts_no_outliers))
tfs_expressed %>% dim()
```

```R
# how many genes expressed?
fit_qualitative %>% dim()
```

```R

```

```R
spring_pastels = c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7")
plot_selection = fit_qualitative %>% dplyr::filter(padj < 0.05, ensembl %in% tfs_expressed$Ensembl.ID) %>% arrange(acrophase)
```

```R
plot_selection$is_kzfp = plot_selection$ensembl %in% fit_kzfps$ensembl
summary(plot_selection$is_kzfp)
```

```R
tfs_expressed$DBD %>% unique()
```

```R
plot_selection$is_c2h2 = (tfs_expressed %>% column_to_rownames('Ensembl.ID'))[plot_selection$ensembl, 'DBD'] %>% grepl('C2H2 ZF', .)
summary(plot_selection$is_c2h2)
```

```R
plot_selection$is_tf = plot_selection$ensembl %in% (tfs_expressed$Ensembl.ID)
summary(plot_selection$is_tf)
```

```R
plot_selection$tf_category = ifelse(plot_selection$is_kzfp, yes = 'KZFP', 
                                    no = ifelse(plot_selection$is_c2h2, yes = 'C2H2 ZF', 
                                                no = ifelse(plot_selection$is_tf, yes = 'TF', no = 'not TF')))
table(plot_selection$tf_category)
```

```R
# FRACTION OF RHYTHMIC TFs vs detected TFs
plot_selection %>% nrow() # rhythmic TFs
tfs_expressed %>% nrow() # detected TFs
```

```R
fit_kzfps %>% dplyr::filter(padj < 0.05) %>% nrow()
```

```R
fit_kzfps %>% nrow()
plot_selection %>% dplyr::filter(tf_category == "KZFP") %>% nrow()

(plot_selection %>% dplyr::filter(tf_category == "KZFP") %>% nrow()) / (fit_kzfps %>% nrow())
```

```R
(fit_qualitative %>% nrow())
```

```R
(fit_qualitative %>% dplyr::filter(padj < 0.05) %>% nrow()) / (fit_qualitative %>% nrow())
```

```R
(plot_selection %>% nrow()) / (tfs_expressed %>% nrow())
```

```R
1:length(plot_selection$tf_category %>% unique())
```

```R
replitiming_colors = c(viridis(length(plot_selection$RepliTiming %>% unique())))
replitiming_colors
```

```R
plot_selection$RepliTiming %>% unique()
```

```R
names(replitiming_colors) = plot_selection$RepliTiming %>% unique() %>% levels()
```

```R
replitiming_colors
```

```R
plot_selection[which(plot_selection$RepliTiming == "NA"), ]
```

```R
plot_selection$rif1 %>% summary()
```

```R
plot_selection$tf_category = factor(plot_selection$tf_category, levels = c("TF", "C2H2 ZF", "KZFP"))
```

```R
c("a", "b")[c(1, 2)]
```

```R
(viridis::turbo(n = 10))[2]
```

```R
gene_category = viridis::turbo(n = 10)[c(2, 5, 7)]
names(gene_category) = levels(plot_selection$tf_category)
```

```R
cols = list(`Repl. timing` = replitiming_colors, 
            RTindex = colorRamp2(breaks = seq(from = 2.5, to = 4.5, length.out = 4), colors = viridis::viridis(4)),
            MESOR = colorRamp2(breaks = c(-2.5, 2.5, 4, 6, 10), colors = viridis::mako(5)),
             RIF1 =  colorRamp2(c(-0.5, -0.2, 0, 0.2, 0.5), hcl.colors(n = 5, palette = "Cividis", rev = F)),
           `Gene category` = gene_category
)
```

```R
# genes to highlight
genes_to_label = c("ZNF783", "ZNF786", "ZNF519", "ZNF274")               

genes_to_label_idx = match(genes_to_label, (plot_selection %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase) %>% pull(symbol)))

label_annotation = ComplexHeatmap::rowAnnotation(foo = anno_mark(at = genes_to_label_idx, side = "left", padding = 2,
    labels = genes_to_label))
```

```R
# KZFP rich regions in M/G1

idx_kzfp_1 = 1
idx_kzfp_2 = plot_selection %>% 
    dplyr::select(tf_category, phase_assigned) %>% 
    dplyr::mutate(idx = 1:nrow(plot_selection)) %>% 
    dplyr::filter(tf_category == "KZFP", phase_assigned == "lG1") %>% 
    dplyr::pull(idx) %>% 
    max

idx_kzfp_2
```

```R
idx_kzfp_3 = plot_selection %>% 
    dplyr::select(tf_category, phase_assigned) %>% 
    dplyr::mutate(idx = 1:nrow(plot_selection)) %>% 
    dplyr::filter(tf_category == "KZFP", phase_assigned == "M") %>% 
    dplyr::pull(idx) %>% 
    min
idx_kzfp_4 = nrow(plot_selection)

idx_kzfp_3
```

```R
c("a", "b", "a")[c("a", "b")]
```

```R
split(1:4, c("a", "b", "a", "b"))[c("a")]
```

```R
# indicating KZFP-rich regions
panel_fun = function(index, nm) {
    #pushViewport(viewport(xscale = c(-0.1, 0.1), yscale = c(-1, 1)))
    #grid.rect()
    #grid.xaxis(gp = gpar(fontsize = 8))
    t = ifelse(index < idx_kzfp_3, yes = "G1 KZFPs", no = "M KZFPs")
    grid.text(label = t, x = 2.5, y = 0.5, gp = gpar(fontsize = 9))
    #popViewport()
}

align_to = list(KZFP_g1 = idx_kzfp_1:idx_kzfp_2,
                KZFP_m = idx_kzfp_3:idx_kzfp_4)
```

```R
anno = anno_link(align_to = align_to, which = "row", panel_fun = panel_fun, 
    size = unit(0, "cm"), gap = unit(0, "cm"), width = unit(0, "cm"), side = "left")
label_annotation = ComplexHeatmap::rowAnnotation(foo = anno_mark(at = genes_to_label_idx, side = "left", padding = 2,
    labels = genes_to_label, labels_gp = gpar(fontsize = 9)), bar = anno)
```

```R
row_anno = plot_selection %>% dplyr::select(RTindex_avg, rif1, mesor, tf_category) %>% 
    dplyr::rename(RTindex = RTindex_avg, RIF1 = rif1, MESOR = mesor, `Gene category` = tf_category)
ha = rowAnnotation(df = row_anno, col = cols)
plot_order = plot_selection$ensembl
```

```R
expr_color = colorRamp2(breaks = c(-2, -1, 0, 1, 2), colors = hcl.colors(5, palette = "RdBu", rev = T))
```

```R
# compressed, no TF names
# all cycling TFs
plot_order = plot_selection$ensembl
n_rhythmic = length(plot_order)
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = (gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol']
colnames(to_plot) = new_phases
htmp = Heatmap(to_plot, 
        name = 'expr., z-scored log2', 
        show_row_names = F,
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        row_title = paste0(n_rhythmic, ' rhythmic TFs (padj < 0.05)'),
        heatmap_legend_param = list(
            title = "expr., z-scored log2",
            direction = "horizontal",
            title_position = "topcenter"),
        right_annotation = ha,
        left_annotation = label_annotation,
       row_names_gp = grid::gpar(fontsize = 7),
              col = expr_color) %v% NULL

pdf('../out/cell_cycle_figures/cycling_TFs_no_names_heatmap.pdf', height = 5, width = 5)
draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()

svg('../out/cell_cycle_figures/cycling_TFs_no_names_heatmap.svg', height = 5, width = 5)
draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()
```

```R
p = plot_selection %>% dplyr::rename(`Gene category` = tf_category) %>% ggplot(aes(phase_assigned, fill = `Gene category`)) + geom_bar() + ggplot2::scale_fill_manual(values = cols$`Gene category`) + theme_classic() + ylab("Number of TFs") + xlab("acrophase")


ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_peak_phases.pdf", p, pdf, height = 2, width = 4)
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_peak_phases.svg", p, svg, height = 2, width = 4)


```

Computing the enrichment using Fisher's exact test: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/

The two gene categories are (among all significantly rhythmic genes):
- KZFP
- TF excluding KZFPs

The two acrophase categories are:
- phi in [S1-G2]
- phi in [M-G1/S]



```R
k = plot_selection %>% dplyr::filter(tf_category == "KZFP") %>% nrow() # number of KZFPs

acr = c( "S2", "S/G2", "G2")

m = plot_selection %>% dplyr::filter(phase_assigned %in% acr) %>% nrow() # number of TFs in S-G2

n = plot_selection %>% dplyr::filter(!phase_assigned %in% acr) %>% nrow() # number of TFs in M-G1/S

x = 0:m # is the variable tested

x_observed = plot_selection %>% dplyr::filter(tf_category == "KZFP", phase_assigned %in% acr) %>% nrow() # number of KZFPs in S-G2,

probs <- dhyper(x, m, n, k, log = FALSE)

# we compute the probability of observing a more extreme depletion, therefore using a one sided test. 

pval_kzfp_vs_TFs_one_sided = sum(probs[x<=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_kzfp_vs_TFs_two_sided = sum(probs[probs <= pval_kzfp_vs_TFs_one_sided])
pval_kzfp_vs_TFs_two_sided
x_observed
k
```

Same test but with KRAB-less C2H2 ZFs vs other TFs

```R
k = plot_selection %>% dplyr::filter(tf_category == "C2H2 ZF") %>% nrow() # number of KRAB-less C2H2


acr = c("S2", "S/G2", "G2")

m = plot_selection %>% dplyr::filter(phase_assigned %in% acr) %>% nrow() # number of TFs in S-G2

n = plot_selection %>% dplyr::filter(!phase_assigned %in% acr) %>% nrow() # number of TFs in M-G1/S

x = 0:m # is the variable tested

x_observed = plot_selection %>% dplyr::filter(tf_category == "C2H2 ZF", phase_assigned %in% acr) %>% nrow() # number of KRAB-less C2H2 in S-G2,

probs <- dhyper(x, m, n, k, log = FALSE)

# we compute the probability of observing a more extreme depletion, therefore using a one sided test. 

pval_c2h2_vs_TFs_one_sided = 1-sum(probs[x<=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_c2h2_vs_TFs_two_sided = sum(probs[probs <= pval_c2h2_vs_TFs_one_sided])
pval_c2h2_vs_TFs_two_sided
```

```R
# same for rhythmic TFs amongst rhythmic genes

k = plot_selection %>% nrow() # number of rhythmic TFs


acr = c( "S2", "S/G2", "G2")

m = fit_qualitative %>% dplyr::filter(padj < 0.05, phase_assigned %in% acr) %>% nrow() # number of genes in S-G2

n = fit_qualitative %>% dplyr::filter(padj < 0.05, !phase_assigned %in% acr) %>% nrow() # number of genes in M-G1/S

x = 0:m # is the variable tested

x_observed = plot_selection %>% dplyr::filter(phase_assigned %in% acr) %>% nrow() # number of TFs in S-G2,

probs <- dhyper(x, m, n, k, log = FALSE)

# we compute the probability of observing a more extreme depletion, therefore using a one sided test. 

pval_TFs_vs_genes_one_sided = 1-sum(probs[x<=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_TFs_vs_genes_two_sided = sum(probs[probs <= pval_TFs_vs_genes_one_sided])
pval_TFs_vs_genes_two_sided

x_observed
k
m
```

```R
n_tf_left = plot_selection %>% dplyr::filter(phase_assigned %in% acr, tf_category == "TF") %>% nrow()
n_kzfp_left = plot_selection %>% dplyr::filter(phase_assigned %in% acr, tf_category == "KZFP") %>% nrow()
n_C2H2_left = plot_selection %>% dplyr::filter(phase_assigned %in% acr, tf_category == "C2H2 ZF") %>% nrow()

n_tf_right = plot_selection %>% dplyr::filter(!phase_assigned %in% acr, tf_category == "TF") %>% nrow()
n_kzfp_right = plot_selection %>% dplyr::filter(!phase_assigned %in% acr, tf_category == "KZFP") %>% nrow()
n_C2H2_right = plot_selection %>% dplyr::filter(!phase_assigned %in% acr, tf_category == "C2H2 ZF") %>% nrow()
```

```R
p = plot_selection %>% dplyr::rename(`Gene category`=tf_category) %>% ggplot(aes(ifelse(phase_assigned %in% acr, yes = "S2-G2", no = "M-S1"), fill = `Gene category`)) + 
    geom_bar() + 
    ggplot2::scale_fill_manual(values = cols$`Gene category`) + 
    theme_classic() + 
    ylab("Number of rhythmic TFs (p-adj < 0.05)") + 
    xlab("acrophase") + 
    annotate(geom="text", x=2, y=n_tf_right+n_kzfp_right/2, label=paste0("KZFP p-val=", formatC(pval_kzfp_vs_TFs_two_sided, format = 'e', digits=2))) + 
    annotate(geom="text", x=2, y=n_tf_right+n_kzfp_right + n_C2H2_right/2, label=paste0("C2H2-ZF p-val=", formatC(pval_c2h2_vs_TFs_two_sided, format = 'e', digits=2)))
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_enrichment.pdf", p, pdf, height = 5, width = 5)
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_enrichment.svg", p, svg, height = 4, width = 5)

```

KZFPs are significantly depleted in G1/S-G2 TFs. 

Now what among C2H2 ZFs?

```R
k = plot_selection %>% dplyr::filter(tf_category == "KZFP") %>% nrow() # number of KZFPs


m = plot_selection %>% dplyr::filter(tf_category %in% c("KZFP", "C2H2 ZF"), phase_assigned %in% acr) %>% nrow() # number of C2H2 in S-G2

n = plot_selection %>% dplyr::filter(tf_category %in% c("KZFP", "C2H2 ZF"), !phase_assigned %in% acr) %>% nrow() # number of C2H2 in M-G1/S

x = 0:m # is the variable tested

x_observed = plot_selection %>% dplyr::filter(tf_category == "KZFP", phase_assigned %in% acr) %>% nrow() # number of KZFPs in S-G2, 

probs <- dhyper(x, m, n, k, log = FALSE)

pval = sum(probs[x<=x_observed])
pval
```

KZFPs are significantly depleted from G1/S to G2, even among C2H2 zinc fingers.


C2H2 ZFs are not depleted in G1/S-G2 compared to other TFs.

```R
p = plot_selection %>% dplyr::rename(`Gene category`=tf_category) %>% ggplot(aes(RepliTiming, fill = `Gene category`)) + geom_bar() + ggplot2::scale_fill_manual(values = cols$`Gene category`) + theme_classic() + ylab("Number of TFs") + xlab("RepliTiming") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
```

```R
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_RepliTiming.pdf", p, pdf, height = 2, width = 3.7)
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_RepliTiming.svg", p, svg, height = 2, width = 3.7)
```

```R
# with continuous RT index values
p = plot_selection %>% dplyr::rename(`Gene category`=tf_category, RTindex = RTindex_avg) %>% 
    ggplot(aes(x = `Gene category`, y = RTindex)) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, outlier.size = 0.5) + 
    ggpubr::stat_compare_means(comparisons = list(c(1, 2),
                                                 c(1,3),
                                                 c(2,3)))+
    theme_classic()+ 
    xlab("") +
    theme(axis.line.x = element_blank(),
         axis.ticks.x = element_blank())
p
```

```R
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_RepliTiming_continuous.svg", p, svg, height = 4, width = 3.7)
```

```R
plot_selection %>% colnames()
```

```R
# with continuous RT index values
p = plot_selection %>% dplyr::rename(`Gene category`=tf_category, RIF1 = rif1) %>% 
    ggplot(aes(x = `Gene category`, y = h3k9me3)) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, outlier.size = 0.5) + 
    ggpubr::stat_compare_means(comparisons = list(c(1, 2),
                                                 c(1,3),
                                                 c(2,3)))+
    theme_classic()+ 
    xlab("") +
    theme(axis.line.x = element_blank(),
         axis.ticks.x = element_blank())
p
```

```R
# S2-G2 TFs and ZNFs vs M-S1 TFs and ZNFs, w.r.t RTindex and RIF1
# ranking replication timing per cluster
p = plot_selection %>% 
    dplyr::filter(tf_category == "C2H2 ZF") %>%
    ggplot(aes(x = phase_assigned %in% acr,y = RTindex_avg)) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(RTindex_avg) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("RT index") + 
    xlab("") + 
    theme_classic()
p

```

```R
# S2-G2 TFs and ZNFs vs M-S1 TFs and ZNFs, w.r.t RTindex and RIF1
# ranking replication timing per cluster
p = plot_selection %>% 
    dplyr::filter(tf_category != "KZFP") %>%
    ggplot(aes(x = phase_assigned %in% acr,y = RTindex_avg)) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(RTindex_avg) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("RT index") + 
    xlab("") + 
    theme_classic()
p

```

```R
# S2-G2 TFs and ZNFs vs M-S1 TFs and ZNFs, w.r.t RTindex and RIF1
# ranking replication timing per cluster
p = plot_selection %>% 
    # dplyr::filter(tf_category == "KZFP") %>%
    ggplot(aes(x = phase_assigned %in% acr,y = rif1)) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
ggpubr::stat_compare_means(label.y = 1.15) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("RT index") + 
    xlab("") + 
    theme_classic() + 
    facet_wrap(facets = "tf_category")
p
```

```R
ggsave("../out/cell_cycle_figures/rif1_across_KZFPs_ZFPs_and_TFs.svg", p, svg, width = 6, height = 3)
```

```R
# mesor vs peak expression

p = plot_selection %>% 
    dplyr::filter(tf_category == "KZFP") %>%
    ggplot(aes(x = phase_assigned %in% acr,y = mesor)) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(mesor) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 7) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("MESOR") + 
    xlab("") + 
    theme_classic()
p

```

```R
ggsave("../out/cell_cycle_figures/KFZPs_rhythmic_MESOR_vs_acrophase.svg", p, svg, height = 3, width = 2.5)
```

```R
# rif1 vs RT

p = plot_selection %>% 
    ggplot(aes(x = rif1,y = RTindex_avg)) + 
    geom_point(size = 0.5) +
    ggpubr::stat_cor(method = "spearman", label.y = 5.5) +
    theme_classic() + 
    ylab("RTindex")  +
    xlab("RIF1") + 
    facet_wrap(facets = "tf_category")
p

```

```R
ggsave("../out/cell_cycle_figures/rif1_vs_RT_across_rhythmic_KZFPs_ZFPs_and_TFs.svg", p, svg, width = 6, height = 2.5)
```

```R
# age vs acrophase, all KZFPs, via a permutation test for the difference between mean ages.

# too many ties for a wilcoxon

age_means= fit_kzfps %>% 
    dplyr::filter(!is.na(age_MA)) %>%
    dplyr::mutate(in_S2_to_M = phase_assigned %in% acr) %>%
    dplyr::group_by(in_S2_to_M) %>% 
    dplyr::summarize(age_mean = mean(age_MA), n_kzfps = n())
age_means
```

```R
age_means_diff_observed = diff(age_means$age_mean)
age_means_diff_observed
```

Now, running permutations: we split the KZFPs in groups of the same size as the one observed and compute the difference in mean age. 

```R
to_permute = fit_kzfps %>% 
    dplyr::filter(!is.na(age_MA)) %>% dplyr::pull(age_MA)
stopifnot(length(to_permute) == (age_means$n_kzfps %>% sum()))
```

```R
to_permute %>% head()
```

```R
set.seed(5)
k_draws = 1e6
diff_means_permuted = rep(NA, times = k_draws)
for (i in c(1:k_draws)) {
    s = sample(1:length(to_permute), size = min(age_means$n_kzfps), replace = F) 
    diff_means_permuted[i] = (to_permute[s] %>% mean()) - (to_permute[-s] %>% mean())
    }
```

```R
diff_means_permuted %>% hist()
```

```R
# two sided pval
permute_pval = sum(abs(diff_means_permuted) > age_means_diff_observed)/k_draws
permute_pval
```

```R
library(grid)
# Create a text
grob <- grobTree(textGrob(paste0(as.character(k_draws), " permutations, p = ", as.character(permute_pval)), x=0.1,  y=0.95, hjust=0, gp=gpar(fontsize=9)))
# Plot
p = fit_kzfps %>% dplyr::filter(!is.na(age_MA)) %>% 
    ggplot(aes(x = phase_assigned %in% acr, y = age_MA)) + 
    geom_violin() + 
    geom_point(position = position_jitter(width = 0.1, height = NULL), size = 0.6) + 
    stat_compare_means(label.x = 0.9, label.y = 330, step.increase = 1) + # wilcoxon is not correct when there are many ties
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("age (MY)") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    coord_cartesian(ylim = c(0, 350),
                      clip = 'off') +
    annotation_custom(grob)
p
```

```R
ggsave('../out/cell_cycle_figures/KZFPs_acrophase_vs_age.svg', p, svg, width = 3, height = 3)
```

Avoiding using the age in a Wilcoxon test since it is rank-based and won't like the presence of too many equalities.

We can use a permutation test, to see how often we find a difference in mean more extreme than the one we observed.

```R
age_means= fit_kzfps %>% 
    dplyr::filter(!is.na(age_MA), padj < 0.05) %>%
    dplyr::mutate(in_S2_to_M = phase_assigned %in% acr) %>%
    dplyr::group_by(in_S2_to_M) %>% 
    dplyr::summarize(age_mean = mean(age_MA), n_kzfps = n())
age_means
```

```R
age_means_diff_observed = diff(age_means$age_mean)
age_means_diff_observed
```

Now, running permutations: we split the KZFPs in groups of the same size as the one observed and compute the difference in mean age. 

```R
to_permute = fit_kzfps %>% 
    dplyr::filter(!is.na(age_MA), padj < 0.05) %>% dplyr::pull(age_MA)
stopifnot(length(to_permute) == (age_means$n_kzfps %>% sum()))
```

```R
to_permute %>% head()
```

```R
set.seed(5)
k_draws = 1e6
diff_means_permuted = rep(NA, times = k_draws)
for (i in c(1:k_draws)) {
    s = sample(1:length(to_permute), size = min(age_means$n_kzfps), replace = F) 
    diff_means_permuted[i] = (to_permute[s] %>% mean()) - (to_permute[-s] %>% mean())
    }
```

```R
diff_means_permuted %>% hist()
```

```R
permute_pval = sum(abs(diff_means_permuted) > age_means_diff_observed)/k_draws
permute_pval
```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval), padj < 0.05) %>% dim()
```

```R
library(grid)
# Create a text
grob <- grobTree(textGrob(paste0(as.character(k_draws), " permutations, p = ", as.character(permute_pval)), x=0.1,  y=0.95, hjust=0, gp=gpar(fontsize=9)))
# Plot
p = fit_kzfps %>% dplyr::filter(!is.na(age_MA), padj < 0.05) %>% 
    ggplot(aes(x = phase_assigned %in% acr, y = age_MA)) + 
    geom_violin() + 
    geom_point(position = position_jitter(width = 0.1, height = NULL), size = 0.6) + 
    #stat_compare_means(label.x = 0.9, label.y = 330, step.increase = 1) + # wilcoxon is not correct when there are many ties
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("age (MY)") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    coord_cartesian(ylim = c(0, 350),
                      clip = 'off') +
    annotation_custom(grob)
p
```

```R
ggsave('../out/cell_cycle_figures/KZFPs_acrophase_vs_age_signif.svg', p, svg, width = 3, height = 3)
```

```R
clusters_rif1high = c("chr6.1", "chr19.8", "chr19.11", "chr19.9")
```

```R
p = fit_kzfps %>% 
    dplyr::mutate(in_rif1_high = ifelse(cluster %in% clusters_rif1high, yes = "RIF1-high", no = "RIF1-low")) %>%
    ggplot(aes(x = phase_assigned %in% acr, y = log(h3k9me3_znfs+1))) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, outlier.size = 0.5) + 
    geom_jitter(width = 0.1, size = 0.5) +
    stat_compare_means(label.x = 0.9, label.y = 9, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("H3K9me3 ZNFs") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    coord_cartesian(ylim = c(3.9, 9.5),
                      clip = 'off') + 
    facet_wrap("in_rif1_high")
p
```

```R
ggsave('../out/cell_cycle_figures/KZFPs_h3k9me3_znf_vs_acrophase.svg', p, svg, width = 4, height = 3)
```

```R
clusters_rif1high = c("chr6.1", "chr19.8", "chr19.9", "chr19.11")
```

```R
p = fit_kzfps %>% dplyr::mutate(in_rif1_high = ifelse(cluster %in% clusters_rif1high, yes = "RIF1-high", no = "RIF1-low")) %>%
    ggplot(aes(x = phase_assigned %in% acr, y = log(znf274_znfs+1))) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, outlier.size = 0.5) + 
    geom_jitter(width = 0.1, size = 0.5) +
    stat_compare_means(label.x = 0.9, label.y = 10.5, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("ZNF274 ChIP-seq") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    coord_cartesian(ylim = c(2, 11),
                      clip = 'off') + 
    facet_wrap("in_rif1_high")

p
```

```R
ggsave('../out/cell_cycle_figures/KZFPs_znf274_znf_vs_acrophase_rif1split.svg', p, svg, width = 4, height = 3)
```

```R
p = fit_kzfps %>%
    ggplot(aes(x = phase_assigned %in% acr, y = log(znf274_znfs+1))) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, outlier.size = 0.5) + 
    geom_jitter(width = 0.1, size = 0.5) +
    stat_compare_means(label.x = 0.9, label.y = 10.5, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("ZNF274 ChIP-seq") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    coord_cartesian(ylim = c(2, 11),
                      clip = 'off')
p
```

```R
ggsave('../out/cell_cycle_figures/KZFPs_znf274_znf_vs_acrophase_rif1split.svg', p, svg, width = 2, height = 3)
```

```R
p = fit_kzfps %>% dplyr::filter(!is.na(TEs.pval), !cluster %in% clusters_rif1high, padj < 0.05) %>% 
    ggplot(aes(x = phase_assigned %in% acr, y = log(znf274_znfs))) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, outlier.size = 0.5) + 
    stat_compare_means(label.x = 0.9, label.y = 11, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("ZNF274 ChIP-seq") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    coord_cartesian(ylim = c(2, 11),
                      clip = 'off')
p
```

Still signif

```R
fit_kzfps %>% nrow()
```

<!-- #raw -->
# same for rhythmic TFs amongst rhythmic genes

k = fit_kzfps %>% nrow() # number of kzfpS

#acr = c("G1/S", "S1", "S2", "S/G2", "G2")

acr = c( "S2", "S/G2", "G2")

m = fit_qualitative %>% dplyr::filter(phase_assigned %in% acr) %>% nrow() # number of kzfps in S-G2

n = fit_qualitative %>% dplyr::filter(!phase_assigned %in% acr) %>% nrow() # number of kzfps in M-G1/S

x = 0:m # is the variable tested

x_observed = plot_selection %>% dplyr::filter(phase_assigned %in% acr) %>% nrow() # number of kzfps in S-G2, with median replitiming minimum

probs <- dhyper(x, m, n, k, log = FALSE)

# we compute the probability of observing a more extreme depletion, therefore using a one sided test. 

pval_TFs_vs_genes_one_sided = 1-sum(probs[x<=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_TFs_vs_genes_two_sided = sum(probs[probs <= pval_TFs_vs_genes_one_sided])
pval_TFs_vs_genes_two_sided

x_observed
k
m
<!-- #endraw -->

<!-- #raw -->
plot_selection$tf
<!-- #endraw -->

<!-- #raw -->
n_tf_left = plot_selection %>% dplyr::filter(phase_assigned %in% acr, tf_category == "TF") %>% nrow()
n_kzfp_left = plot_selection %>% dplyr::filter(phase_assigned %in% acr, tf_category == "KZFP") %>% nrow()
n_C2H2_left = plot_selection %>% dplyr::filter(phase_assigned %in% acr, tf_category == "C2H2 ZF") %>% nrow()

n_tf_right = plot_selection %>% dplyr::filter(!phase_assigned %in% acr, tf_category == "TF") %>% nrow()
n_kzfp_right = plot_selection %>% dplyr::filter(!phase_assigned %in% acr, tf_category == "KZFP") %>% nrow()
n_C2H2_right = plot_selection %>% dplyr::filter(!phase_assigned %in% acr, tf_category == "C2H2 ZF") %>% nrow()
<!-- #endraw -->

<!-- #raw -->
p = plot_selection %>% dplyr::rename(`Gene category`=tf_category) %>% ggplot(aes(ifelse(phase_assigned %in% acr, yes = "S2-G2", no = "M-S1"), fill = `Gene category`)) + 
    geom_bar() + 
    ggplot2::scale_fill_manual(values = cols$`Gene category`) + 
    theme_classic() + 
    ylab("Number of rhythmic TFs (p-adj < 0.05)") + 
    xlab("acrophase") + 
    annotate(geom="text", x=2, y=n_tf_right+n_kzfp_right/2, label=paste0("KZFP p-val=", formatC(pval_kzfp_vs_TFs_two_sided, format = 'e', digits=2))) + 
    annotate(geom="text", x=2, y=n_tf_right+n_kzfp_right + n_C2H2_right/2, label=paste0("C2H2-ZF p-val=", formatC(pval_c2h2_vs_TFs_two_sided, format = 'e', digits=2)))
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_enrichment.pdf", p, pdf, height = 5, width = 5)
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_enrichment.svg", p, svg, height = 4, width = 5)

<!-- #endraw -->

```R
# ranking replication timing per cluster
p = fit_kzfps %>% dplyr::filter(!is.na(TEs.pval)) %>% 
    group_by(cluster) %>% dplyr::select(RepliTiming, phase_assigned, cluster)  %>% 
    dplyr::mutate(min_RepliTiming_per_cluster = min(as.integer(RepliTiming))) %>% dplyr::mutate(repl_timing_is_cluster_min = as.integer(RepliTiming) == min_RepliTiming_per_cluster) %>%

    ggplot(aes(x = phase_assigned %in% acr,fill = repl_timing_is_cluster_min)) + geom_bar() + scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + ylab("# KZFPs") + xlab("") + theme_classic() + scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "lightgrey"))
p
```

```R
ggsave("../out/cell_cycle_figures/KFZPs_min_replitiming_vs_phase.svg", p, svg, height = 3, width = 4)
```

```R
fit_kzfps %>% dplyr::filter(phase_assigned %in% acr) %>% pull(RTindex_avg) %>% mean()
```

```R
# ranking replication timing per cluster
p = fit_kzfps %>% 
    ggplot(aes(x = phase_assigned %in% acr,y = RTindex_avg)) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr) %>% pull(RTindex_avg) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("RT index") + 
    xlab("") + 
    theme_classic()
p
```

```R
ggsave("../out/cell_cycle_figures/KFZPs_replitiming_vs_phase_continuous.svg", p, svg, height = 3, width = 3)
```

```R
# ranking replication timing per cluster
p = fit_kzfps %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = phase_assigned %in% acr,y = RTindex_avg)) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(RTindex_avg) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("RT index") + 
    xlab("") + 
    theme_classic()
p
```

```R
ggsave("../out/cell_cycle_figures/KFZPs_rhythmic_replitiming_vs_phase_continuous.svg", p, svg, height = 3, width = 2.5)
```

Still signif

```R
# ranking replication timing per cluster
p = fit_kzfps %>% 
    ggplot(aes(x = phase_assigned %in% acr,y = rif1)) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr) %>% pull(RTindex_avg) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("RT index") + 
    xlab("") + 
    theme_classic()
p
```

Promoter activity of rhythmic KZFPs

```R
# promoter h3k4me3
p = fit_kzfps %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = phase_assigned %in% acr,y = log(znf274_znfs + 1))) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(znf274_znfs) %>% (function(x) log(x+1)) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("RT index") + 
    xlab("") + 
    theme_classic()
p
```

```R
# promoter h3k4me3
p = fit_kzfps %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = phase_assigned %in% acr,y = log(h3k4me3 + 1))) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(h3k4me3) %>% (function(x) log(x+1)) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("H3K4me3") + 
    xlab("") + 
    theme_classic()
p
```

```R
# promoter h3k4me3
p = fit_kzfps %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = phase_assigned %in% acr,y = log(atac + 1))) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(atac) %>% (function(x) log(x+1)) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("H3K4me3") + 
    xlab("") + 
    theme_classic()
p
```

```R
# promoter h3k4me3
p = fit_kzfps %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = phase_assigned %in% acr,y = log(h3k27ac + 1))) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(h3k27ac) %>% (function(x) log(x+1)) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("H3K4me3") + 
    xlab("") + 
    theme_classic()
p
```

```R
# promoter h3k4me3
p = fit_kzfps %>% 
    dplyr::filter(padj < 2) %>%
    ggplot(aes(x = phase_assigned %in% acr,y = log(h3k9me3_znfs + 1))) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 2) %>% pull(h3k9me3_znfs) %>% (function(x) log(x+1)) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("H3K4me3") + 
    xlab("") + 
    theme_classic()
p
```

```R
# promoter h3k4me3
p = fit_kzfps %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = phase_assigned %in% acr,y = log(h3k9me3_znfs + 1))) + 
    geom_violin() + 
    geom_jitter(size = 0.5, height = 0, width = 0.1) + 
geom_boxplot(width=0.1, outlier.size = 0.5, alpha = 0) +
geom_hline(yintercept = fit_kzfps %>% dplyr::filter(!phase_assigned %in% acr, padj < 0.05) %>% pull(h3k9me3_znfs) %>% (function(x) log(x+1)) %>% median(), lty = 2) + 
ggpubr::stat_compare_means(label.y = 5.1) +
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2")) + 
    ylab("H3K4me3") + 
    xlab("") + 
    theme_classic()
p
```

```R
# rif1 vs RT

p = fit_kzfps %>% 
    ggplot(aes(x = log(znf274_znfs+1), y = RTindex_avg)) + 
    geom_point() +
    ggpubr::stat_cor(method = "spearman") +
    theme_classic()
p

```

```R
p = fit_kzfps %>% 
    dplyr::filter(!cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1")) %>%
    ggplot(aes(x = log(znf274_znfs+1), y = RTindex_avg)) + 
    geom_point() +
    ggpubr::stat_cor(method = "spearman") +
    theme_classic()
p
```

```R
# rif1 vs RT

p = fit_kzfps %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = log(znf274_znfs+1), y = RTindex_avg)) + 
    geom_point() +
    ggpubr::stat_cor(method = "spearman") +
    theme_classic()
p

```

```R
# rif1 vs RT

p = fit_kzfps %>% 
    dplyr::filter(padj < 0.05) %>%
    dplyr::filter(!cluster %in% c("chr19.8", "chr19.9", "chr19.11", "chr6.1")) %>%
    ggplot(aes(x = log(znf274_znfs+1), y = RTindex_avg)) + 
    geom_point() +
    ggpubr::stat_cor(method = "spearman") +
    theme_classic()
p

```

```R

```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval)) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(promoter.pval < 0.05, yes = T, no = ifelse(exons.pval < 0.05, yes = T, no = F)))) + geom_bar(position="fill")

```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval)) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(promoter.pval < 0.05, yes = T, no = F))) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval), padj < 0.05) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(promoter.pval < 0.05, yes = T, no = F))) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval)) %>% dplyr::arrange(promoter.pval) %>% dplyr::select(phase_assigned, padj, symbol) %>% head(20)
```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval), padj < 0.05) %>% dplyr::arrange(promoter.pval) %>% head(10) %>% ggplot(aes(phase_assigned %in% acr)) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval)) %>% ggplot(aes(ifelse(padj < 0.05, yes = T, no = F), fill = ifelse(promoter.pval < 0.05, yes = T, no = F))) + geom_bar(position = "fill")
```

```R
fit_kzfps %>% dplyr::group_by(promoter.pval <0.0005, padj < 0.05) %>% summarize(count = n())
```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval), padj < 0.05) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(TEs.pval < 0.05, yes = T, no = F))) + geom_bar(position="fill")

```

```R

```

```R
fit_kzfps %>% dplyr::filter(padj < 0.05, !is.na(age_MA)) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(age_MA >= 105, yes = T, no = F))) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% dplyr::filter(!is.na(age_MA)) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(age_MA >= 105, yes = T, no = F))) + geom_bar(position="fill")

```

Old KZNFs are enriched with KZFNs peaking in S phase

```R
fit_kzfps %>% dplyr::filter(age_MA < 105, !is.na(promoter.pval)) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(promoter.pval < 0.05, yes = T, no = F))) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% dplyr::filter(age_MA < 105, !is.na(promoter.pval)) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(promoter.pval < 0.05, yes = T, no = ifelse(exons.pval < 0.05, yes = T, no = F)))) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(age_MA >= 105, yes = T, no = F))) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval), padj < 0.05) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(promoter.pval < 0.05, yes = T, no = F))) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval)) %>% ggplot(aes(phase_assigned %in% acr, fill = ifelse(promoter.pval < 0.05, yes = T, no = ifelse(exons.pval < 0.05, yes = T, no = ifelse(introns.pval < 0.05, yes = T, no = F))))) + geom_bar(position = "fill")

```

```R
fit_kzfps %>% dplyr::filter(!is.na(TEs.pval)) %>% ggplot(aes(phase_assigned %in% acr, fill = rif1>0)) + geom_bar(position = "fill")
```

```R
plot_selection %>% dplyr::filter(tf_category == "KZFP")  %>% ggplot(aes(phase_assigned %in% acr, fill = rif1>0)) + geom_bar()
```

```R
plot_selection %>% dplyr::filter(tf_category == "TF")  %>% ggplot(aes(phase_assigned %in% acr, fill = rif1>0)) + geom_bar(position = "fill")
```

```R
plot_selection %>% dplyr::filter(tf_category == "C2H2 ZF")  %>% ggplot(aes(phase_assigned %in% acr, fill = rif1>0)) + geom_bar()
```

```R
acr = c( "S2", "G2", "G2/M")
plot_selection %>% dplyr::filter(tf_category == "C2H2 ZF")  %>% ggplot(aes(x = phase_assigned %in% acr, y = rif1)) + geom_boxplot() + geom_point()
```

```R
plot_selection %>% dplyr::filter(tf_category == "KZFP")  %>% ggplot(aes(x = phase_assigned %in% acr, y = rif1)) + geom_boxplot() + geom_point()
```

```R
plot_selection %>% dplyr::filter(tf_category == "TF")  %>% ggplot(aes(x = phase_assigned %in% acr, y = rif1)) + geom_boxplot() + geom_point()
```

```R
plot_selection %>% dplyr::filter(tf_category != "KZFP")  %>% ggplot(aes(x = phase_assigned %in% acr, y = rif1)) + geom_boxplot() + geom_point()
```

```R
repliTimingFactor = factor(plot_selection$RepliTiming, levels = (plot_selection$RepliTiming %>% unique() %>% sort()))
```

```R
plot(plot_selection$phase_assigned, repliTimingFactor)
```

```R
plot_selection$acrophase_rad = plot_selection$acrophase/8*(2*pi)
plot_selection$acrophase_rad %>% range()
```

```R
plot_selection %>% dplyr::group_by(tf_category) %>% summarize(avg_sin = mean(sin(acrophase_rad)))
```

```R
plot_selection %>% dplyr::group_by(tf_category) %>% summarize(avg_cos = mean(cos(acrophase_rad)))
```

```R
row_anno = plot_selection %>% dplyr::select(mesor)
ha = rowAnnotation(mesor = anno_points(row_anno))

plot_order = plot_selection$ensembl
length(plot_order)
to_plot = expr_zscore[plot_order, ]
rownames(to_plot) = (gene_metadata %>% column_to_rownames('ensembl'))[plot_order, 'symbol']
colnames(to_plot) = new_phases
Heatmap(to_plot, 
        name = 'z-scored_log2_norm.adj.counts.plus.1', 
        cluster_rows = FALSE, 
        cluster_columns = FALSE,
        right_annotation = ha,
       row_names_gp = grid::gpar(fontsize = 7))
```

KZFP seem even more restricted to a peak of expression in eG1 than C2H2, more of which (relatively) peak in S phase.

If the argument of the level of expression to explain everything comes in, we can do a permutation test, looking at the chance of obtaining a distribution of acrophases centered on S just by sampling similar MESOR distributions as the KZFPs.

```R

```

```R

```
