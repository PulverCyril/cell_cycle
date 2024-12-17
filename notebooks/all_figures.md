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

# Evolutionarily recent transcription factors partake in human cell cycle regulation

## All figures

```R
library(tidyverse)
library(ComplexHeatmap)
library(magick) # rasterization of ComplexHeatmap matrices
library(ggpubr)
library(circlize)
library(viridis)
library(spatstat) # to smooth the H3K9me3 signal over rhythmic promoters
library(EBImage) # to smooth the H3K9me3 signal over rhythmic promoters
library(clusterProfiler) # GO enrichment analysis
library(GOSemSim) # GO enrichment analysis
library(org.Hs.eg.db) # for enrichment to work
library(cowplot) # multi-panel plots
library(ggnewscale) # two different color codes
```

```R
rhythm = read.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', header = T)
rhythm %>% dim()
rhythm %>% colnames()

# we remove RTindex_avg, RepliTiming (discrete, not needed anymore) and rif1, which we replace using the gene_metadata_augmented table
rhythm = rhythm %>% dplyr::select(-RepliTiming, -RTindex_avg, -rif1)
```

```R
gene_metadata = read.table('../out/tables/gene_metadata_augmented.tsv', sep = '\t', header = T, quote = '')
gene_metadata %>% dim()
gene_metadata %>% colnames()
```

```R
# sanity check on ZNF695 rhythmicity
rhythm %>% dplyr::filter(symbol == "ZNF695") %>% dplyr::select(stars)
```

# Figure 1

```R
rhythm %>% colnames()
```

```R
rhythm %>% 
    dplyr::left_join(., gene_metadata) %>%
    dplyr::filter(padj < 0.05) %>% 
    arrange(acrophase) %>% head()
```

```R
# exporting the coordinates of rhythmic genes TSS sorted by acrophase as a bed file,
# to retrieve high definition H3K9me3 signal 
rhythm %>% 
    dplyr::left_join(., (gene_metadata %>% dplyr::select(ensembl, start_TSS, end_TSS))) %>%
    dplyr::filter(padj < 0.05) %>% 
    arrange(acrophase) %>% 
    dplyr::select(chr, start_TSS, end_TSS, ensembl) %>%
    write.table('../out/temp/tss_rhythmic_genes_sorted_acrophase.bed', col.names = F, row.names = F, sep = '\t', quote = F)
```

```R
# Computing and adding the width or intensity of the H3K9me3 trough at the promoter of rhythmic genes
h3k9me3_promoters_hires = read.table('../out/temp/tss_rhythmic_genes_sorted_acrophase_h3k9me3_domains.tab', sep = '\t', skip = 3)
```

```R
which(is.na(h3k9me3_promoters_hires))
```

```R
h3k9me3_promoters_hires %>% head(10)
```

```R
h3k9me3_promoters_hires %>% dim()
```

```R
stopifnot(nrow(h3k9me3_promoters_hires)==(nrow(rhythm %>% dplyr::filter(padj < 0.05))))
```

```R
# plotting the raw H3K9me3 signal at promoters
log(h3k9me3_promoters_hires+1) %>% as.matrix() %>% ComplexHeatmap::Heatmap(cluster_rows = F, 
                                                                    cluster_columns = F,
                                                                    col= colorRamp2(c(0, 0.01, 0.5, 1), viridis::magma(4)))
```

```R
# let's attempt smoothing the signal to observe a more general trend
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
cutoff = median(mat_smoothed[, 85:215] %>% as.vector())
cutoff
```

```R
h3k9me3_signal_tss = mat_smoothed[, 85:215] %>% apply(., 1, function(x) sum(x > cutoff))/length(85:215)
```

```R
plot(1:nrow(mat_smoothed), h3k9me3_signal_tss)
```

```R
plot_selection = rhythm %>% dplyr::filter(padj < 0.05) %>% 
    dplyr::arrange(acrophase) %>% 
    dplyr::left_join(., gene_metadata)
```

```R
# setting up the legend of the plot
n_rhythmic = plot_selection %>% nrow()
```

```R
paste0(n_rhythmic, ' rhythmic genes (padj < 0.05)')
```

```R
plot_order = (plot_selection %>% dplyr::arrange(acrophase))$ensembl
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "expr z-scored log2", direction = "horizontal")
```

```R
plot_selection$RIF1_tss_body %>% summary()
```

```R
plot_selection$RTindex_K562 %>% summary()
```

```R
cols = list(`bas. expr.` = colorRamp2(breaks = c(-2.5, 2.5, 4, 6, 10), colors = viridis::mako(5)),
    RIF1 =  colorRamp2(c(-0.5, -0.2, 0, 0.2, 0.5), hcl.colors(n = 5, palette = "Cividis", rev = F)),
    RT = colorRamp2(breaks = seq(from = 2.5, to = 4.5, length.out = 4), colors = viridis::viridis(4)),
    H3K9me3 = colorRamp2(breaks = seq(from = 0.4, to = 0.8, length.out = 10), colors = viridis::magma(10)))
```

```R
# legacy: RepliTiming was added to row_anno in the past `Repl. timing` = RepliTiming,

row_anno = plot_selection %>% dplyr::select(RTindex_K562, RIF1_tss_body, mesor) %>% 
                        dplyr::rename(RIF1 = RIF1_tss_body, `bas. expr.` = mesor, RT = RTindex_K562)
row_anno$H3K9me3 = h3k9me3_signal_tss
```

```R
row_anno %>% head()
```

```R
# annotating the top 10 rhythmic gene per phase merged with the top 10 rhythmic genes overall
genes_to_label = c(rhythm %>% 
    dplyr::filter(padj < 0.05, !phase_assigned %in% c()) %>% 
    dplyr::group_by(phase_assigned) %>%
    dplyr::arrange(padj) %>% dplyr::slice_head(n=1) %>% dplyr::pull(symbol),
  rhythm %>%
      dplyr::arrange(padj) %>% head(10) %>% pull(symbol) %>% unique())
```

```R
genes_to_label_idx = match(genes_to_label, (rhythm %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase) %>% pull(symbol)))
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
# loading the expression values for the heatmap
expr_zscore = read.table('../out/tables/expr_zscore_RNAseq_cellcycle_romain.tsv', sep = '\t', header = T) %>% as.matrix()
expr_zscore %>% head()
```

```R
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
        annotation_height = max_text_width(cn)),
              use_raster = T,
               raster_quality = 5,
               raster_by_magick = TRUE) %v% NULL
```

```R
svg('../out/cell_cycle_figures/cycling_genes_heatmap_with_RepliTiming.svg', height = 5.5, width = 6)
draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()

```

## Adding genes to label for GeneCards

```R
rhythm %>% dplyr::filter(symbol == "PIK3C3") %>% dplyr::pull(acrophase)
```

```R
acr_threshold_min = rhythm %>% dplyr::filter(symbol == "PIK3C3") %>% dplyr::pull(acrophase)
acr_threshold_max = rhythm %>% dplyr::filter(symbol == "METTL7A") %>% dplyr::pull(acrophase)


```

```R
supplementary_genes_to_label = c(rhythm %>% 
    dplyr::filter(padj < 0.05, acrophase < acr_threshold_min) %>% 
    dplyr::group_by(phase_assigned) %>%
    dplyr::arrange(padj) %>% dplyr::slice_head(n=3) %>% dplyr::pull(symbol),
  rhythm %>% 
    dplyr::filter(padj < 0.05, acrophase > acr_threshold_max) %>% 
    dplyr::group_by(phase_assigned) %>%
    dplyr::arrange(padj) %>% dplyr::slice_head(n=3) %>% dplyr::pull(symbol))
```

```R
genes_to_label_idx = match(c(genes_to_label, supplementary_genes_to_label), (rhythm %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase) %>% pull(symbol)))
```

```R
genes_to_label_idx
```

```R
ha = rowAnnotation(df = row_anno, col = cols)
label_annotation = ComplexHeatmap::rowAnnotation(foo = anno_mark(at = genes_to_label_idx, side = "left", padding = 1,
    labels = c(genes_to_label, supplementary_genes_to_label)))
plot_order = plot_selection$ensembl
```

```R
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
        annotation_height = max_text_width(cn)),
              use_raster = T,
               raster_quality = 5,
               raster_by_magick = TRUE) %v% NULL
```

```R
svg('../out/cell_cycle_figures/cycling_genes_heatmap_with_RepliTiming_supplementary_labels_genecards.svg', height = 5.5, width = 6)
draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.off()

```

## Baseline expression, RIF1, RT and H3K9me3 at rhythmic promoters

```R
# RIF1 M-S1 vs. S2-G2/M
acr = c("S2", "G2", "G2/M")
acr_excluded = c()

y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(RIF1_tss_body) %>% 
    median(na.rm = T)

p = plot_selection %>%
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::select(phase_assigned, RIF1_tss_body) %>% 
    dplyr::rename(RIF1 = RIF1_tss_body) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "RIF1", add = "boxplot") + 
    geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 2.4, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("RIF1 [log2 FC over input]") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2/M")) + 
    coord_cartesian(ylim = c(-2, 2.5),
                      clip = 'off')

svg('../out/cell_cycle_figures/RIF1_stat_test.svg', height = 2.5, width = 2)
p
dev.off()
```

```R
# MESOR (baseline expression) M-S1 vs. S2-G2/M
acr = c("S2", "G2", "G2/M")
acr_excluded = c()

y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(mesor) %>% 
    median(na.rm = T)

p = plot_selection %>%
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::select(phase_assigned, mesor) %>% 
    dplyr::rename(`bas. expr.` = mesor) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "bas. expr.", add = "boxplot") + 
    geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 14, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("bas. expr. [log norm. counts]") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2/M")) + 
    coord_cartesian(ylim = c(-5, 15),
                      clip = 'off')

svg('../out/cell_cycle_figures/bas_expr_stat_test.svg', height = 2.5, width = 2)
p
dev.off()
```

```R
plot_selection$RT
```

```R
# RT K562 M-S1 vs. S2-G2/M
acr = c("S2", "G2", "G2/M")
acr_excluded = c()

y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(RTindex_K562) %>% 
    median(na.rm = T)

p = plot_selection %>%
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::select(phase_assigned, RTindex_K562) %>% 
    dplyr::rename(`RT` = RTindex_K562) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "RT", add = "boxplot") + 
    geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 7, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("RT [RTindex]") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2/M")) + 
    coord_cartesian(ylim = c(-6, 8),
                      clip = 'off')

svg('../out/cell_cycle_figures/RT_stat_test.svg', height = 2.5, width = 2)
p
dev.off()
```

```R
# H3K9me3 TSS + body K562 M-S1 vs. S2-G2/M
acr = c("S2", "G2", "G2/M")
acr_excluded = c()

y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(H3K9me3_tss_body) %>% 
    median(na.rm = T)

p = plot_selection %>%
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::select(phase_assigned, H3K9me3_tss_body) %>% 
    dplyr::mutate(`H3K9me3` = log(H3K9me3_tss_body+1)) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "H3K9me3", add = "boxplot") + 
   # geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 7, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("H3K9me3") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2/M"))
    #coord_cartesian(ylim = c(5, 10),
     #                 clip = 'off')

svg('../out/cell_cycle_figures/H3K9me3_tss_body_stat_test.svg', height = 2.5, width = 2)
p
dev.off()
```

```R
# H3K9me3 at promoters, counting the number of bases where the signal is more than a certain cutoff, K562 M-S1 vs. S2-G2/M

# mat contains the raw log(H3K9me3+1) values as a matrix on the promtoers
# im.med contains the median-filtered values
```

```R
mat_median_filtered = im.med %>% as.matrix() %>%.[(nrow_padding+1):(nrow(mat)+nrow_padding), ]
```

```R
cutoff = mean((mat_median_filtered[, 85:215] %>% as.vector()))
cutoff
```

```R
1:4
```

```R
h3k9me3_signal_tss = mat_median_filtered[, 85:215] %>% apply(., 1, function(x) sum(x > cutoff))/(length(85:215))
```

```R
plot(1:nrow(mat), h3k9me3_signal_tss)
```

```R
plot_selection$H3K9me3_tss_processed = h3k9me3_signal_tss
```

```R
acr = c("S2", "G2", "G2/M")
acr_excluded = c()

y_dashed = plot_selection %>% 
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>% 
    dplyr::filter(!phase_in_S2_M) %>% 
    dplyr::pull(H3K9me3_tss_processed) %>% 
    median(na.rm = T)

p = plot_selection %>%
    dplyr::filter(!phase_assigned %in% acr_excluded) %>%
    dplyr::select(phase_assigned, H3K9me3_tss_processed) %>% 
    dplyr::mutate(`H3K9me3` = H3K9me3_tss_processed) %>% 
    dplyr::mutate(phase_in_S2_M = phase_assigned %in% acr) %>%
    ggpubr::ggviolin(x = "phase_in_S2_M", y = "H3K9me3", add = "boxplot") + 
   geom_hline(yintercept = y_dashed, linetype = "dashed") +
    stat_compare_means(label.x = 0.9, label.y = 1.1, step.increase = 1) +
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("H3K9me3 [frac. of prom. covered]") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("M-to-S1", "S2-to-G2/M"))+
    coord_cartesian(ylim = c(0, 1.2),
                      clip = 'off')

svg('../out/cell_cycle_figures/H3K9me3_tss_processed_stat_test.svg', height = 2.5, width = 2.5)
p
dev.off()
```

```R
# exporting for romain
plot_selection %>% 
    dplyr::select(ensembl, symbol, H3K9me3_tss_processed) %>% 
    write.table('../out/tables/H3K9me3_tss_median_filtered_frac_covered.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

## Age of TFs assessed for enrichment at rhythmic promoters

```R
gene_metadata %>% head()
```

```R
# assigning the gene age to 1) Imbeault, 2) genTree (vertebrates), 3) genOrigin (for older than vertebrates)
```

```R
# loading enrichment data
in_phase_vs_all_genes_enrichment_long = read.table('../out/tables/TF_binding_enrichment_in_rhythmic_genes_from_each_phase.tsv', sep = '\t', header = T, quote = '')
in_phase_vs_all_genes_enrichment_long %>% head()
in_phase_vs_all_genes_enrichment_long %>% dim()
in_phase_vs_all_genes_enrichment_long %>% dplyr::select(TF_symbol) %>% unique() %>% dim()
```

```R
gene_metadata_TF = gene_metadata %>% dplyr::filter(symbol %in% in_phase_vs_all_genes_enrichment_long$TF_symbol)
```

```R
in_phase_vs_all_genes_enrichment_long = in_phase_vs_all_genes_enrichment_long %>% 
    dplyr::mutate(promoters_active_in_phase = factor(promoters_active_in_phase, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")))
```

```R
# Exporting enriched DBPs per phase of the promoters they bind, for exporting signatures
dir.create('../out/gene_lists_database_submissions/')
for (p in (in_phase_vs_all_genes_enrichment_long$promoters_active_in_phase %>% unique)) {
    in_phase_vs_all_genes_enrichment_long %>% 
        dplyr::filter(padj_enrich < 0.05, promoters_active_in_phase == p) %>% 
        dplyr::select(TF_symbol) %>% 
        unique() %>%
        dplyr::arrange(TF_symbol) %>%
        write.table(., paste0('../out/gene_lists_database_submissions/', gsub("/", "_", p), '_enriched_DBPs_symbol.txt'), quote = F, sep = '\t', eol = ',', row.names=F, col.names = F)
    }

# all TFs enriched, concatenated
in_phase_vs_all_genes_enrichment_long %>% 
    dplyr::filter(padj_enrich < 0.05) %>% 
    dplyr::select(TF_symbol) %>% 
    unique() %>%
    dplyr::arrange(TF_symbol) %>%
    write.table(., paste0('../out/gene_lists_database_submissions/all_enriched_DBPs_symbol.txt'), quote = F, sep = '\t', eol = ',', row.names=F, col.names = F)
```

```R
#gene_metadata_TF$age_category = cut(gene_metadata_TF$age_combined, breaks = c(-Inf, 25, 50, 100, 250, 500, 1000, 2000, Inf))
```

```R
gene_metadata_TF$age_category = cut(gene_metadata_TF$age_combined, breaks = c(-Inf, 50, 250, 1000, Inf))

```

```R
gene_metadata_TF$age_category %>% levels()
```

```R
gene_metadata_TF %>% dplyr::filter(symbol == "ZNF8")
```

```R
gene_metadata_TF %>% dim()
gene_metadata_TF = gene_metadata_TF %>% dplyr::filter(!((symbol == "ZNF8") & (ensembl == "ENSG00000273439")))
gene_metadata_TF %>% dim()
```

```R
to_plot = in_phase_vs_all_genes_enrichment_long %>% dplyr::left_join(., gene_metadata_TF %>% 
    dplyr::select(symbol, age_combined, tf_category, isC2H2, isKZFP, age_category) %>%
    unique(), by = c("TF_symbol" = "symbol")) %>% 
    dplyr::filter(!age_category == "NA") %>%
    dplyr::group_by(TF_symbol) %>% 
   slice_max(specificity, n = 1, with_ties = F)
```

```R
to_plot %>% head()
```

```R
to_plot$block = ifelse(to_plot$age_combined <=1000, yes = "A", no = "B")

to_plot %>%
ggplot(aes(x = age_combined, fill = tf_category)) + 
geom_histogram(alpha = 0.8) +
#+ facet_grid(.~block, scales = "free_x", space = "free_x", ) +
theme_classic() + 
labs(x = "Emergence [MYA]", y = "DNA-binding proteins", title = paste0(as.character(nrow(to_plot)), " DNA-binding proteins\nwith time of emergence")) + 
scale_fill_manual(values = c("green3", "#FF681F", "grey", "black"))
dev.copy(svg, '../out/cell_cycle_figures/KZFP_age_histogram.svg', height = 3, width = 4)
dev.off()
```

## Over-representation of old TFs in phase specific TFs

```R
for_enrichment = in_phase_vs_all_genes_enrichment_long %>% dplyr::left_join(., gene_metadata_TF %>% 
    dplyr::select(symbol, age_combined, tf_category, isC2H2, isKZFP, age_category) %>%
    unique(), by = c("TF_symbol" = "symbol")) %>% 
    dplyr::filter(!age_category == "NA") %>%
    dplyr::group_by(TF_symbol) %>% 
    slice_min(pval_enrich_min, n = 1, with_ties = F)
for_enrichment %>% head()
for_enrichment %>% dim()
```

```R
age_threshold = 450
```

```R
x_observed = for_enrichment %>% dplyr::filter(age_combined > age_threshold, pval_enrich_min <= 0.05) %>% nrow()
m = for_enrichment %>% dplyr::filter(pval_enrich_min <= 0.05) %>% nrow()
n = for_enrichment %>% dplyr::filter(pval_enrich_min > 0.05) %>% nrow()

k = for_enrichment %>% dplyr::filter(age_combined > age_threshold) %>% nrow()

x = 0:m # is the variable tested
probs <- dhyper(x, m, n, k, log = FALSE)


# we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
pval_one_sided = sum(probs[x>=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])

pval_two_sided
```

```R
in_phase_vs_all_genes_enrichment_long %>% dplyr::left_join(., gene_metadata_TF %>% 
    dplyr::select(symbol, age_combined, tf_category, isC2H2, isKZFP, age_category) %>%
    unique(), by = c("TF_symbol" = "symbol")) %>% 
    dplyr::filter(!age_category == "NA") %>%
    dplyr::group_by(TF_symbol) %>% 
    slice_min(pval_enrich_min, n = 1, with_ties = F) %>%
ggplot(aes(x = age_combined > age_threshold, fill = pval_enrich_min < 0.05)) + 
geom_bar() +
# legend
labs(title = paste0("p = ", format(pval_two_sided, scientific = T, digits = 2), ", Fisher's Exact Test"), x = "Emergence", y = "Number of TFs", fill = "TF enriched at rhythm.\nprom. of \u2265 1 phase\n") +
scale_fill_manual(labels = c("FALSE", "TRUE"), values = c("grey", "black")) +
scale_x_discrete(labels=c(paste0("<", age_threshold,"MYA"), paste0(">", age_threshold, "MYA"))) + 
theme_classic()
dev.copy(svg, '../out/cell_cycle_figures/TF_enrichment_at_rhythmic_proms_vs_age_Fisher_test.svg', height = 2.2, width = 3.8)
dev.off()
```

## Age of TFs vs enrichment for binding promoters of rhythmic genes

```R
to_plot %>% head()
```

```R
to_plot$age_category = cut(to_plot$age_combined, breaks = c(-Inf, 50, 250, 1000, Inf))
```

```R
to_plot %>% ggplot(aes(x = age_category, y = log2(specificity))) + 
    geom_point(aes(col = isKZFP), alpha = 0.8) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, outliers = F) + 
    ggpubr::stat_compare_means(comparisons = list(c(1, 2), c(2, 3), c(3, 4)), method = "wilcox.test") + 
    theme_classic() + 
    scale_color_manual(values = c("black", "#FF681F")) + 
    ggrepel::geom_text_repel(aes(x = age_category, y = log2(specificity), label = TF_symbol, col = is_kzfp),
                              data = to_plot %>% group_by(age_category) %>% slice_max(specificity, n = 6),
                              label.size = 0, fill = alpha(c("white"),0), max.overlaps = 10, seed = 1, force = 5) +
    labs(x = "Emergence", y = "Phase specificity [log]", col = "KZFP") +
    scale_x_discrete(labels=c("0-50MYA", "51-250MYA", "251-1000MYA", "\u2265 1001MYA"))
   

dev.copy(svg, '../out/cell_cycle_figures/specificity_vs_age.svg', height = 4, width = 5.5)
dev.off()
```

## TF binding enrichment at promoters of rhythmic genes


### Focus on G1/S

```R
## x axis: -log10(padj) instead of rank
to_plot = in_phase_vs_all_genes_enrichment_long %>% 
    dplyr::filter(promoters_active_in_phase %in% c("G1/S")) %>%
    dplyr::filter(padj_enrich < 0.05) %>%
    # for labelling KZFPs better
    dplyr::group_by(promoters_active_in_phase) %>% 
    dplyr::mutate(thresh_label_padj = sort(padj_enrich)[8],
                 thresh_label_spec = sort(specificity, decreasing = T)[4])  %>%
    dplyr::mutate(TF_symbol_label = ifelse((padj_enrich<=thresh_label_padj | specificity >= thresh_label_spec | promoters_active_in_phase == "eG1" | is_kzfp | TF_symbol %in% c("RB1", "TFDP1", "E2F1", "E2F2", "E2F3", "E2F4")), TF_symbol, "")) %>%
    dplyr::ungroup()
```

```R
p = to_plot %>%
    ggplot(aes(x = -log10(padj_enrich), y = log2(specificity), label = TF_symbol, col = is_kzfp, size = n_promoters_bound)) + 
    #ggrastr::geom_point_rast(alpha = 0.8) + 
    geom_point(alpha = 0.8) + 
    ggrepel::geom_text_repel(aes(x = -log10(padj_enrich), y = log2(specificity), label = TF_symbol_label, col = is_kzfp),
                              min.segment.length = 0, segment.alpha = 0.8, size = 3.5, label.size = 0, fill = alpha(c("white"),0), max.overlaps = 50, box.padding = 0.3, seed = 2, force = 30) + 
    ylab("Phase specificity [log]") + 
    xlab("-log10(adj. p)") +
    theme_classic() +
    #theme(legend.position = "none") + 
    guides(color = "none") +
    scale_size(name = "Nr. bound prom.", breaks = c(10, 100, 1000, 10000), range = c(1, 3.5)) +
    scale_color_manual(values = c("black", "#FF681F"))
```

```R
p
```

```R
ggsave('../out/cell_cycle_figures/specialization_vs_log10padj_all_TFs_G1_S.svg', p, svg, width = 4.5, height = 3.5)
```

```R
to_plot %>% arrange(desc(specificity)) %>% head(10)
```

```R
in_phase_vs_all_genes_enrichment_long %>% dplyr::filter(promoters_active_in_phase == "G1/S", TF_symbol == "ZNF519")
```

# Fig. 1 supp


## CCNA FACS vs expr.: Romain

## E2F1 expression throughout the eight gates

```R
rhythm_for_plot = read.table('../out/tables/rhythm_for_plot.tsv', sep = '\t', header = T, quote = '')
rhythm_for_plot %>% head()
rhythm_for_plot %>% dplyr::select(symbol) %>% unique() %>% dim()
```

```R
plot_rhythmic_expression <- function(g) {
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
    return(p)
}
```

```R
plot_rhythmic_expression("E2F1")
dev.copy(svg, filename = '../out/cell_cycle_figures/E2F1_rhythmic_expr.svg', width = 4, height = 4)
dev.off()
```

## Correlating phases estimated by Velocycle and the cosinor

```R
sc = read.csv('../data/velocycle_output//cyril_cycle_gene_harmonics_metadata.csv', header = T)
colnames(sc)[1] = "symbol"
sc %>% head()
sc %>% dim()
```

```R
merged = dplyr::full_join(rhythm %>% 
                              dplyr::mutate(acrophase_angular = acrophase*(2*pi/8)), 
                          sc, by = c("symbol"), suffix = c("_cosinor", "_velocycle"))
merged %>% dim()
```

```R
amplitude_threshold = merged$amplitude_velocycle %>% quantile(., probs = 0.90, na.rm = TRUE)
amplitude_threshold
```

Function to compute the direction (mean angle) of several angles: 
![image.png](attachment:74c231d8-2c3f-437e-897f-de253c449990.png)

```R
mean_direction <- function(alphas) {
    return(atan2(mean(sin(alphas)), mean(cos(alphas))))
    }
```

```R
mean_direction(c(2*pi-0.1, 0.1))
```

```R
merged %>% head()
```

```R
phase_deviations = merged %>% 
    dplyr::filter(amplitude_velocycle > amplitude_threshold, padj < 0.05) %>% 
    dplyr::mutate(a = sin(acrophase_angular - mean_direction(acrophase_angular)), b = sin(peak_phase-mean_direction(peak_phase))) %>% 
    dplyr::select(a, b)
```

```R
phase_test = cor.test(x = phase_deviations$a, y = phase_deviations$b, method = "spearman")
phase_test$p.value
```

```R
phase_test$estimate
```

```R
phase_test$estimate[["rho"]]
```

```R
p = merged %>% dplyr::filter(amplitude_velocycle > amplitude_threshold, padj < 0.05) %>% 
    ggplot(aes(x = acrophase, y = peak_phase)) +
    #ggrastr::geom_point_rast(size = 0.5) + 
    geom_point(size = 0.5) +
    theme_classic() + 
    xlab("acrophase\nbulk RNA-seq, cosinor") + 
    ylab("peak phase\nscRNA-seq, VeloCycle") +
    ggtitle(paste0("Circ. Spearman's\nrho = ", 
                   phase_test$estimate[["rho"]] %>% format(., digits = 2),
                   ", p ",
                   ifelse(phase_test$p.value <=2.2e-16,
                         yes = " < 2.2e-16",
                         no = paste0(" = ", phase_test$pvalue))))
p

ggsave("../out/cell_cycle_figures/phase_comparisons_sc_vs_bulk.svg", p, svg, height = 2.5, width = 2.5)
```

## GO term enrichment in rhythmic genes

```R
rhythm = rhythm %>% dplyr::mutate(phase_assigned = factor(phase_assigned, levels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M")))
```

```R
GOcluster_per_phase <- clusterProfiler::compareCluster(symbol~phase_assigned, 
                            data = rhythm %>% dplyr::filter(padj < 0.05) %>% dplyr::select(symbol, phase_assigned), 
                            ont = "BP", 
                            keyType = 'SYMBOL',
                            universe = rhythm %>% dplyr::select(symbol) %>% unique() %>% pull(), 
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
p = clusterProfiler::dotplot(GOcluster_per_phase, showCategory = 5) + 
scale_colour_gradientn(colours = viridis(100), limits=c(0, 0.05))
p
```

```R
ggsave('../out/cell_cycle_figures/rhythmic_genes_enrichment_per_phase_shortened.svg', p, svg, height = 10, width = 6)
```

```R
GOcluster_per_phase %>% write.table('../out/tables/TableSupp_GOBP_FigS1E.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

## MESOR, RT, RIF1, H3K9me3 per bin

```R
plot_selection = rhythm %>% 
    dplyr::filter(padj < 0.05) %>%
    dplyr::left_join(., gene_metadata)
```

```R
p = plot_selection %>% 
    dplyr::select(phase_assigned, mesor) %>% 
    dplyr::rename(`bas. expr.` = mesor) %>% 
    ggplot(aes(x = phase_assigned, y = `bas. expr.`)) + 
    geom_violin() + 
    geom_boxplot(width=0.1, outlier.size = 0.5) + 
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("bas. expr.") + 
    scale_y_continuous(limits = c(-5,10), expand = c(0, 0))
p
dev.copy(svg, '../out/cell_cycle_figures/bas_expr_per_phase.svg', height = 2.5, width = 4.5)
dev.off()

```

```R
p = plot_selection %>% 
    dplyr::select(phase_assigned, RTindex_K562) %>% 
    dplyr::rename(RT = RTindex_K562) %>% 
    ggplot(aes(x = phase_assigned, y = RT)) + 
    geom_violin() + 
#    geom_jitter(width = 0.2, height = 0, size = 0.3, alpha = 0.3) + 
    geom_boxplot(width=0.1, outlier.size = 0.5) + 
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("RT") + 
    scale_y_continuous(limits = c(-6,6), expand = c(0, 0))
p
dev.copy(svg, '../out/cell_cycle_figures/RTindex_per_phase.svg', height = 2.5, width = 4.5)
dev.off()

```

```R
# RIF1
p = rhythm %>% 
    dplyr::filter(padj < 0.05) %>%
    dplyr::left_join(., gene_metadata) %>%
    dplyr::select(phase_assigned, RIF1_tss_body) %>% 
    dplyr::rename(RIF1 = RIF1_tss_body) %>% 
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

dev.copy(svg, '../out/cell_cycle_figures/RIF1_per_phase.svg', height = 2.5, width = 4.5)
dev.off()
```

```R
# H3K9me3 coverage at promoter
p = rhythm %>% 
    dplyr::filter(padj < 0.05) %>%
    dplyr::arrange(acrophase) %>%
    dplyr::mutate(H3K9me3 = h3k9me3_signal_tss) %>%
    dplyr::left_join(., gene_metadata) %>%
    dplyr::select(phase_assigned, H3K9me3) %>% 
    ggplot(aes(x = phase_assigned, y = H3K9me3)) + 
    geom_hline(yintercept = 0) +
    geom_violin() + 
    #geom_jitter(width = 0.2, height = 0, size = 0.3, alpha = 0.3) + 
    geom_boxplot(width=0.1, outlier.size = 0.5) + 
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("H3K9me3 [frac. of prom. covered]")
p
dev.copy(svg, '../out/cell_cycle_figures/H3K9me3_signal_per_phase.svg', height = 2.5, width = 4.5)
dev.off()
```

## GO BP term enriching within TFs enriched at promoters of specific phases

```R
e = clusterProfiler::compareCluster(TF_symbol ~ promoters_active_in_phase,
                              data = in_phase_vs_all_genes_enrichment_long %>% dplyr::filter(padj_enrich < 0.05),
                              OrgDb = "org.Hs.eg.db", 
                              keyType = "SYMBOL",
                              ont = "BP",
                              universe = in_phase_vs_all_genes_enrichment_long %>% pull(TF_symbol) %>% unique(),
                            maxGSSize = 100,
                             )
```

```R
library(viridis)
p = enrichplot::dotplot(e) + 
    scale_colour_gradientn(colours = viridis(100), limits=c(0, 0.05)) + 
    xlab(NULL) +
    xlab("Acrophase of rhythmic promoters") + 
    ggtitle(paste0("Enrichments amongst ", in_phase_vs_all_genes_enrichment_long %>% pull(TF_symbol) %>% unique() %>% length(), " DNA binding\nproteins at rhythmic promoters"))

p
```

```R
ggsave('../out/cell_cycle_figures/TF_binding_in_rhythmic_promoters_all_TFs_padj.svg', p, svg, height = 9.5, width = 7)
```

```R
e %>% as.data.frame() %>% write.table('../out/tables/TableSupp_GOBP_FigS1J.tsv', col.names = T, row.names = F, sep = '\t', quote = F)
```

## TF binding enrichment significance vs specificity for promoters of rhythmic genes in each phase

```R
## x axis: -log10(padj) instead of rank
to_plot = in_phase_vs_all_genes_enrichment_long %>% 
    #dplyr::filter(!promoters_active_in_phase %in% c("eG1", "lG1", "M")) %>%
    dplyr::filter(padj_enrich < 0.05) %>%
    # for labelling KZFPs better
    dplyr::group_by(promoters_active_in_phase) %>% 
    dplyr::mutate(thresh_label_padj = sort(padj_enrich)[8],
                 thresh_label_spec = sort(specificity, decreasing = T)[8])  %>%
    dplyr::mutate(TF_symbol_label = ifelse((padj_enrich<=thresh_label_padj | specificity >= thresh_label_spec | promoters_active_in_phase == "eG1" | is_kzfp), TF_symbol, "")) %>%
    dplyr::ungroup()
```

```R
p = to_plot %>%
    ggplot(aes(x = -log10(padj_enrich), y = log2(specificity), label = TF_symbol, col = is_kzfp, size = n_promoters_bound)) + 
    #ggrastr::geom_point_rast(alpha = 0.8) + 
    geom_point(alpha = 0.8) + 
    ggrepel::geom_text_repel(aes(x = -log10(padj_enrich), y = log2(specificity), label = TF_symbol_label, col = is_kzfp),
                              min.segment.length = 0, segment.alpha = 0.8, size = 3.5, label.size = 0, fill = alpha(c("white"),0), max.overlaps = 50, box.padding = 0.3, seed = 1) + 
    xlim(c(0, 45)) +
    ylim(c(0, 9)) +
    ylab("Phase specificity [log]") + 
    xlab("-log10(adj. p)") +
    theme_classic() +
    guides(color = "none") +
    scale_size(name = "Nr. bound prom.", breaks = c(10, 100, 1000, 10000), range = c(1, 3.5)) + 
    scale_color_manual(values = c("black", "#FF681F")) + 
    facet_wrap( ~ promoters_active_in_phase, nrow = 2, dir = "h", scales = "free", drop = T)
```

```R
p
```

```R
ggsave('../out/cell_cycle_figures/specialization_vs_log10padj_all_TFs.svg', p, svg, width = 12, height = 7)
```

## Visualizing the acrophase of TF targets

```R
TF_to_targets = read.table('../out/tables/TFs_to_target_promoters.tsv', sep = '\t', header = T, quote = '')
TF_to_targets %>% head()
TF_to_targets %>% dim()
TF_to_targets$TF_symbol %>% unique() %>% length()

```

```R
# density of acrophases in rhythmic targets
plot_target_acrophases <- function(gene_of_interest, padj_threshold = 0.05, bw_adjust = 0.5, padding = T) {


    target_genes = TF_to_targets %>% dplyr::filter(TF_symbol == gene_of_interest)
    
    T = 8
    central_df = rhythm
    left_df = central_df %>% dplyr::mutate(acrophase = acrophase-T)
    right_df = central_df %>% dplyr::mutate(acrophase = acrophase+T)

    padded_df = rbind(left_df, central_df, right_df)
    padding_factor = 3
    bw_adjust_padded = bw_adjust/3
    
    if(!padding) {
        padding_factor = 1
        padded_df = central_df
    bw_adjust_padded = bw_adjust}

    p_dens_line = rbind(padded_df %>% 
        dplyr::filter(padj < padj_threshold, (ensembl %in% target_genes$promoter_ensembl | symbol %in% target_genes$promoter_symbol)) %>%
                        dplyr::mutate(`acr. dens. of` = paste0(gene_of_interest, " targets")),
                        padded_df %>% dplyr::filter(padj < padj_threshold) %>%
                        dplyr::mutate(`acr. dens. of` = "all genes")) %>%
        dplyr::mutate(`acr. dens. of` = factor(`acr. dens. of`, levels = c(paste0(gene_of_interest, " targets"), "all genes"))) %>%
                        
        ggplot(aes(x = acrophase, lty = `acr. dens. of`)) + 
        geom_density(adjust = bw_adjust_padded, key_glyph = draw_key_path) + 
        #geom_density(data = padded_df %>% dplyr::filter(padj < padj_threshold), lty = 2, adjust = bw_adjust_padded) + 
        theme_classic() + 
        theme(axis.ticks.x = element_blank(), 
              axis.text.x = element_blank(), 
              axis.line.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank(), 
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(), 
              axis.line.x = element_blank(), 
              legend.position = "right", 
              legend.direction = "horizontal",
             plot.margin = unit(c(0, 0, 0, 0), 
                                    "inches")) +
            scale_x_continuous(name = "acrophase", breaks = seq(from = 0.5, to = 7.5, by = 1), labels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M"), expand = c(0, 0)) + 
        coord_cartesian(xlim=c(0, 8)) + 
        scale_linetype_manual(values = c(1, 2))


    p_dens_col = padded_df %>% 
        dplyr::filter(padj < padj_threshold, (ensembl %in% target_genes$promoter_ensembl | symbol %in% target_genes$promoter_symbol)) %>%
        dplyr::mutate(`acr. dens.` = approxfun(density(acrophase, adjust = bw_adjust_padded))(acrophase)*padding_factor)  %>%
        ggplot(aes(x = acrophase, y = 0)) + 
        geom_point(aes(color = `acr. dens.`), shape = 108, size = 45) + 
        scale_colour_viridis_c() +
        theme_classic() + 
        theme(axis.text.x = element_text(size = 20), 
              axis.line.y = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank(), 
              axis.title.y = element_blank(), 
              axis.title.x = element_blank(), 
              axis.line.x = element_blank(), 
              legend.position = "right", 
              legend.direction = "horizontal",
             plot.margin = unit(c(0, 0, 0, 0), 
                                    "inches")) +
            scale_x_continuous(name = "acrophase", breaks = seq(from = 0.5, to = 7.5, by = 1), labels = c("eG1", "lG1", "G1/S", "S1", "S2", "G2", "G2/M", "M"), expand = c(0, 0)) + 
        scale_linetype_manual(c(1, 2)) +
        coord_cartesian(xlim=c(0, 8), expand=0)


    # combining plots
    pcol = cowplot::plot_grid(p_dens_line + theme(legend.position = "none"), 
                              p_dens_col + theme(legend.position = "none"), 
                              ncol = 1, nrow = 2, rel_heights = c(0.8, 0.2))
    legend1 <- get_legend(p_dens_col)
    legend2 <- get_legend(p_dens_line)
    legends = plot_grid(legend1, legend2, nrow = 1, ncol = 2, rel_widths = c(0.4, 0.6))
    return(plot_grid(pcol, legends, ncol = 1, nrow = 2, rel_heights = c(0.8, 0.2)))
    #return(legend2)
}
```

```R
p = plot_target_acrophases("E2F3")
p
dev.copy(svg, '../out/cell_cycle_figures/E2F3_target_acrophase_density.svg', height = 3, width = 5)
dev.off()
```

```R
p = plot_target_acrophases("MYC")
p
dev.copy(svg, '../out/cell_cycle_figures/MYC_target_acrophase_density.svg', height = 3, width = 5)
dev.off()
```

```R
p = plot_target_acrophases("ZNF519")
p
dev.copy(svg, '../out/cell_cycle_figures/ZNF519_target_acrophase_density.svg', height = 3, width = 5)
dev.off()
```

```R
p = plot_target_acrophases("ZNF786")
p
dev.copy(svg, '../out/cell_cycle_figures/ZNF786_target_acrophase_density.svg', height = 3, width = 5)
dev.off()
```

```R
p = plot_target_acrophases("MYBL2")
p
dev.copy(svg, '../out/cell_cycle_figures/MYBL2_target_acrophase_density.svg', height = 3, width = 5)
dev.off()
```

```R
p = plot_target_acrophases("E2F4")
p
dev.copy(svg, '../out/cell_cycle_figures/E2F4_target_acrophase_density.svg', height = 3, width = 5)
dev.off()
```

## Figure 2
### The UMAP positions perturbations leading to similar imbalances together

```R
# loading imbalance data
res_stats =  read.table('../out/tables/K562_perturbseq_imbalances_statistics.tsv', quote = '', header = T, sep = '\t')
res_stats %>% head()
res_stats %>% dim()
res_stats$condition %>% unique() %>% length()
```

```R
# loading UMAP coordinates
```

```R
umap_plot = read.table('../out/tables/imbalances_umap_plot.tsv', sep = '\t', header = T)
umap_plot %>% head()
umap_plot %>% dim()
```

```R
umap_plot = umap_plot %>% dplyr::left_join(., y = res_stats %>% 
                                                    dplyr::select(bins, condition, p_smaller, p_greater, tf_category) %>% 
                                                    dplyr::arrange(condition, bins) %>% 
                                                    dplyr::group_by(condition) %>% 
                                                    dplyr::summarize(p_smaller_list = list(p_smaller),
                                                    p_greater_list = list(p_greater),
                                                                    p_min = min(c(p_smaller, p_greater)),
                                                                    tf_category = unique(tf_category)) %>% 
                                           dplyr::rename(symbol = condition)) %>% 
                            dplyr::filter(!symbol == "non-targeting")
```

```R
breaks = c(0,
            pi/4+pi/8, # end of M/G1
           pi-pi/8, # end of G1
           pi+pi/6, # end of G1/S
           3*pi/2+pi/8, # end of S
           2*pi # end of G2M
          )
b = length(breaks)-1
b
```

```R
sc <- scale_colour_gradientn(colours = rev(viridis::viridis(100)), limits=range(0, 2), oob = scales::squish)

# list of plots
plot_list = list(p_smaller_list = list(), p_greater_list = list())
```

```R
# we're arranging the plots in two rows and four columns. Rows are greater or smaller. Columns are phases, from M/eG1 to G2/M.
# We only add x label to the first row of plots with the phases, and y labels to the two first plots of each column.


direction_labels = list("p_greater_list"="accumulation", "p_smaller_list"="attrition")
phase_assignment = c("M/eG1", "lG1", "G1/S", "S", "G2/M")

```

```R
for (direction in names(plot_list)) {
    for (i in 1:(length(breaks)-1)) {
        to_plot = umap_plot %>% dplyr::filter(symbol != "non-targeting") %>%
        dplyr::select(UMAP1, UMAP2) %>% 
        dplyr::mutate(pval = (umap_plot %>% 
                      apply(., 1, function(x) x[direction][[1]][[i]], simplify = T) %>% unlist())) %>% 
                            dplyr::mutate(signif = -log10(pval))
        p = to_plot %>% dplyr::rename(`-log10(p-val)`=signif) %>% 
                            ggplot(aes(x = UMAP1, y = UMAP2, color = `-log10(p-val)`)) + 
                            #ggrastr::geom_point_rast(size = 0.2, raster.dpi = 250, alpha = 0.8) + 
                            geom_point(size = 0.2, alpha = 0.8) + 
                            sc + 
                            theme_classic()
        if (direction == "p_greater_list") {
            if (i == 1) {
                # leftmost plot: x and y axis titles
                p = p + xlab(phase_assignment[i]) + ylab(direction_labels[direction]) + 
                theme(axis.line=element_blank(),
                                  axis.text.x=element_blank(),
                                  axis.text.y=element_blank(),
                                  axis.ticks=element_blank(),
                                  #axis.title.x=element_blank(),
                                  #axis.title.y=element_blank()
                             )
                } else {
                    # second to fourth plots on the first row: only x axis title
                    p = p + xlab(phase_assignment[i]) + 
                        theme(axis.line=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks=element_blank(),
                                      #axis.title.x=element_blank(),
                                      axis.title.y=element_blank()
                                 )
                }
            } else if (direction == "p_smaller_list") {
            if (i == 1) {
                # leftmost plot on the second row: only y axis label
                 p = p + ylab(direction_labels[direction]) + 
                theme(axis.line=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks=element_blank(),
                                      axis.title.x=element_blank(),
                                      #axis.title.y=element_blank()
                              )
                } else{
                p = p + theme(axis.line=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks=element_blank(),
                                      axis.title.x=element_blank(),
                                      axis.title.y=element_blank()
                              )
                
            }
            }
                            
        plot_list[[direction]][[i]] = p
    }
}
```

```R
prow = cowplot::plot_grid(
       plot_list[["p_greater_list"]][[1]] + theme(legend.position="none"),
    plot_list[["p_greater_list"]][[2]] + theme(legend.position="none"),
    plot_list[["p_greater_list"]][[3]] + theme(legend.position="none"),
    plot_list[["p_greater_list"]][[4]] + theme(legend.position="none"),
    plot_list[["p_greater_list"]][[5]] + theme(legend.position="none"),
    
  plot_list[["p_smaller_list"]][[1]] + theme(legend.position="none"),
  plot_list[["p_smaller_list"]][[2]] + theme(legend.position="none"),
  plot_list[["p_smaller_list"]][[3]] + theme(legend.position="none"),
  plot_list[["p_smaller_list"]][[4]] + theme(legend.position="none"),
  plot_list[["p_smaller_list"]][[5]] + theme(legend.position="none"),
    

  align = 'vh',
  hjust = -1,
  nrow = 2
)

legend <- get_legend(
  # create some space to the left of the legend
  plot_list[["p_greater_list"]][[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
)

ggsave('../out/cell_cycle_figures/umap_cell_cycle_dev_signif.svg', plot_grid(prow, legend, rel_widths = c(3, .4)), device = svg, width = 12, height = 5)
```

### Age of imbalance inducing genes

```R
extreme_imbalances_in_text = c("CCND3", "CCNF", "CCNK", "TFDP1",
                             "MCM2", "MCM3", "RFC4", "GINS4", "WDHD1",
                              "ATRIP", "RAD21", "RINT1", "RNF138", "TTI1",
                                "MAD2L1", "SMC1A")
```

```R
to_plot = res_stats %>% 
    dplyr::left_join(., gene_metadata %>% 
                    dplyr::mutate(age_category = cut(age_combined, breaks = c(-Inf, 50, 250, 1000, Inf))) %>%
                     dplyr::select(symbol, age_combined, age_category, isKZFP) %>% 
                     unique(),
                     by = c("condition" = "symbol")) %>%
    group_by(condition) %>% slice_min(order_by = p_min, n = 1, with_ties = F) %>%
    dplyr::filter(p_min < 0.01, age_category != "NA")

to_plot %>%
    ggplot(aes(x = age_category, y = -log10(p_min))) + 
    geom_violin() + 
    geom_point(aes(col = isKZFP), alpha = 0.8) + 
    scale_color_manual(values = c("black", "#FF681F")) + 
    geom_boxplot(width = 0.1, outliers = F) +
    ylim(c(2, 7)) +
    geom_abline(intercept = 5, slope = 0, lty = 2) +
    ggpubr::stat_compare_means(comparisons = list(c(1, 4), c(2, 4), c(3, 4)), label.y = 6, step.increase = 0.1) + 
    theme_classic() + 
    ggrepel::geom_text_repel(aes(x = age_category, y = -log10(p_min), label = condition, col = isKZFP),
                              data = rbind(to_plot %>% group_by(age_category) %>% slice_min(p_min, n = 8, with_ties = F) %>% dplyr::filter(age_combined <=250),
                                           to_plot %>% dplyr::filter(condition %in% extreme_imbalances_in_text)),
                              label.size = 0, fill = alpha(c("white"),0), max.overlaps = 100, seed = 1, force = 5) + 
    labs(x = "Emergence", y = "Imbalance\n[-log10(p. min.)]", col = "KZFP") +
    scale_x_discrete(labels=c("0-50MYA", "51-250MYA", "251-1000MYA", "\u2265 1001MYA"))
dev.copy(svg, '../out/cell_cycle_figures/imbalance_vs_age.svg', width = 6, height = 5)
dev.off()
```

## Proportion of non TFs, TFs and KZFPs leading to cell cycle imbalances

```R
p = res_stats %>% dplyr::select(condition, p_min, tf_category) %>% 
    dplyr::group_by(condition, tf_category) %>% 
    dplyr::mutate(p_min = min(p_min), tf_category = factor(tf_category, levels = c("KZFP", "ZFP", "other TF", "not a TF"))) %>% 
    unique() %>%
    ggplot(aes(x = tf_category, fill = p_min < 0.05)) + geom_bar(position = "fill") + 
    scale_y_continuous(labels = scales::percent) + 
    scale_fill_manual(values = c("grey90", "black")) +
    theme_classic() + 
    theme(axis.line = element_blank(),
         axis.title = element_blank())
p
```

```R
ggsave('../out/cell_cycle_figures/tf_category_signif_prop.svg', plot = p, device = svg, width = 4, height = 1.5)
```

### Location of KZFPs and other TFs on the perturbation UMAP

```R
umap_plot$tf_category %>% table()
```

```R
umap_plot = umap_plot %>% 
    dplyr::mutate(tf_category_signif = ifelse(umap_plot$p_min < 0.05, yes = umap_plot$tf_category, no = "n.s.")) %>%
    dplyr::mutate(tf_category_signif = ifelse(tf_category_signif == "ZFP", yes = "other TF", no = tf_category_signif)) %>%
    dplyr::mutate(tf_category_signif = factor(tf_category_signif, levels = c("KZFP", "other TF", "not a TF", "n.s.")))
```

```R
umap_plot$tf_category_signif %>% levels()
```

```R
tf_category_signif_colors = c("#FF681F", "black", "grey70", "grey60")
names(tf_category_signif_colors) = umap_plot$tf_category_signif %>% levels()
```

```R
p = umap_plot %>% 
    dplyr::filter(p_min < 0.05) %>%
    dplyr::mutate(tf_category_signif = factor(tf_category_signif, levels = c("KZFP", "other TF", "not a TF"))) %>%
    ggplot(aes(x = UMAP1, y = UMAP2, color = (tf_category_signif))) + 
    #ggrastr::geom_point_rast(alpha = 0.8) + 
    geom_point(alpha = 0.8) + 
    theme_classic() + 
    theme(axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_blank(),
           legend.title=element_blank()) + 
    scale_color_discrete(type = tf_category_signif_colors) + 
    guides(colour = guide_legend(override.aes = list(size=1)))
p
```

```R
ggsave('../out/cell_cycle_figures/UMAP_cell_positions_KZFPs.svg', plot = p, device = svg, width = 4, height = 3.5)
```

## KZFP categories: Romain

```R
all_DNA_binding_proteins = in_phase_vs_all_genes_enrichment_long$TF_symbol %>% unique()
all_DNA_binding_proteins %>% head()
all_DNA_binding_proteins %>% length()
```

```R
rhythmic_DNA_binding_proteins = rhythm %>% dplyr::filter(padj < 0.05, symbol %in% all_DNA_binding_proteins) %>% dplyr::pull(symbol) %>% unique()
rhythmic_DNA_binding_proteins %>% head()
rhythmic_DNA_binding_proteins %>% length()
```

```R
DNA_binding_proteins_enriched_at_ryhthmic_genes = in_phase_vs_all_genes_enrichment_long %>% dplyr::filter(padj_enrich < 0.05) %>% dplyr::pull(TF_symbol) %>% unique()
DNA_binding_proteins_enriched_at_ryhthmic_genes %>% head()
DNA_binding_proteins_enriched_at_ryhthmic_genes %>% length()
```

```R
res_stats %>% colnames()
```

```R
DNA_binding_proteins_inducing_imbalances = res_stats %>% dplyr::filter(p_min < 0.05, condition %in% all_DNA_binding_proteins) %>% dplyr::pull(condition) %>% unique()
DNA_binding_proteins_inducing_imbalances %>% head()
DNA_binding_proteins_inducing_imbalances %>% length()
```

```R
comb_mat = make_comb_mat(list(rhythmic = rhythmic_DNA_binding_proteins, 
                              `enriched at rhythmic promoters` = DNA_binding_proteins_enriched_at_ryhthmic_genes, 
                              `silencing induces imbalances` = DNA_binding_proteins_inducing_imbalances), 
                         mode = "distinct",
                         universal_set = all_DNA_binding_proteins)
ComplexHeatmap::UpSet(comb_mat)
```

## Cell cycle distribution of K562 transduced with non-targeting gRNAs

```R
phi_df = read.csv('../data/velocycle_output/K562_conditions_phases_phis_results_13032024.csv')
cells_per_cond = phi_df %>% dplyr::select(condition) %>% dplyr::group_by(condition) %>% dplyr::mutate(ncells = n()) %>% unique() %>% as.data.frame()
rownames(cells_per_cond) = cells_per_cond$condition
```

```R
to_plot = phi_df %>% dplyr::filter(condition == "non-targeting") %>% dplyr::mutate(bins = cut(phi, breaks = breaks, include.lowest = T, right = F))
to_plot %>% head()
```

```R
# legend
phase_assignment = c("M/eG1", "lG1", "G1/S", "S", "G2/M")
lgd = to_plot %>% dplyr::select(phi, bins) %>% dplyr::group_by(bins) %>% dplyr::summarize(phi_min = min(phi)) %>% dplyr::mutate(`phase assigned` = phase_assignment)
lgd
```

```R
p = to_plot %>% ggplot(aes(x = phi)) + 
        geom_histogram(aes(y = after_stat(count / sum(count))), bins = 100) + 
        scale_y_continuous(labels = scales::percent) + 
        ggtitle(label = paste0("non-targeting",', ', cells_per_cond["non-targeting", "ncells"], " cells")) + 
        ylab("proportion of cells") + xlab("phase (Velocycle), [0, 2pi]") +
        geom_vline(xintercept = breaks[2:5]) + 
        theme_classic() + 
        geom_text(aes(x = phi_min, y = 0.022, label = `phase assigned`), data = lgd, nudge_x = 0.46)
        # scale_fill_manual(breaks = breaks, pal = plasma)   
ggsave('../out/cell_cycle_figures/perturbseq_non_targeting_distrib_breaks.svg', p, svg, width = 3.5, height = 2.5)
p
```

## Plotting the top 10 genes in term of significant imbalance on a UMAP that also has -log10(p) on it


```R
umap_plot = umap_plot %>% dplyr::mutate(signif_sum_p_smaller = map(p_smaller_list, ~sum(-log10(.x)))) %>%
    dplyr::mutate(signif_sum_p_greater = map(p_greater_list, ~sum(-log10(.x))))

umap_plot$signif_sum_p_smaller = unlist(umap_plot$signif_sum_p_smaller)
umap_plot$signif_sum_p_greater = unlist(umap_plot$signif_sum_p_greater)

umap_plot$signif_sum_all = umap_plot$signif_sum_p_smaller + umap_plot$signif_sum_p_greater
umap_plot %>% head()
```

```R
umap_plot %>% head()
umap_plot %>% dim()
```

```R
res_stats %>% head()
res_stats %>% dim()
```

```R
# finding top 3 genes in each bin, up
top_n_genes = 3
top_genes_df = rbind(res_stats %>% dplyr::left_join(., umap_plot %>% dplyr::select(symbol, signif_sum_p_smaller, signif_sum_p_greater, signif_sum_all), by = c("condition" = "symbol")) %>% 
                     dplyr::group_by(bins) %>% 
                         arrange(p_smaller, desc(signif_sum_all)) %>% 
                         slice_head(n=top_n_genes),
                          res_stats %>% dplyr::left_join(., umap_plot %>% dplyr::select(symbol, signif_sum_p_smaller, signif_sum_p_greater, signif_sum_all), by = c("condition" = "symbol")) %>% 
                     dplyr::group_by(bins) %>% arrange(p_greater, desc(signif_sum_all)) %>% slice_head(n=top_n_genes))
```

```R
set.seed(4)
sc <- scale_colour_gradientn(colours = rev(viridis::viridis(100, begin = 0.35)), limits=range(-log10(0.05), -log10(0.001)), oob = scales::squish, na.value = "grey90")
expand_by = 15
p = umap_plot %>% dplyr::mutate(p_min = ifelse(umap_plot$p_min > 0.05, yes = NA, no = umap_plot$p_min)) %>% ggplot(aes(x = UMAP1, y = UMAP2, color = -log10(p_min))) + 
    #ggrastr::geom_point_rast(size = 0.1, alpha = 0.8) + 
    geom_point(size = 0.1, alpha = 0.8) + 
    ggrepel::geom_text_repel(data = umap_plot %>% arrange(p_min) %>% dplyr::filter(symbol %in% top_genes_df$condition), 
                             aes(x = UMAP1, y = UMAP2, label = symbol), 
                             xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
                             color = "black",
                            min.segment.length = unit(0, 'lines'),
                             segment.color = "black",
                             segment.linetype = 1,
                            force = 35,
                            #force_pull = 1,
                            max.overlaps = Inf) + 
                            theme_classic() + theme(axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_blank(),
                           legend.position = "right") + sc + expand_limits(x = c(-expand_by, expand_by),
                y = c(-expand_by, expand_by))
```

```R
p
```

```R
ggsave('../out/cell_cycle_figures/umap_cell_cycle_dev_signif_top_genes.svg', p, device = svg, width = 5, height = 4)
```

### Plotting genes belonging to the same cascade that fall together on the  UMAP

```R
umap_plot = umap_plot %>% dplyr::mutate(signif_sum_p_smaller = map(p_smaller_list, ~sum(-log10(.x)))) %>%
    dplyr::mutate(signif_sum_p_greater = map(p_greater_list, ~sum(-log10(.x))))

umap_plot$signif_sum_p_smaller = unlist(umap_plot$signif_sum_p_smaller)
umap_plot$signif_sum_p_greater = unlist(umap_plot$signif_sum_p_greater)

umap_plot$signif_sum_all = umap_plot$signif_sum_p_smaller + umap_plot$signif_sum_p_greater
umap_plot %>% head()
```

```R
umap_plot %>% head()
umap_plot %>% dim()
```

```R
res_stats %>% head()
res_stats %>% dim()
```

```R
# finding top 3 genes in each bin, up
top_n_genes = 3
top_genes_df = rbind(res_stats %>% dplyr::left_join(., umap_plot %>% dplyr::select(symbol, signif_sum_p_smaller, signif_sum_p_greater, signif_sum_all), by = c("condition" = "symbol")) %>% 
                     dplyr::group_by(bins) %>% 
                         arrange(p_smaller, desc(signif_sum_all)) %>% 
                         slice_head(n=top_n_genes),
                          res_stats %>% dplyr::left_join(., umap_plot %>% dplyr::select(symbol, signif_sum_p_smaller, signif_sum_p_greater, signif_sum_all), by = c("condition" = "symbol")) %>% 
                     dplyr::group_by(bins) %>% arrange(p_greater, desc(signif_sum_all)) %>% slice_head(n=top_n_genes))
```

```R
to_label = list(RFC = c("RFC3", "RFC4"),
               MED = c("MED12", "MED19", "MED30"),
               ATR = c("ATR", "ATRIP"),
               `repl. init.` = c("GINS1", "GINS2", "GINS3", "GINS4",
            "MCM2", "MCM3", "MCM6",
             "ORC1", "ORC4",
            "CCNE2",
             "DBF4", "DBF4B"))
```

```R
label_genes_df = umap_plot %>% dplyr::filter(symbol %in% (to_label %>% unlist() %>% as.vector()))
label_genes_df
```

```R
to_label
```

```R
foo <- function(x) {
    for (n in names(to_label)) {
        if (x %in% to_label[[n]]) {
            return(n)
            }
        }
    }
```

```R

```

```R
label_genes_df$cascade = factor(sapply(label_genes_df$symbol, foo))
label_genes_df %>% head()
```

```R
set.seed(3)
sc <- scale_colour_gradientn(colours = rev(viridis::viridis(100, begin = 0.35)), limits=range(-log10(0.05), -log10(0.001)), oob = scales::squish, na.value = "grey90")
expand_by = 15
p = umap_plot %>% dplyr::mutate(p_min = ifelse(umap_plot$p_min > 0.05, yes = NA, no = umap_plot$p_min)) %>% ggplot(aes(x = UMAP1, y = UMAP2, color = -log10(p_min))) + 
    #ggrastr::geom_point_rast(size = 0.1, alpha = 0.8) + 
    geom_point(size = 0.1, alpha = 0.8) +
    sc +
    ggnewscale::new_scale_colour() + 
    ggrepel::geom_text_repel(data = label_genes_df, 
                             aes(x = UMAP1, y = UMAP2, label = symbol, col = cascade), 
                             xlim = c(-Inf, Inf), ylim = c(-Inf, Inf),
                            min.segment.length = unit(0, 'lines'),
                             segment.linetype = 1,
                            force = 35,
                            #force_pull = 1,
                            max.overlaps = Inf) + 
                            theme_classic() + theme(axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_blank(),
                           legend.position = "right") + 
    scale_color_manual(values = c("#b49711", "#922d0c", "#cb15c9", "#1136b4")) + 
    expand_limits(x = c(-expand_by, expand_by),
                y = c(-expand_by, expand_by))
```

```R
p
```

```R
ggsave('../out/cell_cycle_figures/umap_cell_cycle_cascades.svg', p, device = svg, width = 5, height = 4)
```

## GO BP enrichment amongst imbalance-inducing targets

```R
res_stats$phase = factor(res_stats$phase, levels = phase_assignment)
```

```R
for_go = res_stats %>% dplyr::mutate(across(phase, factor))
```

```R
GOcluster_all_imbalances <- clusterProfiler::compareCluster(geneClusters = list(`all imbalances` = res_stats %>% dplyr::filter(p_min < 0.05) %>% dplyr::select(condition) %>% unique() %>% pull()), 
                                             ont = "BP", 
                            keyType = 'SYMBOL',
                            universe = for_go %>% dplyr::select(condition) %>% unique() %>% pull() %>% as.character(), 
                            OrgDb = org.Hs.eg.db, 
                            fun="enrichGO",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            minGSSize = 10,
                            maxGSSize = 500
)
```

```R
GOcluster_all_imbalances %>% 
    as.data.frame() %>% head()
```

```R
p = GOcluster_all_imbalances %>% 
    as.data.frame() %>% 
    arrange(p.adjust) %>% 
    dplyr::mutate(x_coord = 1:nrow(GOcluster_all_imbalances %>% as.data.frame())) %>%
    ggplot(aes(x = x_coord, y = -log10(p.adjust), label = Description)) + 
    geom_point(size = 2) + 
    xlab("Rank") +
    ggrepel::geom_text_repel(aes(x = x_coord, y = -log10(p.adjust), label = Description),
                             force = 10) +
    theme_classic()
p
```

```R
ggsave('../out/cell_cycle_figures/go_perturbseq_all_genes_leading_to_imbalance.svg', p, svg, width = 3, height = 3)
```

```R
GOcluster_all_imbalances %>% write.table('../out/tables/TableSupp_GOBP_FigS2C.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

## Top imbalance inducing genes per phase: TODO 


## GO BP of CRISPRi targets leading to attritions

```R
GOcluster_per_phase_depletion <- clusterProfiler::compareCluster(condition~phase, 
                            data = res_stats %>%
                                                                 dplyr::filter(p_smaller < 0.05) %>% 
                                                                 dplyr::select(condition, phase) %>%
                                                                 dplyr::mutate(across(condition, as.character)), 
                            ont = "BP", 
                            keyType = 'SYMBOL',
                            universe = for_go %>% dplyr::select(condition) %>% unique() %>% pull() %>% as.character(), 
                            OrgDb = org.Hs.eg.db, 
                            fun="enrichGO",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.1,
                            minGSSize = 10,
                            maxGSSize = 500
)
```

```R
p = clusterProfiler::dotplot(GOcluster_per_phase_depletion) + 
    scale_colour_gradientn(colours = viridis(100), limits=c(0, 0.1)) + xlab(NULL)
p
```

```R
ggsave('../out/cell_cycle_figures/go_perturbseq_genes_leading_to_phase_depletion.svg', p, svg, width = 6, height = 8)
```

```R
GOcluster_per_phase_depletion %>% write.table('../out/tables/TableSupp_GOBP_FigS2F.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

## GO BP of CRISPRi targets leading to accumulations

```R
GOcluster_per_phase_enrichment <- clusterProfiler::compareCluster(condition~phase, 
                            data = res_stats %>% dplyr::filter(p_greater < 0.05) %>% dplyr::select(condition, phase) %>% dplyr::mutate(across(condition, as.character)), 
                            ont = "BP", 
                            keyType = 'SYMBOL',
                            universe = for_go %>% dplyr::select(condition) %>% unique() %>% pull() %>% as.character(), 
                            OrgDb = org.Hs.eg.db, 
                            fun="enrichGO",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.1,
                            minGSSize = 10,
                            maxGSSize = 500
)
```

```R
p = clusterProfiler::dotplot(GOcluster_per_phase_enrichment) + scale_colour_gradientn(colours = viridis(100), limits=c(0, 0.1)) + xlab(NULL)
p
```

```R
ggsave('../out/cell_cycle_figures/go_perturbseq_genes_leading_to_phase_accumulation.svg', p, svg, width = 5, height = 5)
```

```R
GOcluster_per_phase_enrichment %>% write.table('../out/tables/TableSupp_GOBP_FigS2G.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

## GO term enrichment of DNA binding proteins falling in all three cats

```R
category_df = data.frame(
           rhythmic = all_DNA_binding_proteins %in% rhythmic_DNA_binding_proteins,
    enriched_at_rhythmic_genes = all_DNA_binding_proteins %in% DNA_binding_proteins_enriched_at_ryhthmic_genes,
           inducing_imbalances = all_DNA_binding_proteins %in% DNA_binding_proteins_inducing_imbalances) %>% dplyr::mutate(in_cats = rhythmic + enriched_at_rhythmic_genes + inducing_imbalances,
                                                                                                                          TF_symbol = all_DNA_binding_proteins)
category_df %>% head()
```

```R
category_df %>% dplyr::filter(TF_symbol == "ZNF783")
```

```R
category_df %>% dplyr::filter(in_cats == 2) %>% pull(TF_symbol) %>% length()
```

```R
in_all_three_categories = intersect(rhythmic_DNA_binding_proteins, DNA_binding_proteins_inducing_imbalances) %>% intersect(., DNA_binding_proteins_enriched_at_ryhthmic_genes)
```

```R
in_all_three_categories %>% head()
in_all_three_categories %>% length()
in_all_three_categories %>% sort() %>% print()
```

```R
gene_metadata %>% dplyr::filter(symbol %in% in_all_three_categories, isKZFP) %>% dplyr::select(symbol, Tycko_score_avg)
```

```R
e = clusterProfiler::enrichGO(category_df %>% dplyr::filter(in_cats == 3) %>% pull(TF_symbol), 
                              OrgDb = "org.Hs.eg.db", 
                              keyType = "SYMBOL",
                              ont = "BP",
                              universe = category_df %>% pull(TF_symbol),
                             maxGSSize = 100
                             )
```

```R
p = e %>% 
    as.data.frame() %>% 
    arrange(p.adjust) %>% 
    dplyr::mutate(x_coord = 1:nrow(e %>% as.data.frame())) %>%
    ggplot(aes(x = x_coord, y = -log10(p.adjust), label = Description)) + 
    geom_point(size = 2) + 
    xlab("Rank") +
    ggrepel::geom_text_repel(aes(x = x_coord, y = -log10(p.adjust), label = Description),
                             force = 10) +
    theme_classic()
p
```

```R
e %>% arrange(p.adjust) %>% head()
```

```R
e %>% write.table('../out/tables/TableSupp_GOBP_FigS2H.tsv', sep = '\t', quote = F, col.names = T, row.names = F)
```

# Fig. 3 and S3: Romain

## ZNF519 expression is rhythmic

```R
plot_rhythmic_expression("ZNF519")
dev.copy(svg, filename = '../out/cell_cycle_figures/ZNF519_rhythmic_expr.svg', width = 4, height = 4)
dev.off()
```

## Are ZNF519 true targets over-represented in genes that go in the expected direction in the KD and OE?


All genes in the ZNF519 DE analysis

| | genes with ZNF519 peak | genes without ZNF519 peak | total |
|---|:---:|:---:|:---:|
| Genes in wrong dir. | x | m-x | m |
| Genes in consistant dir. | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
x_observed = 2697
k = 2697 + 501
m = 13257 + 2697
n = 501 + 1173

```

```R
x = 0:m # is the variable tested
probs <- dhyper(x, m, n, k, log = FALSE)


# we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
pval_one_sided = sum(probs[x<=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])

pval_two_sided

```

```R

```

# Fig. 4


## Rhythmic expression of human transcription factors

```R
rhythm %>% colnames()
```

```R
spring_pastels = c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7")
plot_selection = rhythm %>% dplyr::left_join(., gene_metadata) %>% dplyr::filter(padj < 0.05, isTF) %>% arrange(acrophase)
```

```R
plot_selection$isKZFP %>% table()
```

```R
plot_selection$isC2H2 %>% table()
```

```R
plot_selection$tf_category %>% table()
```

```R
plot_selection$tf_category_plot = ifelse(plot_selection$tf_category == "C2H2 ZF", yes = "TF", no = plot_selection$tf_category)
```

```R
plot_selection$tf_category_plot %>% table()
```

```R
# FRACTION OF RHYTHMIC TFs vs detected TFs
plot_selection %>% nrow() # rhythmic TFs
rhythm %>% dplyr::left_join(., gene_metadata) %>% dplyr::filter(isTF) %>% nrow() # detected TFs
```

```R
(plot_selection %>% nrow()) / (rhythm %>% dplyr::left_join(., gene_metadata) %>% dplyr::filter(isTF) %>% nrow())
```

```R
plot_selection %>% dplyr::filter(tf_category == "KZFP") %>% nrow()

(plot_selection %>% dplyr::filter(tf_category == "KZFP") %>% nrow()) / (rhythm %>% dplyr::left_join(., gene_metadata) %>% dplyr::filter(isKZFP) %>% nrow())
```

```R
rhythm %>% nrow()
```

```R
(rhythm %>% dplyr::filter(padj < 0.05) %>% nrow()) / (rhythm %>% nrow())
```

```R
plot_selection$tf_category_plot = factor(plot_selection$tf_category_plot, levels = c("TF", "KZFP"))
```

```R
gene_category = c("#000000", "#FF681FCD")
names(gene_category) = levels(plot_selection$tf_category_plot)
```

```R
cols = list(
            RT = colorRamp2(breaks = seq(from = 2.5, to = 4.5, length.out = 4), colors = viridis::viridis(4)),
            `bas. expr.` = colorRamp2(breaks = c(-2.5, 2.5, 4, 6, 10), colors = viridis::mako(5)),
             RIF1 =  colorRamp2(c(-0.5, -0.2, 0, 0.2, 0.5), hcl.colors(n = 5, palette = "Cividis", rev = F)),
           `cat.` = gene_category
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
row_anno = plot_selection %>% dplyr::select(RTindex_K562, RIF1_tss_body, mesor, tf_category_plot) %>% 
    dplyr::rename(RT = RTindex_K562, RIF1 = RIF1_tss_body, `bas. expr.` = mesor, `cat.` = tf_category_plot)
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
colnames(to_plot) = rhythm$phase_assigned %>% levels()
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
      col = expr_color,
      use_raster = T,
       raster_quality = 5,
       raster_by_magick = TRUE) %v% NULL

draw(htmp, heatmap_legend_side="top", annotation_legend_side="right")
dev.copy(svg, '../out/cell_cycle_figures/cycling_TFs_no_names_heatmap.svg', height = 5, width = 5)
dev.off()
```

## Proportions of TFs and KZFPs that peak in M-to-S1 vs. S2-to-G2/M



Computing the enrichment using Fisher's exact test: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/

The two gene categories are (among all significantly rhythmic genes):
- KZFP
- TF excluding KZFPs

The two acrophase categories are:
- phi in [S1-G2]
- phi in [M-G1/S]



```R
k = plot_selection %>% dplyr::filter(tf_category == "KZFP") %>% nrow() # number of KZFPs

acr = c( "S2", "G2", "G2/M")

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
plot_selection$phase_assigned %>% table()
```

```R
k = plot_selection %>% dplyr::filter(tf_category == "C2H2 ZF") %>% nrow() # number of KRAB-less C2H2

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

m = rhythm %>% dplyr::filter(padj < 0.05, phase_assigned %in% acr) %>% nrow() # number of genes in S-G2

n = rhythm %>% dplyr::filter(padj < 0.05, !phase_assigned %in% acr) %>% nrow() # number of genes in M-G1/S

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
p = plot_selection %>% dplyr::rename(`cat.`=tf_category) %>% ggplot(aes(ifelse(phase_assigned %in% acr, yes = "G2-c", no = "G1-c"), fill = `cat.`)) + 
    geom_bar() + 
    ggplot2::scale_fill_manual(values = cols$`cat.`) + 
    theme_classic() + 
    ylab("Number of rhythmic TFs (p-adj < 0.05)") + 
    xlab("acrophase") + 
    annotate(geom="text", x=2, y=n_tf_right+n_kzfp_right/2, label=paste0("KZFP p-val=", formatC(pval_kzfp_vs_TFs_two_sided, format = 'e', digits=2))) + 
    annotate(geom="text", x=2, y=n_tf_right+n_kzfp_right + n_C2H2_right/2, label=paste0("C2H2-ZF p-val=", formatC(pval_c2h2_vs_TFs_two_sided, format = 'e', digits=2)))
p
ggsave("../out/cell_cycle_figures/KFZPs_vs_TFs_enrichment.svg", p, svg, height = 4, width = 5)

```

```R
plot_selection %>% dplyr::filter(phase_assigned %in% c("S2", "G2", "G2/M"), isKZFP)
```

## Age of M-to-S1 vs S2-to-G2/M KZFPs

```R
plot_selection = rhythm %>% 
    dplyr::left_join(., gene_metadata) %>%
    dplyr::filter(!is.na(age_combined), padj < 0.05, isKZFP)
age_means= plot_selection %>%
    dplyr::mutate(in_S2_to_M = phase_assigned %in% acr) %>%
    dplyr::group_by(in_S2_to_M) %>% 
    dplyr::summarize(age_mean = mean(age_combined), n_kzfps = n())
age_means
```

```R
age_means_diff_observed = diff(age_means$age_mean)
age_means_diff_observed
```

Now, running permutations: we split the KZFPs in groups of the same size as the one observed and compute the difference in mean age. 

```R
to_permute = plot_selection %>% dplyr::pull(age_combined)
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
library(grid)
# Create a text
grob <- grobTree(textGrob(paste0(as.character(k_draws), " permutations, p = ", as.character(permute_pval)), x=0.1,  y=0.95, hjust=0, gp=gpar(fontsize=9)))
# Plot
p = plot_selection %>% 
    ggplot(aes(x = phase_assigned %in% acr, y = age_combined)) + 
    geom_violin() + 
    geom_boxplot(width = 0.1, outliers = T) +
    #geom_point(position = position_jitter(width = 0.1, height = NULL), size = 0.6) + 
    #stat_compare_means(label.x = 0.9, label.y = 330, step.increase = 1) + # wilcoxon is not correct when there are many ties
    theme_classic() + 
    theme(axis.ticks.x=element_blank(), axis.line.x = element_blank()) + 
    xlab("") + ylab("age (MY)") + 
    scale_x_discrete(breaks=c(F, T),
        labels=c("G1-c", "G2-c")) + 
    coord_cartesian(ylim = c(0, 350),
                      clip = 'off') +
    annotation_custom(grob)
p
```

```R
ggsave('../out/cell_cycle_figures/KZFPs_acrophase_vs_age.svg', p, svg, width = 3, height = 3)
```

## Fraction of rhythmic KZFPs undergoing an H3K9me3 loss on their zinc finger encoding array upon ZNF274 KD

TODO doublecheck delta ZNF274KO znf assignations, results are strange

```R
plot_selection = rhythm %>% 
    dplyr::left_join(., gene_metadata) %>%
    dplyr::filter(isKZFP, padj < 0.05, !is.na(deltaK9_ZNF274KO_znf))
```

```R
enrich_stats = plot_selection %>% 
    dplyr::arrange(cluster, acrophase) %>% 
    dplyr::select(padj, phase_assigned, cluster, deltaK9_ZNF274KO_znf)
colnames(enrich_stats)
```

On all RIF1-low KZFPs:

| | KZFPs with K9 loss | KZFPs without K9 loss | total |
|---|:---:|:---:|:---:|
| M-to-S1 KZFPs | x | m-x | m |
| S2-to-G2 KZFPs | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
#clusters_rif1high = c("chr6.1", "chr19.8", "chr19.9", "chr19.11")

```

```R
k = nrow(enrich_stats %>% dplyr::filter(deltaK9_ZNF274KO_znf))
k
```

```R
acr = c('S2','G2','G2/M')
```

```R
x_observed = nrow(plot_selection %>% dplyr::filter(deltaK9_ZNF274KO_znf, !phase_assigned %in% acr))
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
dplyr::mutate(H3K9me3 = c("H3K9me3 unaffected", "H3K9me3 loss")[ifelse(test = (deltaK9_ZNF274KO_znf),
                                     yes = T,
                                     no = F)+1],
              acr_group = c("G1-c", "G2-c")[phase_assigned %in% acr+1]) %>%
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

dev.copy(svg, "../out/cell_cycle_figures/deltaK9_enrichment_KZFPs.svg", width = 3.5, height = 2.2)
dev.off()
```

### ZNF274 binding in K562 at M-to-S1 vs S2-to-G2/M KZFPs

```R
acr = c("S2", "G2", "G2/M")
```

```R
plot_selection = rhythm %>% 
    dplyr::left_join(., gene_metadata) %>%
    dplyr::filter(isKZFP, padj < 0.05)

plot_selection %>% ggplot(aes(x = c("G1-c", "G2-c")[phase_assigned %in% acr+1], y = log(ZNF274_znfs+1))) + 
geom_violin(fill = "grey70", col = "white") + 
geom_boxplot(outlier.shape=NA, width = 0.1) + 
geom_hline(yintercept = median(log(plot_selection %>% dplyr::filter(padj < 0.05, !phase_assigned %in% acr) %>% pull(ZNF274_znfs) + 1)), lty = 2) +
theme_classic() + 
ylab("log(ZNF274 ChIP-seq)") +
xlab("") +
ggpubr::stat_compare_means(label.y = 13)
#ggtitle(paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2)))
dev.copy(svg, "../out/cell_cycle_figures/znf274_rhythmic_KZFPs_znfs_s2-g2_vs_m-g1s.svg", height = 3, width = 2)
dev.off()
```

### Fold change upon ZNF274KO in unsorted 293Ts

```R
rhythm %>% 
    dplyr::left_join(., gene_metadata) %>%
    dplyr::filter(isKZFP, padj < 0.05) %>% dplyr::pull(foldChange_ZNF274KO_vs_WT_HEK293T) %>% hist(xlim = c(-5, 5), breaks = 100000)
```

```R
plot_selection = rhythm %>% 
    dplyr::left_join(., gene_metadata) %>%
    dplyr::filter(isKZFP, padj < 0.05) %>% 
    dplyr::mutate(logFC = log(foldChange_ZNF274KO_vs_WT_HEK293T^sign(foldChange_ZNF274KO_vs_WT_HEK293T)+1))
```

```R
plot_selection$logFC %>% hist()
```

```R
plot_selection$logFC
```

```R
plot_selection %>% ggplot(aes(x = c("G1-c", "G2-c")[phase_assigned %in% acr+1], y = logFC)) + 
geom_violin(fill = "grey70", col = "white") + 
geom_boxplot(outlier.shape=NA, width = 0.1) + 
geom_hline(yintercept = median(plot_selection %>% dplyr::filter(padj < 0.05, !phase_assigned %in% acr) %>% pull(logFC)), lty = 2) +
theme_classic() + 
ylab("logFC\nZNF274KO/WT") +
xlab("") +
ggpubr::stat_compare_means()
#ggtitle(paste0("p-val=", formatC(t.274_znfs$p.value, format = "e", digits = 2)))
dev.copy(svg, "../out/cell_cycle_figures/znf274_rhythmic_KZFPs_FCZNF274KO_s2-g2_vs_m-g1s.svg", height = 3, width = 2.5)
dev.off()
```

## standard KRAB summary score:

```R
plot_selection %>% ggplot(aes(x = c("G1-c", "G2-c")[phase_assigned %in% acr+1], y = Tycko_score_avg)) + 
geom_violin(fill = "grey70", col = "white") + 
geom_boxplot(outlier.shape=NA, width = 0.1) + 
geom_hline(yintercept = median(plot_selection %>% dplyr::filter(padj < 0.05, !phase_assigned %in% acr) %>% pull(Tycko_score_avg)), lty = 2) +
theme_classic() + 
ylab("logFC\nZNF274KO/WT") +
xlab("") +
ggpubr::stat_compare_means()
dev.copy(svg, "../out/cell_cycle_figures/znf274_rhythmic_KZFPs_tycko_score_s2-g2_vs_m-g1s.svg", height = 3, width = 2.5)
dev.off()
```

### ZNF274 expression is rhythmic

```R
plot_rhythmic_expression("ZNF274")
dev.copy(svg, filename = '../out/cell_cycle_figures/ZNF274_rhythmic_expr.svg', width = 4, height = 4)
dev.off()
```

## KZFPs vs other genes in G1 and S/G2 upon ZNF274 KO

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
plot_selection_1 = rhythm %>% dplyr::left_join(de_genes$table %>% 
                                dplyr::select(-logCPM, -PValue) %>% 
                               dplyr::rename(ensembl = genes) %>%
                                                dplyr::mutate(znf274_ko_gate = "G1"))
```

```R
# DE in S/G2:
et = exactTest(y, pair=c("WT_S", "ZNF274KO_S"))
de_genes = topTags(et, n = nrow(et$table))
de_genes$table %>% head()
```

```R
plot_selection_2 = rhythm %>% dplyr::left_join(de_genes$table %>% 
                                dplyr::select(-logCPM, -PValue) %>% 
                               dplyr::rename(ensembl = genes) %>%
                                                dplyr::mutate(znf274_ko_gate = "S/G2"))
```

```R
plot_selection = rbind(plot_selection_1, plot_selection_2) %>% dplyr::left_join(gene_metadata)
```

```R
library(rstatix)
plot_selection %>% 
    dplyr::filter(!is.na(isKZFP), !is.na(znf274_ko_gate), padj < 0.05) %>%
    dplyr::group_by(znf274_ko_gate) %>%
    rstatix::wilcox_test(logFC ~ isKZFP) %>% 
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj")
```

```R
plot_selection %>% 
    dplyr::filter(!is.na(isKZFP), !is.na(znf274_ko_gate), padj < 0.05) %>%
    dplyr::group_by(isKZFP) %>%
    rstatix::wilcox_test(logFC ~ znf274_ko_gate, paired = T) %>% 
    rstatix::adjust_pvalue(method = "BH") %>%
    rstatix::add_significance("p.adj")
```

```R
plot_selection %>% 
    dplyr::filter(!is.na(isKZFP), !is.na(znf274_ko_gate), padj < 0.05) %>%
    ggplot(aes(x = znf274_ko_gate, y = logFC)) +
    geom_violin(aes(fill = isKZFP), position = position_dodge(0.9)) + 
    geom_boxplot(aes(fill = isKZFP), position = position_dodge(0.9), width = 0.1) + 
    theme_classic()
```

Shouldn't we show the difference in fold change between the G1 and S/G2 conditions instead?

```R
# DE in G1
et = exactTest(y, pair=c("WT_G1", "ZNF274KO_G1"))
de_genes = topTags(et, n = nrow(et$table))
de_genes$table %>% head()

```

```R
plot_selection = rhythm %>% dplyr::left_join(de_genes$table %>% 
                                dplyr::select(-logCPM, -PValue) %>% 
                               dplyr::rename(ensembl = genes,
                                            logFC_WT_vs_ZNF274KO_G1 = logFC,
                                            FDR_WT_vs_ZNF274KO_G1 = FDR))
```

```R
# DE in S/G2:
et = exactTest(y, pair=c("WT_S", "ZNF274KO_S"))
de_genes = topTags(et, n = nrow(et$table))
de_genes$table %>% head()
```

```R
plot_selection = plot_selection %>% dplyr::left_join(de_genes$table %>% 
                                dplyr::select(-logCPM, -PValue) %>% 
                               dplyr::rename(ensembl = genes,
                                            logFC_WT_vs_ZNF274KO_S = logFC,
                                            FDR_WT_vs_ZNF274KO_S = FDR))
```

```R
plot_selection$logFC_WT_vs_ZNF274KO_G1 %>% head()
plot_selection$logFC_WT_vs_ZNF274KO_S %>% head()
```

```R
plot_selection %>% dplyr::left_join(gene_metadata) %>%
    dplyr::filter(isTF, !is.na(isKZFP), padj < 0.05) %>%
    dplyr::mutate(tf_category_plot = factor(ifelse(isKZFP, yes = "KZFP", no = "TF"), levels = c("TF", "KZFP"))) %>%
    dplyr::group_by(isKZFP) %>%
    dplyr::mutate(deltaLogFC = logFC_WT_vs_ZNF274KO_S-logFC_WT_vs_ZNF274KO_G1) %>%
    ggplot(aes(x = tf_category_plot, y = deltaLogFC)) +
    geom_violin(fill = "grey70", col = "white") + 
    geom_boxplot(outliers = F, width = 0.1) + 
    theme_classic() + 
    xlab("") +
    ylab("\u0394 logFC, S/G2 vs. G1")
    ylim(c(-0.75, 1.2))
dev.copy(svg, '../out/cell_cycle_figures/znf274ko_s_vs_g1_deltaFC_kzfps_tfs.svg', height = 2.5, width = 2)
dev.off()
    
```

```R
# does the deltalogfc differ from zero for KZFPs and other genes?

enrich_data_for_plot = plot_selection %>% dplyr::left_join(gene_metadata) %>%
    dplyr::filter(isTF, !is.na(isKZFP), padj < 0.05) %>%
    dplyr::mutate(deltaLogFC = logFC_WT_vs_ZNF274KO_S-logFC_WT_vs_ZNF274KO_G1) %>%
    dplyr::group_by(isKZFP) %>% 
    summarise(P = wilcox.test(deltaLogFC, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", "ns"),
              MaxWidth = max(deltaLogFC))
enrich_data_for_plot
```

### G1-centered vs G2-centered rhythmic KZFPs, difference in logFC upon ZNF274 KO in G1 vs S/G2

```R
plot_selection %>% dplyr::left_join(gene_metadata) %>%
    dplyr::filter(isKZFP, padj < 0.05) %>%
    dplyr::mutate(deltaLogFC = logFC_WT_vs_ZNF274KO_S-logFC_WT_vs_ZNF274KO_G1) %>%
    ggplot(aes(x = c("G1-c", "G2-c")[phase_assigned %in% acr+1], y = deltaLogFC)) +
    geom_violin(fill = "grey70", col = "white") + 
    geom_boxplot(outliers = F, width = 0.1) + 
    theme_classic() + 
    xlab("") +
    ylab("\u0394 logFC, S/G2 vs. G1") + 
    ylim(c(-0.75, 1.2))
dev.copy(svg, '../out/cell_cycle_figures/znf274ko_s_vs_g1_deltaFC_m-to-s1_vs_s2_to_g2_kzfps.svg', height = 2.5, width = 2)
dev.off()
```

```R
enrich_data_for_plot = plot_selection %>% dplyr::left_join(gene_metadata) %>%
    dplyr::filter(isKZFP, padj < 0.05) %>%
    dplyr::mutate(deltaLogFC = logFC_WT_vs_ZNF274KO_S-logFC_WT_vs_ZNF274KO_G1) %>%
    dplyr::group_by(phase_assigned %in% acr) %>% 
    summarise(P = wilcox.test(deltaLogFC, mu = 0)$p.value,
              Sig = ifelse(P < 0.05, "*", "ns"),
              MaxWidth = max(deltaLogFC))
enrich_data_for_plot
```

```R
plot_selection %>% dplyr::left_join(gene_metadata) %>%
    dplyr::filter(isKZFP, padj < 0.05) %>%
    ggplot(aes(x = c("G1-c", "G2-c")[phase_assigned %in% acr+1], y = logFC_WT_vs_ZNF274KO_S)) +
    geom_violin(fill = "grey70", col = "white") + 
    geom_boxplot(outliers = F, width = 0.1) + 
    theme_classic() + 
    xlab("") +
    ggpubr::stat_compare_means() + 
    ylab("LogFC, S/G2")
```

```R
dev.copy(svg, '../out/cell_cycle_figures/znf274ko_s_vs_g1_deltaFC_kzfps_tfs.svg', height = 2.5, width = 2)
dev.off()
    
```

## Figure 5

### RT of KZFPs vs other TFs

```R
plot_selection = rhythm %>% 
    dplyr::left_join(gene_metadata)

plot_selection = plot_selection %>% 
    dplyr::mutate(tf_category_plot = ifelse(plot_selection$tf_category %in% c("TF", "C2H2 ZF"), yes = "TF", no = plot_selection$tf_category))
plot_selection$tf_category_plot %>% unique()
```

```R
plot_selection %>% 
    dplyr::filter(tf_category_plot %in% c("KZFP", "TF")) %>% 
    ggplot(aes(x = tf_category_plot, y = RTindex_K562)) +
    geom_violin(fill = "grey70") + 
    geom_boxplot(width = 0.1) +
    theme_classic() + 
    ggpubr::stat_compare_means()
```

### RT of rhythmic TFs and KZFPs

```R
plot_selection %>% 
    dplyr::filter(tf_category_plot %in% c("TF")) %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = c("G1-c", "G2-c")[phase_assigned %in% acr+1], y = RTindex_K562)) +
    geom_violin(fill = "grey70") + 
    geom_boxplot(width = 0.1) +
    theme_classic() + 
    ggpubr::stat_compare_means()
```

```R
plot_selection %>% 
    dplyr::filter(tf_category_plot %in% c("KZFP")) %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = c("G1-c", "G2-c")[phase_assigned %in% acr+1], y = RTindex_K562)) +
    geom_violin(fill = "grey70") + 
    geom_boxplot(width = 0.1) +
    theme_classic() + 
    ggpubr::stat_compare_means()
```

### Correlation between RT index and RIF1 binding at TFs and KZFPs

```R
plot_selection %>% 
    dplyr::filter(tf_category_plot %in% c("TF")) %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = RIF1_tss_body, y = RTindex_K562)) +
    geom_point() + 
    theme_classic() +
    xlim(c(-1.5, 1.5)) + 
    ggpubr::stat_cor(method = "spearman")
```

```R
plot_selection %>% 
    dplyr::filter(tf_category_plot %in% c("KZFP")) %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = RIF1_tss_body, y = RTindex_K562)) +
    geom_point() + 
    theme_classic() +
    xlim(c(-1.5, 1.5)) + 
    ggpubr::stat_cor(method = "spearman")
```

## Correlation between differential RT diff upon RIF1 KO and RIF1 binding at KZFPs and TFs

```R
plot_selection %>% 
    dplyr::filter(tf_category_plot %in% c("TF")) %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = RIF1_tss_body, y = RTdiff_RIF1KO)) +
    geom_point() +
    theme_classic() +
    xlim(c(-1.5, 1.5)) + 
    ggpubr::stat_cor(method = "spearman")
```

```R
plot_selection %>% 
    dplyr::filter(tf_category_plot %in% c("KZFP")) %>% 
    dplyr::filter(padj < 0.05) %>%
    ggplot(aes(x = RIF1_tss_body, y = RTdiff_RIF1KO)) +
    geom_point() +
    theme_classic() +
    xlim(c(-1.5, 1.5)) + 
    ggpubr::stat_cor(method = "spearman")
```

### Enrichment of DBPs at late unchanged regions upon RIF1KO undergoing a late-to-early switch upon H3K9me3 KD, i.e. SUV39H1/SUV39H2/SETDB1 KD

```R
LU_K9kdLtE = read.table('../out/tables/LU_K9KdLtE_vs_otherBins_enrichment_ENCODE_tronolabKZFPs.tsv', sep = '\t', header = T)
LU_K9kdLtE %>% head()
```

```R
LU_K9kdLtLr = read.table('../out/tables/LU_K9KdLtLr_vs_otherBins_enrichment_ENCODE_tronolabKZFPs.tsv', sep = '\t', header = T)
LU_K9kdLtLr %>% head()
```

```R
LU_K9kdLtE %>% 
    dplyr::left_join((gene_metadata %>% dplyr::select(symbol, isKZFP)), by = c("TF_symbol" = "symbol")) %>%
    dplyr::filter(!TF_symbol %in% (LU_K9kdLtLr %>% dplyr::filter(padj_enrich < 0.05) %>% pull(TF_symbol))) %>%
    dplyr::filter(padj_enrich < 0.05) %>%
    ggplot(aes(x = ratio_foldChange, y = -log10(padj_enrich), col = isKZFP, label = TF_symbol)) + 
    geom_point() + 
    theme_classic() +
    scale_color_manual(values = c("black", "#FF681F")) +
    ggrepel::geom_text_repel()

```

### Differential RT index upon ZNF274KO vs. differential gene expression upon ZNF274 KO for KZFPs

```R
plot_selection %>% 
    ggplot(aes(x = log2(foldChange_ZNF274KO_vs_WT_HEK293T), y = RTdiff_ZNF274KO, col = isKZFP)) +
    geom_point() +
    theme_classic()
```

## Exporting supplementary tables, aka TableSupp

```R
gene_metadata$RT %>% table()
```

```R
gene_metadata %>% 
    dplyr::select(chr, start, end, ensembl, strand, symbol, age_combined, isKZFP, isC2H2, isTF, Tycko_score_avg,
                                ZNF274_znfs, 
                                foldChange_ZNF274KO_vs_WT_HEK293T,
                               deltaK9_ZNF274KO_znf, 
                                RTindex_K562, RTdiff_RIF1KO, RTdiff_ZNF274KO) %>% 
    arrange(chr, start, end, strand) %>%
    write.table('../out/tables/TableSupp_gene_metadata.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

Column descriptions:

age_combined: estimated time of emergence of the gene, in million years.

isKZFP: TRUE iff the gene is a KZFP.

isC2H2: TRUE iff the gene is a C2H2 Zinc finger protein. C2H2 ZFPs encompass KZFPs.

isTF: TRUE iff the gene is a TF. TFs encompass both C2H2 ZFPs and KZFPs.

Tycko_score_avg: predicts the silencer power of the KRAB domain for KZFPs. Range: ]-inf; 1] where 1 means max. silencing.

ZNF274_znfs: ZNF274 ChIP-seq signal at the 3' zinc finger array-encoding exon of C2H2 zinc finger proteins in K562.

foldChange_ZNF274KO_vs_WT_HEK293T: fold change gene expression upon ZNF274 KO in HEK293T cells, as measured by RNA-seq.

deltaK9_ZNF274KO_znf: TRUE iff the 3' zinc finger array-encoding exon of C2H2 zinc finger proteins and thus KZFPs overlaps a region that has undergone H3K9me3 loss upon ZNF274 KO  in HEK293T.

RTindex_K562: replication timing index in K562. Larger (more positive) values indicate an earlier replication timing.

RTdiff_RIF1KO: differential replication timing index induced by RIF1 KO in HCT116 cells. Larger (more positive) values indicate a that RIF1 KO led to an earlier replication timing.

RTdiff_ZNF274KO: differential replication timing index induced by ZNF274 KO in HEK293T cells. Larger (more positive) values indicate a that ZNF274 KO led to an earlier replication timing.


```R
rhythm %>% head()
```

```R
rhythm %>% colnames()
```

```R
to_merge = data.frame(rhythm %>% dplyr::filter(padj < 0.05) %>% arrange(acrophase) %>% pull(ensembl), h3k9me3_signal_tss) %>% dplyr::rename(symbol = 1, H3K9me3_promoter_coverage = 2)
to_merge %>% head()
to_merge %>% dim()
```

```R
rhythm %>% dplyr::select(chr, start, end, strand, ensembl, symbol,
                        mesor, amplitude, acrophase, pvalue, padj, phase_assigned) %>% 
    dplyr::left_join(to_merge) %>% 
    write.table('../out/tables/TableSupp_cosinor.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

Column descriptions:

mesor: Midline Estimated Statistic Of Rhythm, i.e. baseline expression in voom-normalized counts, as estimated by the cosinor.

amplitude: max. estimated expression minus baseline expression, as estimated by the cosinor.

acrophase: phase of peak expression, re-scaled to the range [0-8], as estimated by the cosinor.

pvalue: of the F-test for rhythmicity derived from the cosinor.

padj: Benjamini-Hochberg adjusted p-value.

phase_assigned: factor in [eG1, lG1, G1/S, S1, S2, G2, G2/M, M] corresponding to the acrophase, rounded to the closest corresponding cell cycle bin.

H3K9me3_promoter_coverage: median-filtered H3K9me3 promoter coverage of significantly rhythmic genes arranged by acrophase (padj < 0.05). Should not be used out of context due to the median filtering.

```R
rhythm_for_plot %>% colnames()
```

```R
rhythm_for_plot %>% dplyr::select(chr, start, end, strand, ensembl, symbol,
                                  mesor, amplitude, acrophase, padj,
                                  sample_id, gate, phase_corrected, pca_outlier,
                                  norm.counts, batch, delta) %>%
write.table('../out/tables/TableSupp_rhythm_for_chronogram_plot.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

Column descriptions:

mesor: Midline Estimated Statistic Of Rhythm, i.e. baseline expression in voom-normalized counts, as estimated by the cosinor.

amplitude: max. estimated expression minus baseline expression, as estimated by the cosinor.

delta: estimated linear regression coefficient for batch adjustment.

acrophase: phase of peak expression, re-scaled to the range [0-8], as estimated by the cosinor.

padj: Benjamini-Hochberg adjusted p-value.

sample_id: bulk RNA-seq sample ID in cell cycle-sorted K562

gate: integer from 1 to 8, corresponds to the FACS gates used for sorting prior to bulk RNA-seq in cell cycle-sorted K562. 1 is for early G1, and 8 for M.

phase_corrected: final name of the gate, named after GOBP enrichment of rhythmic genes

pca_outlier: TRUE iff the correpsonding bulk RNA-seq sample is a PCA outlier, and was thus excluded from further analyses

norm.counts: voom-normalized RNA-seq counts, on the log2 scale.

batch: sequencing batch

delta: linear regression coefficient controlling the extent of batch correction

```R
to_export = expr_zscore %>% as.data.frame() %>% tibble::rownames_to_column("ensembl")
colnames(to_export) = gsub("\\.", "/", colnames(to_export))
to_export %>% head()
```

```R
to_export %>% write.table('../out/tables/TableSupp_zscored_mean_expr.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

Matrix cell content:

Values: zscored mean batch-adjusted norm counts.

```R
in_phase_vs_all_genes_enrichment_long %>% head()
```

```R
in_phase_vs_all_genes_enrichment_long %>% 
    dplyr::rename(symbol = TF_symbol) %>% 
    dplyr::select(symbol, n_promoters_bound,
                  promoters_active_in_phase, pval_enrich, padj_enrich, rank, median_rank_signif_padj, specificity) %>%
    write.table('../out/tables/TableSupp_DBP_enrichment_at_rhythmic_promoters.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

Column names:

symbol: DNA binding protein (DBP) whose binding at promoters was checked for enrichment at rhythmic promoters.

n_promoters_bound: total number of promoters bound by the DBP.

promoters_active_in_phase: phase in which the rhythmic genes whose promoters is checked for binding are peaking (i.e. acrophase of the target genes).

pval_enrich: p-value resulting from Fisher's Exact Test.

padj_enrich: Benjamini-Hochberg adjusted p-value.

rank: rank of the adj. p-value amongst DBPs enriching for binding in promoters of the phase at hand, encoded in promoters_active_in_phase.

median_rank_signif_padj: median adj. p-value rank across the 8 cell cycle phases.

specificity: median_rank_signif_padj/rank.








```R
umap_plot %>% head()
```

```R
umap_plot %>% dplyr::select(symbol, 
                            UMAP1, UMAP2,
                            signif_sum_all) %>%
    write.table('../out/tables/TableSupp_imbalance_UMAP.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
                            
```

Column names:

symbol: gene targeted for silencing by CRISPRi in the K562 perturb-seq of Replogle et al., Cell 2022

UMAP1: x coordinate of the UMAP projection based on the cell abundances relative to the unperturbed K562 distribution across cell cycle phases

UMAP2: y coordinate

signif_sum_all: sum of -log10(p) across attrition (p_smaller) and accumulation (p_greater) p-values.

```R
ls()
```

```R
res_stats %>% head()
```

```R
res_stats %>% dplyr::rename(symbol = condition, bin = bins) %>% 
    dplyr::select(symbol, bin, phase, 
                            count, count_expected_median, 
                            p_smaller, p_greater) %>%
    write.table('../out/tables/TableSupp_perturbseq_imbalances_statistics.tsv', sep = '\t', col.names = T, row.names = F, quote = F)
```

Column names:

symbol: gene targeted for silencing by CRISPRi in the K562 perturb-seq of Replogle et al., Cell 2022.

bin: boundaries of the cell cycle phase bin in which imbalance is assessed.

phase: cell cycle phase corresponding to the bin.

count: number of cells bearing the CRISPRi guides specified in symbol found in that cell cycle phase bin.

count_expected_median: median (1e5 trials) expected number of cells found in that cell cycle bin for an unperturbed population of K562 cells.

p_smaller: proportion of times count < count_expected in 1e5 trials

p_greater: proportion of times count > count_expected in 1e5 trials

```R
LU_K9kdLtE %>% 
    dplyr::left_join((gene_metadata %>% dplyr::select(symbol, isKZFP)), by = c("TF_symbol" = "symbol")) %>%
    dplyr::filter(!TF_symbol %in% (LU_K9kdLtLr %>% dplyr::filter(padj_enrich < 0.05) %>% pull(TF_symbol))) %>% 
    dplyr::filter(padj_enrich < 0.05) %>%
    dplyr::mutate(symbol = TF_symbol, DBP_peaks_total = TF_peaks_total, DBP_peaks_in_LU_K9KdLtE = TF_peaks_in_LU_K9KdLtE) %>%
    dplyr::select(symbol,
                  p_enrich, padj_enrich, 
                  ratio_observed, ratio_expected, ratio_foldChange, 
                  DBP_peaks_total, DBP_peaks_in_LU_K9KdLtE) %>% arrange(p_enrich) %>% 
    write.table('../out/tables/TableSupp_DBPs_enriched_in_LU_RIF1KO_LtE_H3K9me3KD.tsv', sep = '\t', row.names = F, col.names = T, quote = F)
```

Column names:

symbol: DNA binding protein (DBP) whose peaks were checked for enrichment in 50kB regions late unchanged upon RIF1 KO, but switching from late to early upon additional H3K9me3 KD, i.e. SUV39H1/SUV39H2/SETDB1 KD.

p_enrich: Fisher's Exact Test p-value

padj_enrich: Benjamini Hochberg-adjusted p-value

ratio_observed: fraction of peaks overlapping said regions

ratio_expected: expected fraction of peaks overlapping said regions in the case of random binding

ratio_foldChange: ratio_observed/ratio_expected

DBP_peaks_total: total number of DBP peaks

DBP_peaks_in_LU_K9KdLtE: number of DBP peaks falling within said regions
