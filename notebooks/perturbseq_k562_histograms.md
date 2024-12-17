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

# Cell cycle distribution of K562 perturb-seq

Perturb-seq single cell RNA seq data was obtained from https://doi.org/10.1016/j.cell.2022.05.013

Cell phases were estimated using VeloCycle (https://doi.org/10.1101/2024.01.18.576093 ) by Alex Lederer. 

## Goal

Computing the deviations in cell cycle distribution of each perturbation in the repli-seq.

```R
library(tidyverse)
library(boot) # BOOTSTRAP
library(cowplot) # multi-panel plots
library(clusterProfiler) # gene set enrichment analysis
library(org.Hs.eg.db)
install.packages("devEMF") # exporting plots in EMF for Romain
```

## Loading data

```R
# Velocycle resultss
phi_df = read.csv('../data/velocycle_output/K562_conditions_phases_phis_results_13032024.csv')
```

```R
phi_df %>% head()
```

```R
# total number of cells: 
phi_df %>% nrow()
```

```R
# number of cells per condition
cells_per_cond = phi_df %>% dplyr::select(condition) %>% dplyr::group_by(condition) %>% dplyr::mutate(ncells = n()) %>% unique() %>% as.data.frame()
rownames(cells_per_cond) = cells_per_cond$condition
```

```R
cells_per_cond$ncells %>% quantile(probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99))
```

```R
cells_per_cond %>% head()
```

```R
cells_per_cond %>% write.table('../out/tables/cells_per_cond_perturbseq.tsv', sep = '\t', col.names = T, row.names = F)
```

```R
# number of conditions in the perturb-seq: 
cells_per_cond %>% dim()
```

## Plotting raw cell cycle distributions for KZFPs and cell cycle marker genes, as well as genes of interest

```R
dir.create('../out/velocycle_phase_plots/non_targeting/', recursive = T)
bin_sizes = c(20, 30, 50, 100)
for (b in bin_sizes) {
    dir.create(paste0('../out/velocycle_phase_plots/kzfps/bins',b), recursive = T)
}

```

```R
conds = phi_df$condition %>% unique() %>% sort()
```

```R
# KZFP rhythmicity
fit_kzfps = read.csv('../out/tables/kzfps_cycling_info.tsv', sep = '\t', header = T)
fit_kzfps %>% head()
```

```R
kzfps_olga = read.table('../data/from_Olga/Table4Cyril.txt', sep = '\t', header = T)
kzfps_olga$summaryKZFP %>% table()
kzfps_olga %>% head(1)
```

```R
olga_kzfps_symbols = (kzfps_olga %>% dplyr::filter(Common.Name == "Human", summaryKZFP != "No.KZFP"))$KRAB.Gene.id %>% unique() %>% as.character()
olga_kzfps_symbols
```

```R
kzfps_rhythmic_and_olga = c(fit_kzfps$symbol, olga_kzfps_symbols) %>% unique() %>% sort()
kzfps_rhythmic_and_olga
kzfps_rhythmic_and_olga %>% length()
```

```R
# kzfps
conds_subset = conds[conds %in% kzfps_rhythmic_and_olga]
conds_subset %>% length()
```

```R
for (b in bin_sizes) {
    for (c in conds_subset) {
        p = phi_df %>% dplyr::filter(condition == c) %>% ggplot(aes(x = phi)) + 
            geom_histogram(aes(y = after_stat(count / sum(count))), bins = b) + 
            scale_y_continuous(labels = scales::percent) + 
            ggtitle(label = paste0(c,', ', cells_per_cond[c, "ncells"], " cells, ", b, " equally spaced bins")) + 
            ylab("frequency") + 
            geom_point(aes(x = phi, y = -0.01), size = 0.01)
        ggsave(plot = p, filename=paste0('../out/velocycle_phase_plots/kzfps/bins', b, '/', c, '.pdf'), height = 5, width = 5)
        }
}


```

```R
# non-targeting
c = "non-targeting"
for (b in bin_sizes) {
    p = phi_df %>% dplyr::filter(condition == c) %>% ggplot(aes(x = phi)) + 
        geom_histogram(aes(y = after_stat(count / sum(count))), bins = b) + 
        scale_y_continuous(labels = scales::percent) + 
        ggtitle(label = paste0(c,', ', cells_per_cond[c, "ncells"], " cells, ", b, " equally spaced bins")) + 
        ylab("frequency") + 
        geom_point(aes(x = phi, y = -0.01), size = 0.01)
    ggsave(plot = p, filename=paste0('../out/velocycle_phase_plots/non_targeting/bins', b,'.pdf'), height = 5, width = 5)
}
```

```R
# cell cycle marker genes romain
for (b in bin_sizes) {
    dir.create(paste0('../out/velocycle_phase_plots/cell_cycle_romain_gene_set/bins',b), recursive = T)
}
```

```R
cc_genes = read.csv('../data/cc_romain_genes_16_09_22.csv', header = T)$Gene %>% sort()
cc_genes %>% head()
```

```R
conds_subset = conds[conds %in% cc_genes]
```

```R
for (b in bin_sizes) {
    for (c in conds_subset) {
        p = phi_df %>% dplyr::filter(condition == c) %>% ggplot(aes(x = phi)) + 
            geom_histogram(aes(y = after_stat(count / sum(count))), bins = b) + 
            scale_y_continuous(labels = scales::percent) + 
            ggtitle(label = paste0(c,', ', cells_per_cond[c, "ncells"], " cells, ", b, " equally spaced bins")) + 
            ylab("frequency") + 
            geom_point(aes(x = phi, y = -0.01), size = 0.01)
        ggsave(plot = p, filename=paste0('../out/velocycle_phase_plots/cell_cycle_romain_gene_set/bins', b, '/', c, '.pdf'), height = 5, width = 5)
        }
}


```

## Bootstrapping from the distribution of of non-perturbed cells

Conditions sometimes have low number of cells, which leads to high variability. To get comparable statistics on cell cycle distribtutions across conditions of the perturb-seq, we compute the likelihood that the observed fraction of cells in each bin is x standard deviations away from what would be estimated by drawing the same number of cells sampled from the non-perturbed distribution.

Bins are defined as left-closed and right-open intervals, e.g. x in 0 <= x <1 -> [0, 1)


### Defining bins

We aim at having wider bins where there are less cells, and thinner bins where there are more cells. We use the non-perturbed cell distribution as a guide, looking at the phi estimates for well known cell cycle genes drawn from the unperturbed condition.

```R
breaks = c(0,
            pi/4+pi/8, # end of M/G1 = 3pi/8
           pi-pi/8, # end of G1 = 7pi/8
           pi+pi/6, # end of G1/S = 7pi/6
           3*pi/2+pi/8, # end of S = 13pi/8
           2*pi # end of G2M = 2pi
          )
b = length(breaks)-1
b
```

TODO: plotting the phi plot for genes, with seurat marker genes as was done for the clel lines. Problem: we don't have those values, we need to ask Lederer. 


### Plotting the distribution of cells transduced with non-targeting (i.e control) guide RNAs (gRNAs)

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
ggsave('../out/cell_cycle_figures/perturbseq_non_targeting_distrib_breaks.svg', p, svg, width = 4, height = 3)
p
```

How do the untreated cells look like ?

```R
p = phi_df %>% dplyr::filter(condition == "non-targeting") %>% ggplot(aes(x = phi)) + 
        geom_histogram(aes(y = after_stat(count / sum(count))), breaks = breaks) + 
        scale_y_continuous(labels = scales::percent) + 
        ggtitle(label = paste0("non-targeting",', ', cells_per_cond["non-targeting", "ncells"], " cells\n", b, " bins")) + 
        ylab("proportion of cells") + 
        geom_vline(xintercept = pi/4+pi/8)+ # end of M/G1
        geom_vline(xintercept = pi-pi/8)+ # end of G1 +
        geom_vline(xintercept = pi+pi/6) + # end of G1/S
        geom_vline(xintercept = 3*pi/2+pi/8) +  # end of S
        geom_vline(xintercept = 2*pi)  + # end of G2M
        theme_classic()

        
p
```

```R
ggsave('../out/cell_cycle_figures/perturbseq_non_targeting_distrib_binned_breaks.svg', p, svg, width = 4, height = 3)

```

Showing a few perturbed conditions

```R
g = "ZNF519"
p = phi_df %>% dplyr::filter(condition == g) %>% ggplot(aes(x = phi)) + 
        geom_histogram(aes(y = after_stat(count / sum(count))), breaks = breaks) + 
        scale_y_continuous(labels = scales::percent) + 
        ggtitle(label = paste0(g,', ', cells_per_cond[g, "ncells"], " cells\n", b, " bins")) + 
        ylab("proportion of cells") + 
        geom_vline(xintercept = pi/4+pi/8)+ # end of M/G1
        geom_vline(xintercept = pi-pi/8)+ # end of G1 +
        geom_vline(xintercept = pi+pi/6) + # end of G1/S
        geom_vline(xintercept = 3*pi/2+pi/8) +  # end of S
        geom_vline(xintercept = 2*pi)  + # end of G2M
        theme_classic()

        
p
```

```R
g = "TFDP1"
p = phi_df %>% dplyr::filter(condition == g) %>% ggplot(aes(x = phi)) + 
        geom_histogram(aes(y = after_stat(count / sum(count))), breaks = breaks) + 
        scale_y_continuous(labels = scales::percent) + 
        ggtitle(label = paste0("non-targeting",', ', cells_per_cond[g, "ncells"], " cells\n", b, " bins")) + 
        ylab("proportion of cells") + 
        geom_vline(xintercept = pi/4+pi/8)+ # end of M/G1
        geom_vline(xintercept = pi-pi/8)+ # end of G1 +
        geom_vline(xintercept = pi+pi/6) + # end of G1/S
        geom_vline(xintercept = 3*pi/2+pi/8) +  # end of S
        geom_vline(xintercept = 2*pi)  + # end of G2M
        theme_classic()
        
p
```

## Generating the bootstrapped distributions for all number of cells, by drawing from the non-targeting cells.

```R
# number of bootrap draws
k_draws = 1e3

# limits on the bootstrap draws
nCells_min = min(cells_per_cond$ncells)
nCells_max = max((cells_per_cond %>% dplyr::filter(condition!= "non-targeting") %>% dplyr::select(ncells)))
```

```R
# range of the number of cells per conditions
nCells_min
nCells_max
```

```R
# subsetting on non-targeted cells
non_targeting_df = phi_df %>% dplyr::filter(condition == "non-targeting") %>% dplyr::select(phi)
non_targeting_bins = cut(non_targeting_df$phi, breaks = breaks, include.lowest = T, right = F)
non_targeting_bins %>% head()
```

```R
# number of cells per bin for non-targeted cells
table(non_targeting_bins)
```

```R
# number of cells in the non-targeting condition
non_targeting_bins %>% length()
```

```R
count_bins <- function(data, indices, nCells) {
    # function to generate draws from the distribution of non-perturbed cells
    return(table(data[indices[1:nCells]]))
    }
```

```R
boot_test = boot(data = non_targeting_bins, statistic = count_bins, R = k_draws, nCells = 100)
```

```R
obs_counts = phi_df %>% dplyr::filter(condition == c) %>% dplyr::select(phi) %>% 
    dplyr::mutate(bins = cut(phi, breaks = breaks, include.lowest = T, right = F)) %>% 
    dplyr::group_by(bins, .drop = F) %>% dplyr::summarize(count = n())
```

```R
obs_counts$count_expected_median = boot_test$t %>% apply(., 2, median)

obs_counts$p_smaller = map(1:b, function(i) sum(boot_test$t[, i] <= obs_counts$count[i])/k_draws) # how many times in the bootstrapped data do we find a value lower than the obsered one?
obs_counts$p_greater = map(1:b, function(i) sum(boot_test$t[, i] >= obs_counts$count[i])/k_draws)# how many times in the bootstrapped data do we find a value greater than the obsered one?
```

```R
obs_counts
```

<!-- #raw -->
TODO repair

obs_counts$count_expected_median_multinorm = test %>% apply(., 2, median)

obs_counts$p_smaller_multinorm = map(1:b, function(i) sum(test[, i] <= obs_counts$count[i])/k_draws) # how many times in the bootstrapped data do we find a value lower than the obsered one?
obs_counts$p_greater_multinorm = map(1:b, function(i) sum(test[, i] >= obs_counts$count[i])/k_draws)# how many times in the bootstrapped data do we find a value greater than the obsered one?
<!-- #endraw -->

```R
obs_counts
```

The results are very similar to each other, whether we boostrap or sample from a multinomial distribution with probs defined empirically using the non-targeting cells. The second option is much cheaper computationally, so we opt for it instead. This will greatly facilitate the computations.


## Sampling from the multinomial to approximate a random draw, for all genes

For each condition, we draw N cells from the multinomial based on the non-targeting distribution, where N is the number of cells present in the condition. We repeat the process 10'000 times to get stable estimates. We then compute how frequently we find frequencies more extreme (less than, more than) than what we observe.

```R
# computing the probability distribution for the non-targeting condition
non_targeting_df = phi_df %>% dplyr::filter(condition == "non-targeting") %>% dplyr::select(phi)
non_targeting_bins = cut(non_targeting_df$phi, breaks = breaks, include.lowest = T, right = F)
non_targeting_bins %>% head()
probs = table(non_targeting_bins)/length(non_targeting_bins)
probs
```

```R
# number of cells in each cell cycle phase bin for each condition (i.e. the observed distribution)
counts_per_bin = phi_df %>% dplyr::filter(condition != "non-targeting") %>% 
    dplyr::select(condition, phi) %>% 
    dplyr::mutate(across(where(is.character), as.factor)) %>%
    dplyr::group_by(condition, .drop = F) %>% 
    dplyr::mutate(bins = cut(phi, breaks = breaks, include.lowest = T, right = F)) %>% 
    dplyr::group_by(bins, condition, .drop = F) %>% 
    dplyr::summarize(count = n()) %>% 
    dplyr::ungroup() %>%
    dplyr::arrange(condition, bins)
counts_per_bin %>% head() 
```

```R
# initializing the result df
res_df_list = list()
k_draws = 1e5
conds_subset = conds[conds != "non-targeting"]
set.seed(10) # drawing from the multinomial is random, we set a seed for reproducibility
```

```R
for (i in 1:length(conds_subset)) {
    obs_counts = counts_per_bin %>% dplyr::filter(condition == conds_subset[[i]])
    nCells = sum(obs_counts$count)
    # drawing from the multinomial k times
    multinom_draws = rmultinom(n = k_draws, 
                                size = nCells, 
                               prob = probs) %>% t()
    
    obs_counts$count_expected_median = multinom_draws %>% apply(., 2, median)
    obs_counts$p_smaller = map(1:b, function(i) sum(multinom_draws[, i] <= obs_counts$count[i])/k_draws) # how many times in the sampled data do we find a value lower than the obsered one?
    obs_counts$p_greater = map(1:b, function(i) sum(multinom_draws[, i] >= obs_counts$count[i])/k_draws)# how many times in the sampled data do we find a value greater than the obsered one?
    res_df_list[[i]] = obs_counts
}
```

```R
# gathering results
res_stats = do.call(rbind, res_df_list)
res_stats
```

```R
res_stats %>% saveRDS('../out/tables/perturbseq_cell_cycle_deviations.RDS')
```

```R
res_stats = readRDS('../out/tables/perturbseq_cell_cycle_deviations.RDS')
```

```R
# setting the minimum p-vals to 1/kdraws
res_stats$p_smaller = ifelse(res_stats$p_smaller == 0, yes = 1/k_draws, no = res_stats$p_smaller)
res_stats$p_greater = ifelse(res_stats$p_greater == 0, yes = 1/k_draws, no = res_stats$p_greater)

```

```R
res_stats$phase = phase_assignment[res_stats$bins %>% as.integer()]
```

```R
res_stats = res_stats %>% dplyr::mutate(across(where(is.list), as.double))

```

```R
res_stats$p_smaller_adj = NULL
```

```R
res_stats %>% write.table('../out/tables/perturbseq_cell_cycle_deviations.tsv', sep = '\t', col.names = T, row.names = F)
```

## Data analysis

```R
res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga, p_smaller < 0.05) %>% arrange(p_smaller) %>% head()
```

```R
res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga, p_greater < 0.05) %>% arrange(p_greater) %>% head()
```

```R
# ranking by statistical significance of the accumulation / attrition
res_stats$p_min = apply(res_stats, 1, function(x) min(x["p_smaller"], x["p_greater"])) %>% as.double()
```

```R
res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga, (p_smaller < 0.05 | p_greater < 0.05)) %>% arrange(p_min) %>% head()
```

```R
res_stats$stars = ifelse(res_stats$p_min < 0.001, yes = "***", no = ifelse(
        res_stats$p_min < 0.01, yes = "**", no = ifelse(
            res_stats$p_min < 0.05, yes = "*", no = "")))
```

Do we get some KZFPs that are not present in K562 but in Olga's table?


### Generating plots for KZFPs and cell cycle marker genes

```R
dir.create('../out/velocycle_imbalance_plots/')
```

```R
kzfps_with_imbalance = res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga) %>% dplyr::filter(p_min < 0.05) %>% pull(condition) %>% as.character() %>% unique()
```

```R
kzfps_with_imbalance %>% length()
```

```R
126/348
```

A third of KZFPs affect the cell cycle when silenced. What about other genes, TFs and KRAB-devoid C2H2 TFs?

```R
tfs = read.csv('../data/DatabaseExtract_v_1.01.csv', sep = ',', header = T, check.names = T)
tfs %>% head(1)
```

```R
# keeping only tfs assayed in the perturbseq
tfs %>% dim()

tfs_in_conds = tfs %>% dplyr::filter(Is.TF. == 'Yes', HGNC.symbol %in% (res_stats$condition %>% unique()))
tfs_in_conds %>% dim()
```

```R
res_stats %>% dplyr::filter(condition %in% tfs_in_conds$HGNC.symbol) %>% pull(condition) %>% unique() %>% length()
```

```R
res_stats %>% dplyr::filter(condition %in% tfs_in_conds$HGNC.symbol, p_min < 0.05) %>% pull(condition) %>% unique() %>% length()
```

```R
585/1519
```

KZFPs are exactly like TFs in term of enrichment

```R
tfs_in_conds$DBD %>% unique()
```

```R
c2h2 = tfs_in_conds %>% dplyr::filter(grepl('C2H2 ZF', DBD)) %>% pull(HGNC.symbol) %>% unique()
c2h2 %>% length()
```

```R
krabless_c2h2 = c2h2[!c2h2 %in% kzfps_rhythmic_and_olga]
krabless_c2h2 %>% length()
krabless_c2h2
```

```R
res_stats %>% dplyr::filter(condition %in% krabless_c2h2) %>% pull(condition) %>% unique() %>% length()
```

```R
res_stats %>% dplyr::filter(condition %in% krabless_c2h2, p_min < 0.05) %>% pull(condition) %>% unique() %>% length()
```

```R
137/356
```

Exactly the same proportion as KRABless C2H2 genes. 


What about TFs that are neither KZFPs or ZFPs?

```R
res_stats %>% dplyr::filter(condition %in% tfs_in_conds$HGNC.symbol, !condition %in% kzfps_rhythmic_and_olga, !condition %in% krabless_c2h2) %>% pull(condition) %>% unique() %>% length()
```

```R
res_stats %>% dplyr::filter(condition %in% tfs_in_conds$HGNC.symbol, !condition %in% kzfps_rhythmic_and_olga, !condition %in% krabless_c2h2, p_min < 0.05) %>% pull(condition) %>% unique() %>% length()
```

```R
322/817
```

Really exactly the same.


Annotating these TF categories:

```R
res_stats = res_stats %>% dplyr::mutate(is_kzfp = res_stats$condition %in% kzfps_rhythmic_and_olga,
                            is_c2h2 = res_stats$condition %in% krabless_c2h2,
                            is_other_tf = res_stats$condition %in% tfs_in_conds$HGNC.symbol & 
                                !res_stats$condition %in% kzfps_rhythmic_and_olga & 
                                !res_stats$condition %in% krabless_c2h2)
```

```R
res_stats$tf_category = ifelse(res_stats$is_kzfp, yes = "KZFP", no = 
                                ifelse(res_stats$is_c2h2, yes = "ZFP", no = 
                                        ifelse(res_stats$is_other_tf, yes = "other TF", no = "not a TF")))
```

```R
res_stats$tf_category %>% table()
```

```R
res_stats %>% dplyr::filter(p_min < 0.05) %>% pull(condition) %>% as.character() %>% unique() %>% length()
```

```R
res_stats$condition %>% unique() %>% length()
```

```R
4026/(res_stats$condition %>% unique() %>% length())
```

```R
res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga) %>% dplyr::filter(p_min < 0.05) %>% arrange(p_min) %>% head(20)
```

```R
# saving table for imbalance plots
res_stats %>% write.table('../out/tables/K562_perturbseq_imbalances_statistics.tsv', quote = F, row.names = F, col.names = T, sep = '\t')
```

```R
# plotting imbalances as a percentage deviation from the non-targeting cells
require(tidyverse)

res_stats_for_plot = read.table('../out/tables/K562_perturbseq_imbalances_statistics.tsv', header = T, sep = '\t')

imbalance_plot <- function(c, # condition e.g "TP53"
                           imbalance_results # res_stats
                          ) {
    p = imbalance_results %>% dplyr::filter(condition == c) %>% 
    dplyr::mutate(imbalance = (count-count_expected_median)/count_expected_median) %>% 
    dplyr::rename(p = p_min) %>%
    ggplot(aes(x = bins, y = imbalance, fill = p < 0.05)) + 
    scale_y_continuous(labels = scales::percent) + 
    geom_col() + 
    ggtitle(label = paste0(c,', ', cells_per_cond[c, "ncells"], " cells")) + ylab("Imbalance (% expected)") + 
    theme_classic() + scale_fill_manual(values = c("FALSE" = "lightgrey", "TRUE" = "black")) + 
    geom_text(aes(x = bins, y = imbalance, label = stars), col = "red") + 
    scale_x_discrete(labels = phase_assignment)
    
    return(p) ##ggplot object
    }
```

```R
for (c in c("RFC4", "MED12", "CASP8AP2", "CCNF", "ATRIP", kzfps_with_imbalance)) {
    p = imbalance_plot(c, res_stats)
    ggsave(paste0('../out/velocycle_imbalance_plots/', c, '.svg'), p, svg, width = 3.5, height = 2)
    }
    
```

## Representing the cell cycle distributions of perturbed cells on a UMAP plot

```R
# representing the conditions on a UMAP plot
library(reticulate)
library(umap)
```

```R
phi_df %>% head()
```

```R
binned_phi_df = phi_df %>% dplyr::select(condition, phi) %>% dplyr::mutate(bins = cut(phi, breaks = breaks, include.lowest = T, right = F))
binned_phi_df %>% head()
```

```R
bin_counts = binned_phi_df %>%
    dplyr::group_by(condition, bins, .drop = F) %>% dplyr::summarize(count = n()) %>% dplyr::ungroup() %>% dplyr::arrange(condition, bins)
```

```R
library(reshape2)
bin_matrix = reshape2::acast(bin_counts, bins~condition, value.var="count") %>% scale() %>% t()
```

```R
bin_matrix %>% head()
```

```R
set.seed(1)
umap_bins = umap::umap(bin_matrix)
```

```R
umap_bins %>% saveRDS('../out/umap_perturbseq.RDS')
```

```R
umap_bins$layout %>% as.data.frame() %>% ggplot(aes(x = V1, y = V2)) + geom_point()
```

```R
umap_plot = umap_bins$layout %>% as.data.frame() %>% tibble::rownames_to_column("symbol") %>% dplyr::rename(UMAP1 = V1, UMAP2 = V2)
```

```R
umap_plot$is_kzfp = umap_plot$symbol %in% kzfps_rhythmic_and_olga
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
umap_plot %>% head()
```

```R
gene_category = viridis::turbo(n = 10)[c(2, 5, 7)]
names(gene_category) = c("other TF", "ZFP", "KZFP")
```

```R
kzfp_color = gene_category[[3]]
```

```R
kzfps_highlighted = c("ZNF695", "ZNF850", "ZNF669", "ZNF519")
```

```R
library(ggrepel)
set.seed(5)
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = is_kzfp)) + geom_point(size = 0.4) + 
    ggrepel::geom_text_repel(data = umap_plot %>% dplyr::filter(symbol %in% kzfps_highlighted), 
                             aes(x = UMAP1, y = UMAP2, label = symbol), 
                             color = "black",
                            min.segment.length = unit(0, 'lines')) + 
    theme_classic() + theme(axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_blank(),
                           legend.position = "none") + scale_color_manual(values = c("grey90", "grey50")  )
```

```R
ggsave('../out/cell_cycle_figures/UMAP_cell_positions_KZFPs_highlighted.pdf', plot = p, device = pdf, width = 4, height = 4)
ggsave('../out/cell_cycle_figures/UMAP_cell_positions_KZFPs_highlighted.svg', plot = p, device = svg, width = 4, height = 4)
```

```R
umap_plot %>% head()
```

```R
# exporting data for Romain:
umap_plot %>% 
    dplyr::select(-p_smaller_list, -p_greater_list) %>% 
    write.table('../out/tables/imbalances_umap_plot.tsv', sep = '\t', col.names = T, row.names = F)
```

There are many cases where we don't find any number of cells more extreme than what we observe, i.e. there are ties. To choose top genes to show, we order by the second smallest p_val among those.

```R
umap_plot %>% arrange(p_min) %>% head(10) %>% pull(p_greater_list)
```

```R
umap_plot = umap_plot %>% 
    dplyr::mutate(p_all_list = map2(p_greater_list, p_smaller_list, ~ c(.x, .y))) %>% 
    dplyr::mutate(p_all_sorted = map(p_all_list, sort)) %>% 
    dplyr::mutate(top_2_pval = map(p_all_sorted, ~head(.x, n=2))) %>%
    dplyr::mutate(p_min_second = map(top_2_pval, ~.x[[2]])) %>% dplyr::select(-top_2_pval, -p_all_sorted, -p_all_list)
umap_plot$p_min_second = unlist(umap_plot$p_min_second)
umap_plot %>% head()
```

```R
umap_plot %>% arrange(p_min, p_min_second) %>% head(20)
```

There are still equalities at second pval... Could we just sum the -log10(pval)?


```R
umap_plot = umap_plot %>% dplyr::mutate(signif_sum_p_smaller = map(p_smaller_list, ~sum(-log10(.x)))) %>%
    dplyr::mutate(signif_sum_p_greater = map(p_greater_list, ~sum(-log10(.x))))

umap_plot$signif_sum_p_smaller = unlist(umap_plot$signif_sum_p_smaller)
umap_plot$signif_sum_p_greater = unlist(umap_plot$signif_sum_p_greater)

umap_plot$signif_sum_all = umap_plot$signif_sum_p_smaller + umap_plot$signif_sum_p_greater
umap_plot %>% head()
```

```R
umap_plot %>% arrange(p_min, p_min_second) %>% head(20)
```

CASP8AP2: required for S phase transition: https://www.pnas.org/doi/full/10.1073/pnas.0604227103
CIAPIN1: ?
COX15: increase cell survival when expressed: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2164945/
MED12 and MED19: required to initiate gene transcription
NAA20: required for cell cycle progresion https://portlandpress.com/biochemj/article/415/2/325/44493/Identification-of-the-human-N-acetyltransferase
RARS??
TTI1: regulates DNA damage response
CHMP6: required to complete abcission: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4230781/


```R

```

```R
top_perturbing_genes = umap_plot %>% arrange(desc(signif_sum_all)) %>% head(15)
top_perturbing_genes
```

HSPA9: mortalin, controls cell proliferation
SUPT6H: transcription
EIF3H: translation
SMC1A: structural maintenance of chromosomes
MCM3: cell 

Take top 15, contains ATRIP


MED12 and MED19


This is great, I now have a means to summarize the effect on the perturb-seq using a single metric.

```R
umap_plot$signif_sum_all %>% quantile(., probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
```

```R
umap_plot %>% dplyr::filter(p_min < 0.05) %>% pull(signif_sum_all) %>% quantile(., probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
```

```R
umap_plot %>% dplyr::filter(is_kzfp) %>% dplyr::arrange(desc(signif_sum_all))
```

```R
umap_plot %>% dplyr::filter(symbol == "ZNF274")
```

```R
res_stats %>% head()
```

```R
res_stats$signif_sum_all = NULL
res_stats = res_stats %>% dplyr::left_join(., umap_plot %>% dplyr::select(symbol, signif_sum_all) %>% unique(), by = c("condition" = "symbol"))
res_stats %>% head()
```

```R
res_stats %>% dplyr::group_by(bins) %>% arrange(p_smaller, desc(signif_sum_all)) %>% slice_head(n=5)
```

```R
# finding top 3 genes in each bin, up
top_n_genes = 3
top_genes_df = rbind(res_stats %>% dplyr::group_by(bins) %>% arrange(p_smaller, desc(signif_sum_all)) %>% slice_head(n=top_n_genes),
                          res_stats %>% dplyr::group_by(bins) %>% arrange(p_greater, desc(signif_sum_all)) %>% slice_head(n=top_n_genes))
```

```R
top_genes_df %>% dplyr::pull(condition) %>% unique() %>% sort()
```

- ABCE1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4889273/
- ADAM10 si induces cell cycle arrest: https://cancerci.biomedcentral.com/articles/10.1186/s12935-020-01727-5
- AMBRA1: Substrate-recognition component of a DCX (DDB1-CUL4-X-box) E3 ubiquitin-protein ligase complex involved in cell cycle control and autophagy (PubMed:20921139, 23524951, 24587252, 33854232, 33854235, 33854239, 32333458). The DCX(AMBRA1) complex specifically mediates the polyubiquitination of target proteins such as BECN1, CCND1, CCND2, CCND3, ELOC and ULK1 (PubMed:23524951, 33854232, 33854235, 33854239). Acts as an upstream master regulator of the transition from G1 to S cell phase: AMBRA1 specifically recognizes and binds phosphorylated cyclin-D (CCND1, CCND2 and CCND3), leading to cyclin-D ubiquitination by the DCX(AMBRA1) complex and subsequent degradation (PubMed:33854232, 33854235, 33854239) (genecards)
- ATP5F1A: ATP synthase (respiratory chain)
- ATR: DNA repair
- ATRIP: DNA repair
- BGLAP: ?
- BLZF1: Golgi function, ?
- BNIP1: BCL2 interactor, apoptosis, ER golgi ?
- CASP8AP2: required for S phase transition: https://www.pnas.org/doi/full/10.1073/pnas.0604227103
- CCND3: cyclin D
- CCNF: Cyclin F
- CCNK: Cyclin K
- CHAMP6: required to complete abcission: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4230781/
- CIAPIN1: ?
- COX15: increase cell survival when expressed: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2164945/
- DAD1: defender against apoptosis, loss induces apoptosis
- DBR1: mRNA splicing
- DCTN4: part of the dynein/dynactin, required for proper mitosis: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2173046/
- DDX20: splicing
- DENR: translation initiation
- DHDDS: lipid metabolism
- DHX15: mRNA splicing
- DNAJA3: protein folding
- DNAJC19:
- DHODH: De novo pyrimidine synthesis -> growth
- DYNLL2: cytoplasmic microtubules
- ElFs: initiation of translation
- ENTPD6: metabolism, ?
- ERCC2: DNA repair
- ETV4: ETS TF, accelerates G1/S transition: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5805611/
- EXOSC2: exosome?
- EXOSC8:
- FARSA: tRNA processing, translation
- FARSB: 
- FIGLA: ?
- GAB2: adaptor protein, proliferation
- GFER: mitochondrial metabolism
- GFM1: mitochondrial translation
- GTF2H1: general TF, DNA repair
- GTF3C2: 
- GTF3C3: 
- GEMIN5: Part of The SMN complex which catalyzes the assembly of small nuclear ribonucleoproteins (snRNPs), the building blocks of the spliceosome, and thereby plays an important role in the splicing of cellular pre-mRNAs (PubMed:16857593, 18984161, 20513430, 33963192) (genecards)
- GINS4: DNA replication initiation
- HARS: tRNA
- HBS1L: mRNA surveillance and repair
- HSD17B12: lipid metabolism
- INTS2: integrator complex subunit, RNA processing
- KANSL2: NSL complex, cell division
- HSPA9: mortalin, controls cell proliferation
- ITPR2: calcium channel receptor, may induce P53 etc https://www.nature.com/articles/s41467-021-20993-z
- KLHL17?
- LIN28A: let7 miRNA regulator
- MAD2L1: mitotic spindle assembly
- MCM2: DNA replication initiation
- MCM3: DNA replication initiation
- MED12: required to initiate gene transcription, MEDIATOR complex linked to cell cycling: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4386457
- MED19: required to initiate gene transcription
- MRPL19: mitochondrial ribosomal protein -> probably metabolism
- MTBP: MDM2 binding protein, regulates the cell cycle
- MYBBP1A: MYB related TF
- NAA20: required for cell cycle progresion https://portlandpress.com/biochemj/article/415/2/325/44493/Identification-of-the-human-N-acetyltransferase
- NAE1: required for S/M checkpoint
- NCBP: RNA processing
- NELFE: RNA regulation
- NUBP1: centrosome regulation: https://pubmed.ncbi.nlm.nih.gov/23028652/
- NUDT15: DNA damange repair
- NUP133: probably essential (nucleoporin)
- OXA1L: mitochondria
- PARPBP: DNA repair
- PGAM1: glycolysis
- PHB, PHB2: mitochondrial respiration
- PMPCA: mitochondria function
- PSMD13: proteasome -> cell cycle
- RAD21: DNA repair
- QARS: Glutaminyl-tRNA synthetase -> translation
- RARS: ?
- RFC4: DNA replication
- RINT1: DNA repair, cell cycle progression, telomere maintenance
- RPAP2: RNA pol II transcription
- RNF138: DNA damage response
- RPL24: ribosomal subunit -> translation, essential
- SACM1L: golgi and mitotic spindle
- SAMM50: mitochondria
- SIN3A: TF of cell cycle genes
- SLBP: regulator of histone protein expression
- SMC1A: sister chromatid cohesion, mitosis
- SMN2: SMN complex subunit -> splicing
- SPDL1: mitotic spindle formation
- SRP68: SRP = ER and Golgi
- SRRM1: mRNA processing
- STXBP5: vesicle trafficking
- SUPT5H: mRNA processing
- TAF2: GTF, RNA pol II
- TIMM23B: mitochondrial function
- TRNT1: tRNA translation
- VARS: tRNA valine
- VPS29: golgi trafficking
- WDHD1: DNA replication
- WDR82: SET1 subunit, chromatin maintenance cell cycle
- XPO1: nuclear export
- 
- TFDP1: E2F intecator
- TRRAP: TP53, E2F1 and E2F4 interactor for regulation
- TTI1: regulates DNA damage response
- YBX2: ? mRNA stability?
- YEATS4: involved in the cell cycle-regulating chromatin modifying complex NuA4  
- ZCRB1: splicing
- ZNHIT6: RNA processing

Summary by categories:
- 

```R
res_stats %>% dplyr::filter(as.integer(bins) == 1) %>% arrange(p_smaller) %>% head(5)
```

```R
library(ggrastr)
```

```R
#plotting the top 10 genes in term of significant imbalance on a UMAP that also has -log10(p) on it
set.seed(5)
library(viridis)
sc <- scale_colour_gradientn(colours = viridis::viridis(100, begin = 0.35), limits=range(-log10(0.05), -log10(0.001)), oob = scales::squish, na.value = "grey90")

p = umap_plot %>% dplyr::mutate(p_min = ifelse(umap_plot$p_min > 0.05, yes = NA, no = umap_plot$p_min)) %>% ggplot(aes(x = UMAP1, y = UMAP2, color = -log10(p_min))) + 
    ggrastr::geom_point_rast(size = 0.4) + 
    ggrepel::geom_text_repel(data = umap_plot %>% arrange(p_min) %>% dplyr::filter(symbol %in% top_genes_df$condition), 
                             aes(x = UMAP1, y = UMAP2, label = symbol), 
                             color = "black",
                            min.segment.length = unit(0, 'lines'),
                             segment.color = "grey60",
                             segment.linetype = 2,
                            force = 5,
                            force_pull = 1,
                            max.overlaps = Inf) + 
                            theme_classic() + theme(axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_blank(),
                           legend.position = "right") + sc
```

```R
p
```

```R
ggsave('../out/cell_cycle_figures/umap_cell_cycle_dev_signif_top_genes.svg', p, device = svg, width = 5, height = 4)
```

```R
# extreme imbalances: p < 1e-5 (limit of detection)
res_stats %>% dplyr::filter(p_min == 1e-5) %>% dplyr::pull(condition) %>% unique() %>% sort()
```

```R
res_stats %>% dplyr::filter(p_min == 1e-5) %>% dplyr::pull(condition) %>% unique() %>% sort() %>% length()
```

```R
library(Seurat)

umap_plot %>% dplyr::filter(symbol %in% cc_genes, p_min < 0.05) %>% pull(symbol) %>% unique() %>% length()
```

```R
cc_genes %in% umap_plot$symbol %>% sum()
```

```R
umap_plot %>% dplyr::filter(symbol %in% cc.genes$s.genes, p_min < 0.05) %>% pull(symbol) %>% unique() %>% length()
```

```R
cc.genes$s.genes %in% umap_plot$symbol %>% sum()
```

There seems to be an overrepresentation of seurat cell cycle marker genes in those causing imbalances, though not a strong one. More pronounced for S genes.

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = symbol %in% cc_genes)) + geom_point() + 
geom_text(data = umap_plot %>% dplyr::filter(symbol %in% cc_genes), aes(x = UMAP1, y = UMAP2, label = symbol))
ggsave('../out/cell_cycle_figures/UMAP_cell_positions_cc_genes.pdf', plot = p, device = pdf, width = 20, height = 20)
```

```R
# all gene names
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2)) + geom_point() + 
ggrepel::geom_text_repel(data = umap_plot, aes(x = UMAP1, y = UMAP2, label = symbol), max.overlaps = 100, max.iter = 1e6)
ggsave('../out/cell_cycle_figures/UMAP_cell_positions_all_genes.pdf', plot = p, device = pdf, width = 200, height = 200, limitsize = F)
```

```R
umap_plot %>% colnames()
```

### Plotting the position of KZFPs on the UMAP, comparison with age

```R
umap_plot$tf_category %>% table()
```

```R
umap_plot = umap_plot %>% dplyr::mutate(tf_category_signif = ifelse(umap_plot$p_min < 0.05, yes = umap_plot$tf_category, no = "n.s."))
umap_plot$tf_category_signif = factor(umap_plot$tf_category_signif, levels = c("KZFP", "ZFP", "other TF", "not a TF", "n.s."))
```

```R
tf_category_signif_colors = c(gene_category, "grey90", "grey60")
names(tf_category_signif_colors)[[4]] = "n.s."
names(tf_category_signif_colors)[[5]] = "other TF"
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (tf_category_signif))) + 
    ggrastr::geom_point_rast(size = 0.1) + 
    theme_classic() + 
    theme(axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            axis.text = element_blank(),
                            axis.title = element_blank(),
           legend.title=element_blank()) + 
    scale_color_discrete(type = tf_category_signif_colors) + 
    guides(colour = guide_legend(override.aes = list(size=1)))
```

```R
p
```

```R
ggsave('../out/cell_cycle_figures/UMAP_cell_positions_KZFPs.svg', plot = p, device = svg, width = 4, height = 3.5)
```

### Proportion of KZFPs vs. other TFs

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

Age of KZFPs with signif perturbation

```R
kzfps_jonas = read.csv('../data/kzfps_jonas.tsv', sep = '\t')
kzfps_jonas %>% head()
```

```R
kzfps_rhythmic_and_olga[!kzfps_rhythmic_and_olga %in% kzfps_jonas$assigned_gene]
```

```R
kzfps_jonas$age_MA %>% summary()
```

```R
cut(kzfps_jonas$age_MA, breaks = c(0, 43.20, 105, 352))
```

```R
tab <- data.frame(xtabs(~ cut + color, data = diamonds))
head(tab)
```

```R
to_plot = res_stats %>% dplyr::filter(is_kzfp) %>% 
    dplyr::select(condition, p_min) %>% 
    dplyr::group_by(condition) %>% 
    dplyr::mutate(p_min = min(p_min)) %>% unique() %>%
    dplyr::left_join(., kzfps_jonas %>% dplyr::select(assigned_gene, age_MA), by = c("condition" = "assigned_gene")) %>% 
    dplyr::mutate(age_binned = cut(age_MA, breaks = c(0, 21, 76, Inf)))
                                   #breaks = c(0, 25, 30, 60, 105, 352)))
```

```R
to_plot %>% ggplot(aes(x = age_binned, fill = p_min < 0.05)) + geom_bar(position = "fill")
```

TODO: discuss with Romain

```R
# age vs leads to imbalance, all KZFPs, via a permutation test for the difference between mean ages.

# too many ties for a wilcoxon

age_means= to_plot %>% 
    dplyr::filter(!is.na(age_MA)) %>%
    dplyr::mutate(causes_imbalance = p_min <0.05) %>%
    dplyr::group_by(causes_imbalance) %>% 
    dplyr::summarize(age_mean = mean(age_MA), n_kzfps = n())
age_means
```

```R
age_means_diff_observed = diff(age_means$age_mean)
age_means_diff_observed
```

Now, running permutations: we split the KZFPs in groups of the same size as the one observed and compute the difference in mean age. 

```R
to_permute = to_plot %>% 
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
p = to_plot %>% ggplot(aes(x = p_min < 0.05, y = age_MA)) + 
    geom_violin() + 
    geom_jitter_rast(width = 0.1, height = 0, col = "black", size = 0.1) + 
    theme_classic() + 
    ylab("Age (MY)") +
    theme(axis.line = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x = element_blank()) + 
    scale_x_discrete(labels = c("no change", "imbalance"))
p
```

```R
ggsave('../out/cell_cycle_figures/imbalance_vs_age.svg', plot = p, device = svg, width = 2.5, height = 1.5)
```

```R
top_plot %>% dplyr::filter(age_MA < 30, p_min < 0.05)
```

```R
top_plot$age_MA %>% table()
```

```R
umap_plot %>% dplyr::filter(is_kzfp, p_min < 0.05) %>% arrange(desc(signif_sum_all)) %>% head(20)
```

```R
umap_plot %>% dplyr::filter(symbol %in% c("ZNF92", "ZNF850", "ZNF519", "ZNF93", "ZNF695", "ZNF184", "ZNF473", "ZFP69B", "ZNF530", "ZNF714"))
```

KZFPs that are associated to proliferation are overrepresented among those leading to imbalances. (6/10)

```R
umap_plot %>% dplyr::filter(symbol %in% c("ZNF540", "ZFP2", "ZNF688", "ZNF846", "ZNF671", "ZNF25", "ZNF582", "ZNF208", "ZNF415", "ZNF483"))
```

```R
umap_plot %>% dplyr::filter(symbol %in% c("ZNF540", "ZFP2", "ZNF688", "ZNF846", "ZNF671", "ZNF25", "ZNF582", "ZNF208", "ZNF415", "ZNF483"))
```

5/10 KZFPs anticorrelated with proliferation lead to imbalances.

```R
11/20
```

```R
126/348
```

| | Top KZFPs associated to prolif | KZFPs not among top associated to prolif | total |
|---|:---:|:---:|:---:|
| KZFPs leading to imbalances | x | m-x | m |
| KZFPs not leading to imbalances | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
m = umap_plot %>% dplyr::filter(is_kzfp, p_min < 0.05) %>% pull(symbol) %>% unique() %>% length()
m 
```

```R
n = umap_plot %>% dplyr::filter(is_kzfp, p_min >= 0.05) %>% pull(symbol) %>% unique() %>% length()
n
```

```R
k = 20
x_observed = 9
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

It's barely not significant, unlucky


KZFPs from Filipe's figure 1:
| | KZFPs associated with survival | KZFPs not associated with survival | total |
|---|:---:|:---:|:---:|
| KZFPs leading to imbalances | x | m-x | m |
| KZFPs not leading to imbalances | k-x | n-(k-x) | n |
| total | k | (m+n-k) | m+n |

```R
x_observed = 17
k = 41-4
x = 0:m # is the variable tested
probs <- dhyper(x, m, n, k, log = FALSE)


# we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
pval_one_sided = sum(probs[x>=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])
pval_two_sided
```

Not significant either.


What about the telomere / centromere binders found by Evarist by aligning to T2T?

```R
umap_plot %>% dplyr::filter(symbol %in% c("ZNF671", "ZNF189", "ZNF429", "ZNF397", "ZNF449", "ZNF555", "ZNF792", "ZNF582", "ZNF823", "ZNF485", "ZNF805"))
```

```R
umap_plot %>% dplyr::filter(symbol %in% c("ZNF223", "ZNF792", "ZNF823", "ZNF224", "ZNF671", "ZNF189", "ZNF429",  "ZNF530", "ZNF19", "ZNF354C", "ZNF727", "ZNF445", "ZNF749", "ZNF223"))
```

```R
umap_plot %>% dplyr::filter(symbol %in% c("ZNF671", "ZNF429",  "ZNF530", "ZNF19", "ZNF354C", "ZNF727", "ZNF445", "ZNF749")) %>% mutate(signif = p_min < 0.05) %>% pull(signif) %>% table()
```

```R
x_observed = 5
k = 8

x = 0:m # is the variable tested
probs <- dhyper(x, m, n, k, log = FALSE)


# we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
pval_one_sided = sum(probs[x>=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])
pval_two_sided
```

That seems much more promising. If we take outliers, they are enriched. What if we take all the ones evarist highlighted in red?

```R
telo_centro = read.csv('../data/all_telomeric_centromeric_counts.csv') %>% dplyr::rename(sample_name = X)
telo_centro %>% head()
```

```R
telo_centro$symbol = telo_centro$sample_name %>% gsub("_pubM|_lowconf|_n_noid|_rep|rep|_n|oid|_1|_2|ocorr", "",  .)
```

```R
telo_centro$symbol %>% unique() %>% print()
```

```R
telo_centro_kzfps = telo_centro %>% dplyr::filter(log10(telo.frac.u) > 1.85 | log10(centro.frac.u) >4.85, symbol %in% kzfps_rhythmic_and_olga) %>% pull(symbol) %>% unique()
```

```R
telo_centro_kzfps
```

```R
umap_plot %>% dplyr::filter(symbol %in% telo_centro_kzfps) %>% mutate(signif = p_min < 0.05) %>% pull(signif) %>% table()
```

```R
x_observed = 13
k = 26

x = 0:m # is the variable tested
probs <- dhyper(x, m, n, k, log = FALSE)


# we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
pval_one_sided = sum(probs[x>=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])
pval_two_sided
```

No enrichment. What if we take the top n KZFPs that are either centro or telo?

```R
top_k = 5
top_telo_centro = c(telo_centro %>% dplyr::filter(symbol %in% kzfps_rhythmic_and_olga) %>% arrange(desc(log10(telo.frac.u))) %>% head(top_k) %>% pull(symbol),
                    telo_centro %>% dplyr::filter(symbol %in% kzfps_rhythmic_and_olga) %>% arrange(desc(log10(centro.frac.u))) %>% head(top_k) %>% pull(symbol)) %>% unique()
```

```R
top_telo_centro %>% unique()
```

```R
umap_plot %>% dplyr::filter(symbol %in% top_telo_centro) %>% mutate(signif = p_min < 0.05) %>% pull(signif) %>% table()
```

```R
x_observed = 7
k = 10

x = 0:m # is the variable tested
probs <- dhyper(x, m, n, k, log = FALSE)


# we compute the probability of observing a more extreme enrichment, therefore using a one sided test. 
pval_one_sided = sum(probs[x>=x_observed])

# we make the test two-sided, by summing the probabilities that are smaller or equal to our pval
pval_two_sided = sum(probs[probs <= pval_one_sided])
pval_two_sided
```

## Gene set enrichment for genes leading to imbalance in cell population

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
library(GOSemSim)
```

```R
p = dotplot(GOcluster_all_imbalances) + 
    scale_colour_gradientn(colours = rev(viridis(100)), limits=c(0, 0.1)) + xlab(NULL)
p
```

```R
ggsave('../out/cell_cycle_figures/go_perturbseq_all_genes_leading_to_imbalance.svg', p, svg, width = 5, height = 3)
```

```R
# exporting data for Romain
GOcluster_all_imbalances %>% as.data.frame() %>% write.table('../out/tables/go_perturbseq_all_genes_leading_to_imbalance.tsv', col.names = T, row.names = F, sep = '\t')
```

```R
t
```

```R
GOcluster_per_phase_both_directions <- clusterProfiler::compareCluster(condition~phase, 
                            data = res_stats %>% dplyr::filter(p_min < 0.05) %>% dplyr::select(condition, phase) %>% dplyr::mutate(across(condition, as.character)), 
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
clusterProfiler::dotplot(GOcluster_per_phase_both_directions)
```

```R
GOsimpl = clusterProfiler::simplify(GOcluster_per_phase_both_directions, 0.1)
clusterProfiler::dotplot(GOsimpl)
```

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
    scale_colour_gradientn(colours = rev(viridis(100)), limits=c(0, 0.1)) + xlab(NULL)
p
```

```R
ggsave('../out/cell_cycle_figures/go_perturbseq_genes_leading_to_phase_depletion.svg', p, svg, width = 6, height = 8)

# ggsave('../out/cell_cycle_figures/go_perturbseq_genes_leading_to_phase_depletion.emf', p, device = devEMF::emf, width = 6, height = 8)
```

```R
GOsimpl = clusterProfiler::simplify(GOcluster_per_phase_depletion, 0.8)
clusterProfiler::dotplot(GOsimpl)
```

A depletion could mean either: acceleration to the next phase OR death during that phase. 

Depletion in phase i, accumulation in phase i+1: acceleration
depletion in phase i, accumulation in phase i-1: death. Let's check that hypothesis


Genes leading to accumulations in phases upon knockdown

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
p = clusterProfiler::dotplot(GOcluster_per_phase_enrichment) + scale_colour_gradientn(colours = rev(viridis(100)), limits=c(0, 0.1)) + xlab(NULL)
p
```

```R
ggsave('../out/cell_cycle_figures/go_perturbseq_genes_leading_to_phase_accumulation.svg', p, svg, width = 6, height = 5)
```

```R
GOsimpl = clusterProfiler::simplify(GOcluster_per_phase_depletion, 0.1)
clusterProfiler::dotplot(GOsimpl)
```

Conclusion: we don't simplify, it adds spurious terms that have nothing to do with what we want.


## Vicinities:

Are there any position similarities between KZFPs and genes leading to accumulation in particular phases upon KO?

Are there position similarities between transcripts that are rhythmic and genes leading to accumulation during particular phases of the cell cycle?

- ZNF695: DNA damage / replication stress genes (ATR/ATRIP), general TFs and RNA transcription (TAFs)
- ZKSCAN3: proteasome and golgi trafficking
- 

```R
fit_all_genes = read.table('../out/tables/all_genes_cycling_info.tsv', sep = '\t', header = T)
```

```R
fit_all_genes %>% head()
```

```R
umap_plot = umap_plot %>% dplyr::left_join(., fit_all_genes %>% dplyr::filter(padj < 0.05) %>% dplyr::select(symbol, phase_assigned, padj))
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (phase_assigned == "eG1" & padj < 0.05))) + geom_point()
p
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (phase_assigned == "eG1" & padj < 0.01))) + geom_point()
p
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (phase_assigned == "lG1" & padj < 0.01))) + geom_point()
p
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (phase_assigned == "G1/S" & padj < 0.001))) + geom_point()
p
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (phase_assigned == "S1" & padj < 0.001))) + geom_point()
p
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (phase_assigned == "S2" & padj < 0.001))) + geom_point()
p
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (phase_assigned == "G2" & padj < 0.001))) + geom_point()
p
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (phase_assigned == "G2/M" & padj < 0.001))) + geom_point()
p
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = (padj < 0.001))) + geom_point()
p
```

The genes peaking in a certain phase do not seem to accumulate in a specific place of the UMAP...

```R
# is there an effect of the number of cells in each phase?
umap_plot = umap_plot %>% dplyr::left_join(., y = cells_per_cond, by = c("symbol" = "condition"))
```

```R
umap_plot %>% head()
```

```R
library(viridis)
sc <- scale_colour_gradientn(colours = viridis(100), limits=c(50, 500))
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = ncells)) + ggrastr::rasterize(geom_point(size = 0.5), dpi = 300) + sc + theme_classic() + theme(axis.line=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks=element_blank(),
                                      axis.title.x=element_blank(),
                                      axis.title.y=element_blank()
                                 )
ggsave('../out/cell_cycle_figures/umap_cell_cycle_ncells.svg', plot = p, device = svg, width = 4, height = 4)
```

UMAP borders always have low cell numbers, and the center generally has more cells. Again, I bet the center conditions are less damaging to the cell cycle than the ones on the edge.


UMAP plots and deviations from the expected number of cells:

```R
umap_plot %>% head()
```

```R
res_stats %>% dplyr::filter(condition== "NFYB")
```

```R
umap_plot %>% head()
```

```R
res_stats %>% 
    dplyr::select(bins, condition, p_smaller, p_greater) %>% 
    dplyr::arrange(condition, bins) %>% head(10)
```

```R
umap_plot %>% head()
```

```R
umap_plot %>% dplyr::select(p_smaller_list) %>% .[[1]]
```

```R
umap_plot %>% dim()
```

```R
umap_plot %>% apply(., 1, function(x) x["p_smaller_list"][[1]][[1]], simplify = T) %>% unlist() %>% 
```

```R
which((umap_plot %>% apply(., 1, function(x) length(x["p_smaller_list"][[1]]), simplify = T) %>% unlist())<5)
```

```R
umap_plot[4945, ]
```

This one is preventing the plot from being made by failing the mutate. We exclude it.

```R
dir.create('../out/cell_cycle_figures/umaps_cell_cycle_imbalances/')
```

```R
breaks
```

```R
library(viridis)
sc <- scale_colour_gradientn(colours = viridis::viridis(100), limits=range(0, 2), oob = scales::squish)

# list of plots
plot_list = list(p_smaller_list = list(), p_greater_list = list())
```

```R
# we're arranging the plots in two rows and four columns. Rows are greater or smaller. Columns are phases, from M/eG1 to G2/M.
# We only add x label to the first row of plots with the phases, and y labels to the two first plots of each column.

phase_assignment

direction_labels = list("p_greater_list"="accumulation", "p_smaller_list"="attrition")
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
                            geom_point(size = 0.2) + 
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
prow = plot_grid(
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

ggsave('../out/cell_cycle_figures/umap_cell_cycle_dev_signif.pdf', plot_grid(prow, legend, rel_widths = c(3, .4)), device = pdf, width = 12, height = 5)
```

## Does depletion in a phase mean a quicker exit from that phase, or death during that phase?

To answer that question, we test whether the phase i+1 is accumulated (indicating a quicker exit) or whether phase i-1 is depleted (indicating death).

```R
# testing accumulation
phase_assignment
```

```R
accelerated_in_depleted_phase = list()
for (i in 1:length(phase_assignment)) {
    phase_depleted = phase_assignment[i]
    phase_tested_accumulation = phase_assignment[(i+1)]
    if(i == length(phase_assignment)) {
        phase_tested_accumulation = phase_assignment[1]
        }    
    genes_leading_to_depl_i = res_stats %>% 
        dplyr::filter(p_smaller < 0.1, phase == phase_depleted) %>% 
        dplyr::select(condition) %>% pull() %>% as.character()
    
    genes_leading_to_acc_i_plus_1 = res_stats %>% dplyr::filter(p_greater < 0.1, phase == phase_tested_accumulation) %>% dplyr::select(condition) %>% pull() %>% as.character()
    
    accelerated_in_depleted_phase[[phase_assignment[i]]] = genes_leading_to_depl_i[which(genes_leading_to_depl_i %in% genes_leading_to_acc_i_plus_1)]
    }
```

```R
GOcluster_per_phase_depletion_is_acceleration <- clusterProfiler::compareCluster(accelerated_in_depleted_phase, 
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
p = clusterProfiler::dotplot(GOcluster_per_phase_depletion_is_acceleration) + 
    scale_colour_gradientn(colours = rev(viridis(100)), limits=c(0, 0.1)) + xlab(NULL)
ggsave('../out/cell_cycle_figures/go_perturbseq_genes_leading_to_phase_depletion_acceleration.svg', p, svg, width = 5, height = 4)
```

Does this mean that removing transcriptional regulators accelerates going from M/eG1 to lG1? TODO: check the genes by hand.

```R
died_in_depleted_phase = list()
for (i in 1:length(phase_assignment)) {
    phase_depleted = phase_assignment[i]
    phase_tested_accumulation = phase_assignment[(i-1)]
    if(i == 1) {
        phase_tested_accumulation = phase_assignment[length(phase_assignment)]
        }    
    genes_leading_to_depl_i = res_stats %>% 
        dplyr::filter(p_smaller < 0.1, phase == phase_depleted) %>% 
        dplyr::select(condition) %>% pull() %>% as.character()
    
    genes_leading_to_acc_i_minus_1 = res_stats %>% dplyr::filter(p_greater < 0.1, phase == phase_tested_accumulation) %>% dplyr::select(condition) %>% pull() %>% as.character()
    
    died_in_depleted_phase[[phase_assignment[i]]] = genes_leading_to_depl_i[which(genes_leading_to_depl_i %in% genes_leading_to_acc_i_minus_1)]
    }
```

```R
GOcluster_per_phase_depletion_is_death <- clusterProfiler::compareCluster(died_in_depleted_phase, 
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
p = clusterProfiler::dotplot(GOcluster_per_phase_depletion_is_death) + 
scale_colour_gradientn(colours = rev(viridis(100)), limits=c(0, 0.1)) + xlab(NULL)
ggsave('../out/cell_cycle_figures/go_perturbseq_genes_leading_to_phase_depletion_death.svg', p, svg, width = 5, height = 4)
```

Does this mean that removing transcriptional regulators kills cells in G1/S? Also, the terms are very similar to the previous enrichment test. That's fishy.


This will be nicer with the velocycle velocity model, where speed will be compared.

```R
phase_assignment[0]
```

```R
p = umap_plot %>% ggplot(aes(x = UMAP1, y = UMAP2, color = symbol %in% accelerated_in_depleted_phase[[4]])) + geom_point()
p
```

```R
res_stats %>% dplyr::filter(p_smaller < 0.05, phase == "G1/S") %>% dplyr::inner_join(., res_stats %>% dplyr::filter(p_smaller < 0.05, phase == "lG1"))
```

```R

```

Why don't we see anything for genes accumulating in S? We saw a clear enrichment for DNA replication on GOBP via ENRICHR

```R
GOcluster_enrichment_S <- clusterProfiler::compareCluster(geneClusters = list(enriched_in_S = res_stats %>% dplyr::filter(p_greater < 0.05, phase == "S") %>% dplyr::select(condition) %>% dplyr::mutate(across(condition, as.character)) %>% pull()),
                            ont = "BP", 
                            keyType = 'SYMBOL',
                            universe = for_go %>% dplyr::select(condition) %>% unique() %>% pull() %>% as.character(), 
                            OrgDb = org.Hs.eg.db, 
                            fun="enrichGO",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 1,
                            qvalueCutoff  = 1,
                            minGSSize = 10,
                            maxGSSize = 500
)
```

```R
GOcluster_enrichment_S %>% dplyr::arrange(p.adjust) %>% head(5)
```

```R
GOcluster_depleted_S <- clusterProfiler::compareCluster(geneClusters = list(enriched_in_S = res_stats %>% dplyr::filter(p_smaller < 0.05, phase == "S") %>% dplyr::select(condition) %>% dplyr::mutate(across(condition, as.character)) %>% pull()),
                            ont = "BP", 
                            keyType = 'SYMBOL',
                            universe = for_go %>% dplyr::select(condition) %>% unique() %>% pull() %>% as.character(), 
                            OrgDb = org.Hs.eg.db, 
                            fun="enrichGO",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 1,
                            qvalueCutoff  = 1,
                            minGSSize = 10,
                            maxGSSize = 500
)
```

```R
GOcluster_depleted_S %>% dplyr::arrange(p.adjust) %>% head()
```

```R
GOcluster_enrichment_S
```

Indeed the symbols lost by compareCluster, reinput into enrichr, fail to recover DNA replication. That is annoying...

```R
GOcluster <- clusterProfiler::compareCluster(geneClusters = list(`all imbalances` = res_stats %>% dplyr::filter(p_min < 0.05) %>% dplyr::select(condition) %>% unique() %>% pull()), 
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
# converting to entrezid
gene_metadata = read.csv('../data/cell_cycle/2110_hg19_annotation.csv', sep= ',') %>% 
    dplyr::filter(!(is.na(entrezgene_id) | hgnc_symbol == "" | is.na(hgnc_symbol)))
gene_metadata = gene_metadata[!duplicated(gene_metadata$hgnc_symbol), ]
rownames(gene_metadata) = gene_metadata$hgnc_symbol
gene_metadata %>% head()   
gene_metadata %>% dim()
```

```R
GOcluster <- clusterProfiler::compareCluster(geneClusters = list(`all imbalances` = gene_metadata[res_stats %>% dplyr::filter(p_min < 0.05) %>% dplyr::select(condition) %>% unique() %>% pull(), "entrezgene_id"]), 
                            universe = gene_metadata[res_stats %>% dplyr::select(condition) %>% unique() %>% pull(), "entrezgene_id"], 
                            fun="enrichKEGG",
                            keyType = "ncbi-geneid",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.1,
                            minGSSize = 10,
                            maxGSSize = 500
)
```

```R
dotplot(GOcluster)
```

```R
# doesnt yield anything sensible
# enrichR package does not offer the option to define a background set of genes, which sucks.
```

## KZFPs that are both rhythmic and leading to imbalances in phases

```R
fit_kzfps %>% head()
```

```R
kzfps_imbalances_vs_rhythm = fit_kzfps %>% dplyr::left_join(., res_stats, by = c("symbol" = "condition"))
```

```R
kzfps_imbalances_vs_rhythm %>% dim()
```

```R
kzfps_imbalances_vs_rhythm %>% dplyr::filter(p_min < 0.05) %>% dplyr::select(symbol) %>% unique() %>% pull() %>% length()
```

```R
kzfps_imbalances_vs_rhythm %>% dplyr::filter(p_min < 0.05, padj < 0.05) %>% dplyr::select(symbol) %>% unique() %>% pull() %>% length()
```

```R
kzfps_imbalances_vs_rhythm %>% dplyr::filter(p_min < 0.05, padj < 0.05) %>% arrange(phase_assigned) %>% dplyr::select(phase_assigned, symbol) %>% unique()
```

## KZFPs w.r.t. location in the heatmap, and the genes surrounding them. 

```R
umap_plot %>% head()
```

```R
res_stats %>% head()
```

### KZFPs accumulating in M/eG1, located in the striking clusters

```R
res_stats %>% dplyr::filter(phase == "M/eG1", p_greater <=0.05, condition %in% kzfps_rhythmic_and_olga) %>% arrange(p_greater)
```

ZNF671 clusters with ATR and ATRIP and is the top centromere binder KZFP. ![image.png](attachment:0141fe38-ba81-4ea0-aa78-ea1cf8db20cf.png)


ZNF418 clusters with ZNF528, ZNF347, NME6 (inhibitor of P53 apoptosis), TUBGCP6 (centrosome), TAF3 (RNA pol II subunit), YY1AP1 (DNA repair), POLR1D (RNA pol I and III subunit) ![image.png](attachment:9be56aa4-5e2c-49c8-8c01-000a15ec1058.png)


## KZFPs accumulating in lG1, located at the top of the UMAP

```R
res_stats %>% dplyr::filter(phase == "lG1", p_greater <=0.05, condition %in% kzfps_rhythmic_and_olga) %>% arrange(p_greater)
```

ZNF669 colocates with proteins involved in membrane trafficking, golgi, lipid metabolism, protein translation at the membrane. Also NEK2, involved in mitotic regulation


![image.png](attachment:6476b52f-0482-4485-99b3-11a7b1d66996.png)


ZNF585B clusters with POC1B (centriole, mitoticl spindle), CDK2AP2 (NuRD complex, regulates G1/S), PSMD6 (proteasome), CCT3 (telomere maintenance), DDX39B (splicing), MARS (tRNA-peptide loading) ![image.png](attachment:8b373cd4-22ed-4f1c-9886-7613722645a1.png)


ZNF600 clusters with CDC73, POLR2F (transcription, cell cycle), BRD7 (P53-induced senecesence, DNA damage), NCAPH, CDCA5 (chromatin organization for mitosis), ATG5 (autophagy)![image.png](attachment:aaf6b423-6558-4ed2-95ad-3006f8b64562.png)


ZNF695: ACVR2B (signal transduction for survival), DIAPH1 (locates APC and CLASP to the cell membrane), DYNC1H1 (mitrotubule assembly), as well as ZSCAN29![image.png](attachment:f4333d7b-af6d-415c-aebb-0aa70d44369c.png)


ZKSCAN8 clusters with ODF2 (centrosome matrix assembly), EIF4B (translation), ENO1 (cmyc repressor), DDX6 putative overall KZFP regulator via translation regulation ![image.png](attachment:460d1c1e-fef7-445f-9e30-084d10638934.png)


### KZFPs leading to accumulation in G1/S

```R
res_stats %>% dplyr::filter(phase == "G1/S", p_greater <=0.05, condition %in% kzfps_rhythmic_and_olga) %>% arrange(p_greater)
```

ZNF577 clusters with TFDP1 (E2F binding partner, required for G1/S transition), E2F3, COQ10/DUFS1 (oxidative phosphorylation), ![image.png](attachment:8852ec60-dcee-4b14-b642-54626f5559f3.png)


ZNF430: ANO6, APOM, PLSCR3, ERG28 (to check) (membrane lipid metabolism), PPM1G (cell stress response, cell cycle progression), GSTZ1 (glutathione s transferase, limits cell stress), RPRD1B (RNA pol II holoenzyme, reguletes CCND1 and the cyclins) ![image.png](attachment:86789726-e63f-48b4-8452-58daf366065f.png)


ZNF850: ANAPC7 (APC subunit),  KMT2E (transcriptional regulation of the cell cycle), RBBP7 (binds RB), RHOH (RAS homolog, inhibits proteins of the DNA damage response), CENPH (mitotic spindle) ![image.png](attachment:91b9e0a9-270f-42bf-9c71-970b5abe1b41.png)


ZNF90: clusters with ZNF611, but not sure I understand what these genes are about ![image.png](attachment:d7f29875-5047-4833-9782-e4dc37bad98b.png)


ZNF732: close to ZNF471, associates with proteasome activators and subunits FAM192A, PSEM1, POLD4 (DNA replication)![image.png](attachment:4bad2b5a-4728-4f8c-a237-3629443160a1.png)


ZNF789 associates with DMC1 (DNA repair), POGZ (transposase-containing Zinc finger, involved in mitosis), HPF1 (DNA repair) ![image.png](attachment:bc59f156-36e9-4bd7-b521-d1d428dac337.png)


ZNF221 clusters with ZNF827, RECQL5 (helicase, Genome stability), CDK19, BCL6B, CENPM,  ![image.png](attachment:e6abe0aa-d9a0-469f-a5cb-3049441b1713.png)


### KZFPs leading to an accumulation in S (bottom right of UMAP)

```R
res_stats %>% dplyr::filter(phase == "S", p_greater <=0.05, condition %in% kzfps_rhythmic_and_olga) %>% arrange(p_greater)
```

ZNF84: clusters with lysine methyltransferase involved in active chromatin deposition (KMT2B), WDR90/B9D1 (cilia/centriole component), TUBB3 (major tubulin subunit) ![image.png](attachment:e36ee0f6-87fa-4506-a035-aaaaa8a2764d.png)


ZNF726, ZNF304, ZNF7 cluster with heterochromatin enzymes SEDTB1, KDM1B, chromatin modifier ENY2, TRAF3IP1 (negative regulator of microtubule stability), STUB1 (E3 ubibquitin ligase involved in DNA repair), KIF21A (microtubule motor) ![image.png](attachment:05886776-806d-49e2-930f-f99fb1550c3a.png)


ZKSCAN2: TOB1 (antiproliferative protein) ![image.png](attachment:8e8249f3-2a98-4358-a209-701893e60205.png)


ZNF69: DDB1 (UV reponsive DNA repair), TRIM56 (cGAS STING, dsDNA sensing, activates IFNB1), IL10RB (antiviral immunity, JAK/STAT). ![image.png](attachment:43764d67-a4a6-4612-bac0-2174d64e0f0c.png)


ZNF599, ZNF621 are close to RBL1 (RB ligand), ![image.png](attachment:475d1240-6b16-47e3-978c-91ff0e9a2302.png)


ZNF528: quite far on the center right of the UMAP. 


### KZFPs leading to an accumulation in G2/M (left of the UMAP)

```R
res_stats %>% dplyr::filter(phase == "G2/M", p_greater <=0.05, condition %in% kzfps_rhythmic_and_olga) %>% arrange(p_greater)
```

ZNF211: ?


ZNF655: clusters with VKR2 (inhibits BANF1), API1: apoptosis inhibitor ![image.png](attachment:8f153f78-733f-4815-ae5d-65fde7a805a6.png)


ZNF222, ZNF557, ZNF782 cluster with ORC1 (DNA replication origin), CCNE2 (cyclin E2), UBE3A (degrader of cell cycle players), SUGT1 (kinetochore assembly, required for G1/S and G2/M), PRIM1 (DNA primase subunit)![image.png](attachment:0b791f63-b1e4-4874-922e-3053ec3c98ba.png)


ZNF354C: TOPBP1/CDC45 (DNA replication, DNA damage response), CENPT (centromeres), PFKFB2/MRPL47 (metabolism), RPAP2/CWC15, CLP1 (RNA pol II holoenzyme, spliceosome), pretty strong telomere binder (Evaris)![image.png](attachment:dde95563-b5bd-4e2e-95b4-c25863cb00a3.png)


ZNF189: APEX2 (DNA repair), DNASE2, MTA2 (NuRD deacetylase subunit), CCT8 (chaperonin), is a high centromere binder![image.png](attachment:5a8e008e-3a53-426b-b20c-f08d2f124117.png)


ZNF415: ANKLE1 (DNA damage response), UBAP2L (ubiquitin response after DNA damage), DCUN1D2 (cullin),  ![image.png](attachment:dc85bc76-adf8-4b7d-b46d-07f03f4a731f.png)


ZNF257: DDX41 and DDX42 (RNA helicases, RNA splicing, TP53 resistance to apoptosis), NRDE2 (RNA splicing, prevents mRNA degradation in the nucleus), SF3B4 (splicin),  POLR2B (RNA polymerase subunit), POLDIP2 (DNA polymerase subunit, bypasses DNA lesions), TICRR (TOPBP1 interactor, DNA replication, S/M and G2/M checkpoints), CDC7 (cell cycle progression, S phase), PDCD2L (delays cell cycle progression),  ![image.png](attachment:ef6b81e8-aca4-4850-834e-5fca903acb62.png)


ZNF224: ZNF479, TP53BP1 (DNA repair), FAIM (Fas apoptosis) ![image.png](attachment:1942132a-7f33-43e5-9e03-bb1122d31753.png)


ZNF684: PCNA (DNA replication), NPM1 (centrosome, DNA replication, P53), DELE1 (integrated stress response), several proteins in RNA metabolism ![image.png](attachment:7dfba6fd-508e-48e2-861d-90e9caa1e28f.png)


## KZFPs leading to depletions in specific cell phases


### KZFPs leading to depletion in M/eG1

```R
res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga, phase == "M/eG1", p_smaller <= 0.05) %>% arrange(p_smaller)
```

ZNF655: already checked: accumulation in G2/M


ZNF331: BORA (AURKA activator), TAF1C (RNA POL I), HECTD3 (cell cycle progression), MIEN1 (regulation of apoptosis) ![image.png](attachment:7443d816-fa73-4f36-ac55-501c437f1535.png)


ZNF792: CCNDBP1 (Cyclin D1 binder), CENPO (centromere) ![image.png](attachment:8b8fc499-9047-4a23-98ad-cc42e42a19b2.png)


## KZFPs leading to depletions in lG1 (bottom left)

```R
res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga, phase == "lG1", p_smaller <= 0.05) %>% arrange(p_smaller)
```

ZNF684 leads to an accumulation in S, already characterized

```R
res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga, phase == "S", p_smaller <= 0.05) %>% arrange(p_smaller)
```

```R
res_stats %>% dplyr::filter(condition %in% kzfps_rhythmic_and_olga, phase == "G2/M", p_smaller <= 0.05) %>% arrange(p_smaller)
```

ZNF33A: SERTAD1 (E2F1 coregulator)


ZNF565: LIN9 (tumor suppressor, G1/S transition, inhibits DNA synthesis), RIF1, ![image.png](attachment:91933a00-5839-4c45-ab95-78d671dc949c.png)


ZNF528: custers with ZNF347, ZNF418. Of note: 418 and 528 are on the same KZFP cluster (19.9

```R

```
