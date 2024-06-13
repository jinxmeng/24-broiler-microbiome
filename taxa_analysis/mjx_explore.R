#### info ####
# Jinxin Meng, 2023-10-09, 2023-11-11
pacman::p_load(tibble, dplyr, tidyr, purrr, ggplot2, vegan, ggpmisc, ggpubr)
setwd("/share/data1/mjx/proj/01.broiler_metagenome_20230301/analysis/taxa")
source("/share/data1/mjx/Rproj/script/taxa_compos.R")
source("/share/data1/mjx/Rproj/script/profile_process.R")
source("/share/data1/mjx/Rproj/script/diversity.R")
source("/share/data1/mjx/Rproj/script/diff_test.R")
source("/share/data1/mjx/Rproj/script/utilities.R")

#### phylum composition ####
library(ComplexHeatmap)
group <- read.delim("../data/sample_group")
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% 
  data.frame() 
taxonomy <- read.delim("../data/species_taxonomy_rename.tsv")

group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

dat <- taxa_trans(otu, taxonomy, group, to = "feature", out_all = T) %>% 
  apply(., 2, \(x) ifelse(x > 0, 1, 0)) %>% 
  data.frame() %>% 
  rownames_to_column(var = "feature") %>% 
  merge(., taxonomy, by = "feature")
write.table(dat, "tmp_presence_heatmap.txt", sep = "\t", quote = F, row.names = F)

# phylum
plotdat <- select(dat, feature, all_of(group_order)) %>% 
  column_to_rownames(var = "feature")
row_split <- factor(dat$phylum)
pdf(file = "tmp_presence_heatmap.pdf", height = 6, width = 4)
pheatmap(plotdat, 
         cluster_rows = T, cluster_cols = T, treeheight_row = 0, treeheight_col = 0,
         cellwidth = 10, cellheight = 1, show_rownames = F, fontsize = 8,
         color = circlize::colorRamp2(c(1, 0), c("#d53e4f", "white")),
         border_gp = gpar(col = "black"),  
         row_split = row_split, row_gap = unit(0,'mm'), row_title_rot = 0, 
         legend = F)
dev.off()

#### taxa freq ####
dat <- taxa_trans(otu, taxonomy, group, to = "feature", out_all = T) %>% 
  apply(., 2, \(x) ifelse(x > 0, 1, 0)) %>% 
  data.frame() %>% 
  rownames_to_column(var = "feature") %>% 
  merge(., select(taxonomy, feature, family), by = "feature") %>% 
  mutate( d1 = ifelse(d1 == 1, family, NA),
          d3 = ifelse(d3 == 1, family, NA),
          d5 = ifelse(d5 == 1, family, NA),
          d7 = ifelse(d7 == 1, family, NA),
          d14 = ifelse(d14 == 1, family, NA),
          d21 = ifelse(d21 == 1, family, NA),
          d28 = ifelse(d28 == 1, family, NA),
          d35 = ifelse(d35 == 1, family, NA),
          d42 = ifelse(d42 == 1, family, NA)) %>% 
  select(-feature, -family)
 
plotdat <- map2_dfr(dat, colnames(dat), \(x, y) unlist(x) %>% as.character() %>% 
      get_freq(., out_df = T) %>% add_column(group = y)) %>% 
  mutate(group = factor(group, group_order))

p_list <- list()
for (i in unique(plotdat$name)) {
  plotdat_i <- filter(plotdat, name %in% i)
  p_list[[i]] <- ggplot(plotdat_i, aes(x = group, y = freq, fill = group)) +
    geom_bar(stat = "identity", position = position_stack(), show.legend = F, 
             lwd = .4, width = .85, ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = group_color) +
    labs(x = "", y = "", title = i) +
    theme_classic() +
    theme(axis.line = element_line(linewidth = .4, color = "#000000"),
          axis.ticks = element_line(linewidth = .4, color = "#000000"),
          axis.text = element_text(size = 8, color = "#000000"),
          plot.title = element_text(size = 8, color = "#000000", hjust = .5, face = "italic"),
          aspect.ratio = 1)
  }
cowplot::plot_grid(plotlist = p_list, nrow = 8)
ggsave("tmp_taxa_family_barplot.pdf", width = 32, height = 32)


#### cluster analysis p_lvl ####
dat <- read.delim("../../tmp_embryo_dat/bk_p_profile", row.names = 1)
group <- read.delim("../data/sample_group") %>% rbind(., read.delim("../../tmp_embryo_dat/sample_group"))
taxonomy <- read.delim("../../07.gtdbtk/gtdbtk_classify/gtdbtk.bac120.ncbi.summary.tsv") %>% 
  taxa_split(sep = ";", taxonomy_name = c(feature = "Genome.ID", taxonomy = "Majority.vote.NCBI.classification"))
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% data.frame() %>% 
  taxa_trans(taxonomy, group, to = "phylum", out_all = T, smp2grp = F)

rownames(otu) <- gsub("p__", "", x = rownames(otu))
rownames(otu)[rownames(otu)=="Actinobacteria"] = "Actinobacteriota"
rownames(otu)[rownames(otu)=="Bacteroidetes"] = "Bacteroidota"
rownames(otu)[rownames(otu)=="Fusobacteria"] = "Fusobacteriota"
taxa_vec <- rownames(otu)[is.element(rownames(otu), rownames(dat))]

dat <- rownames_to_column(dat, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu <- rownames_to_column(otu, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu2 <- merge(dat, otu, by = "feature") %>% column_to_rownames(var = "feature") %>%
  t() %>% data.frame() %>% mutate(group = group$group[match(rownames(.), group$sample)])

# 层次聚类分析
dist <- dist(otu2[-ncol(otu2)], method = "euclidean")
hclust <- hclust(dist, method = "average")
plot(hclust)

# library(ggtree)
# library(tidytree)
tr <- ape::as.phylo(hclust)
tr_dat <- tr %>% as.treedata() %>% as_tibble()
ggtree(tr, linetype='dashed', color = "#487AA1", layout = "fan")  +
  # geom_tiplab(size = 2) +
  geom_point2(aes(subset=(node==29)), size = 1.2, color = "#fc4e2a") +
  geom_point2(aes(subset=(node==34)), size = 1.2, color = "#fc4e2a") +
  geom_point2(aes(subset=(node==35)), size = 1.2, color = "#fc4e2a") +
  geom_text2(aes(subset=(node==29), label = "d1F"), size = 2, color = "black", nudge_x = .1) +
  geom_text2(aes(subset=(node==34), label = "d1M"), size = 2, color = "black", nudge_x = .1) +
  geom_text2(aes(subset=(node==35), label = "d1N"), size = 2, color = "black", nudge_x = .1) +
  geom_cladelab(node = 121, label = " Embryo microbiome", offset = .02, barsize = .8, 
                size = 1, barcolour = "#fc4e2a", hjust = -.1) +
  geom_cladelab(node = 124, label = " Day 1", offset = .02, barsize = .8, 
                size = 1, barcolour = "#fc4e2a", hjust = -.1)
ggsave("tmp_hclust_tree_p_lvl.pdf", height = 5, width = 5)

#### cluster analysis f_lvl ####
dat <- read.delim("../../tmp_embryo_dat/bk_f_profile", row.names = 1)
group <- read.delim("../data/sample_group") %>% rbind(., read.delim("../../tmp_embryo_dat/sample_group"))
taxonomy <- read.delim("../../07.gtdbtk/gtdbtk_classify/gtdbtk.bac120.ncbi.summary.tsv") %>% 
  taxa_split(sep = ";", taxonomy_name = c(feature = "Genome.ID", taxonomy = "Majority.vote.NCBI.classification"))
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% data.frame() %>% 
  taxa_trans(taxonomy, group, to = "family", out_all = T, smp2grp = F)

rownames(otu) <- gsub("f__", "", x = rownames(otu))
taxa_vec <- rownames(otu)[is.element(rownames(otu), rownames(dat))]

dat <- rownames_to_column(dat, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu <- rownames_to_column(otu, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu2 <- merge(dat, otu, by = "feature") %>% column_to_rownames(var = "feature") %>%
  t() %>% data.frame() %>% mutate(group = group$group[match(rownames(.), group$sample)])

# 层次聚类分析
dist <- dist(otu2[-ncol(otu2)], method = "euclidean")
hclust <- hclust(dist, method = "average")
plot(hclust)

# library(ggtree)
# library(tidytree)
tr <- ape::as.phylo(hclust)
tr_dat <- tr %>% as.treedata() %>% as_tibble()
ggtree(tr, linetype='dashed', color = "#487AA1", layout = "fan")  +
  # geom_tiplab(size = 2) +
  geom_cladelab(node = 127, label = " Embryo microbiome", offset = .02, barsize = .8, 
                size = 1, barcolour = "#fc4e2a", hjust = -.1) +
  geom_cladelab(node = 128, label = " Day 1", offset = .02, barsize = .8, 
                size = 1, barcolour = "#fc4e2a", hjust = -.1)
ggsave("tmp_hclust_tree_f_lvl.pdf", height = 5, width = 5)

#### cluster analysis g_lvl ####
dat <- read.delim("../../tmp_embryo_dat/bk_g_profile", row.names = 1)
group <- read.delim("../data/sample_group") %>% rbind(., read.delim("../../tmp_embryo_dat/sample_group"))
taxonomy <- read.delim("../../07.gtdbtk/gtdbtk_classify/gtdbtk.bac120.ncbi.summary.tsv") %>% 
  taxa_split(sep = ";", taxonomy_name = c(feature = "Genome.ID", taxonomy = "Majority.vote.NCBI.classification"))
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% data.frame() %>% 
  taxa_trans(taxonomy, group, to = "genus", out_all = T, smp2grp = F)

rownames(otu) <- gsub("g__", "", x = rownames(otu))
rownames(dat)[rownames(dat)=="Clostridium sensu stricto 1"] = "Clostridium"
rownames(dat)[rownames(dat)=="Escherichia-Shigella"] = "Escherichia"
taxa_vec <- rownames(otu)[is.element(rownames(otu), rownames(dat))]

dat <- rownames_to_column(dat, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu <- rownames_to_column(otu, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu2 <- merge(dat, otu, by = "feature") %>% column_to_rownames(var = "feature") %>%
  t() %>% data.frame() %>% mutate(group = group$group[match(rownames(.), group$sample)])

# 层次聚类分析
dist <- dist(otu2[-ncol(otu2)], method = "euclidean")
hclust <- hclust(dist, method = "average")
plot(hclust)

# library(ggtree)
# library(tidytree)
tr <- ape::as.phylo(hclust)
tr_dat <- tr %>% as.treedata() %>% as_tibble()
ggtree(tr, linetype='dashed', color = "#487AA1", layout = "fan")  +
  # geom_tiplab(size = 2) +
  geom_point2(aes(subset=(node==32)), size = 1.2, color = "#fc4e2a") +
  geom_point2(aes(subset=(node==34)), size = 1.2, color = "#fc4e2a") +
  geom_text2(aes(subset=(node==32), label = "d1J"), size = 2, color = "black", nudge_x = .1) +
  geom_text2(aes(subset=(node==34), label = "d1M"), size = 2, color = "black", nudge_x = .1) +
  geom_cladelab(node = 125, label = " Embryo", offset = .02, barsize = .8, 
                size = 1, barcolour = "#fc4e2a", hjust = -.1) +
  geom_cladelab(node = 127, label = " Day 1", offset = .02, barsize = .8, 
                size = 1, barcolour = "#fc4e2a", hjust = -.1) +
  geom_cladelab(node = 140, label = " Day 1", offset = .02, barsize = .8, 
                size = 1, barcolour = "#fc4e2a", hjust = -.1) +
  geom_cladelab(node = 153, label = " Day 1", offset = .02, barsize = .8, 
                size = 1, barcolour = "#fc4e2a", hjust = -.1)
ggsave("tmp_hclust_tree_g_lvl.pdf", height = 5, width = 5)

#### corr analysis p_lvl ####
library(ComplexHeatmap)
dat <- read.delim("../../tmp_embryo_dat/bk_p_profile", row.names = 1)
group <- read.delim("../data/sample_group") %>% rbind(., read.delim("../../tmp_embryo_dat/sample_group"))
taxonomy <- read.delim("../../07.gtdbtk/gtdbtk_classify/gtdbtk.bac120.ncbi.summary.tsv") %>% 
  taxa_split(sep = ";", taxonomy_name = c(feature = "Genome.ID", taxonomy = "Majority.vote.NCBI.classification"))
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% data.frame() %>% 
  taxa_trans(taxonomy, group, to = "phylum", out_all = T, smp2grp = F)

rownames(otu) <- gsub("p__", "", x = rownames(otu))
rownames(otu)[rownames(otu)=="Actinobacteria"] = "Actinobacteriota"
rownames(otu)[rownames(otu)=="Bacteroidetes"] = "Bacteroidota"
rownames(otu)[rownames(otu)=="Fusobacteria"] = "Fusobacteriota"
taxa_vec <- rownames(otu)[is.element(rownames(otu), rownames(dat))]

dat <- rownames_to_column(dat, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu <- rownames_to_column(otu, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu2 <- merge(dat, otu, by = "feature") %>% column_to_rownames(var = "feature")

# 相关性分析
library(psych)
corr_res <- corr.test(otu2, method = "spearman")
annotation_col <- data.frame(row.names = colnames(corr_res$r)) %>% 
  mutate(group = group$group[match(rownames(.), group$sample)])
group_color <- structure(c("#f47474","#79b1d3","#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3"),
                 names = c("embryo","d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42"))
color_list <- list(group = group_color)

pdf("tmp_corr_heatmap_p_lvl.pdf", width = 6, height = 6)
pheatmap(corr_res$r, cellwidth = 2, cellheight = 2, show_colnames = F, show_rownames = F,
         annotation_col = annotation_col, annotation_row = annotation_col, 
         treeheight_row = 10, treeheight_col = 10, annotation_names_row = F,
         annotation_colors = color_list,
         color = colorRampPalette(rev(c("#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6")))(100))
dev.off()

#### corr analysis f_lvl ####
library(ComplexHeatmap)
dat <- read.delim("../../tmp_embryo_dat/bk_f_profile", row.names = 1)
group <- read.delim("../data/sample_group") %>% rbind(., read.delim("../../tmp_embryo_dat/sample_group"))
taxonomy <- read.delim("../../07.gtdbtk/gtdbtk_classify/gtdbtk.bac120.ncbi.summary.tsv") %>% 
  taxa_split(sep = ";", taxonomy_name = c(feature = "Genome.ID", taxonomy = "Majority.vote.NCBI.classification"))
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% data.frame() %>% 
  taxa_trans(taxonomy, group, to = "family", out_all = T, smp2grp = F)

rownames(otu) <- gsub("f__", "", x = rownames(otu))
taxa_vec <- rownames(otu)[is.element(rownames(otu), rownames(dat))]

dat <- rownames_to_column(dat, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu <- rownames_to_column(otu, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu2 <- merge(dat, otu, by = "feature") %>% column_to_rownames(var = "feature")

# 相关性分析
library(psych)
corr_res <- corr.test(otu2, method = "spearman")
annotation_col <- data.frame(row.names = colnames(corr_res$r)) %>% 
  mutate(group = group$group[match(rownames(.), group$sample)])
group_color <- structure(c("#f47474","#79b1d3","#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3"),
                         names = c("embryo","d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42"))
color_list <- list(group = group_color)

pdf("tmp_corr_heatmap_f_lvl.pdf", width = 6, height = 6)
pheatmap(corr_res$r, cellwidth = 2, cellheight = 2, show_colnames = F, show_rownames = F,
         annotation_col = annotation_col, annotation_row = annotation_col, 
         treeheight_row = 10, treeheight_col = 10, annotation_names_row = F,
         annotation_colors = color_list,
         color = colorRampPalette(rev(c("#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6")))(100))
dev.off()

#### corr analysis g_lvl ####
library(ComplexHeatmap)
dat <- read.delim("../../tmp_embryo_dat/bk_g_profile", row.names = 1)
group <- read.delim("../data/sample_group") %>% rbind(., read.delim("../../tmp_embryo_dat/sample_group"))
taxonomy <- read.delim("../../07.gtdbtk/gtdbtk_classify/gtdbtk.bac120.ncbi.summary.tsv") %>% 
  taxa_split(sep = ";", taxonomy_name = c(feature = "Genome.ID", taxonomy = "Majority.vote.NCBI.classification"))
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% data.frame() %>% 
  taxa_trans(taxonomy, group, to = "genus", out_all = T, smp2grp = F)

rownames(otu) <- gsub("g__", "", x = rownames(otu))
rownames(dat)[rownames(dat)=="Clostridium sensu stricto 1"] = "Clostridium"
rownames(dat)[rownames(dat)=="Escherichia-Shigella"] = "Escherichia"
taxa_vec <- rownames(otu)[is.element(rownames(otu), rownames(dat))]

dat <- rownames_to_column(dat, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu <- rownames_to_column(otu, var = "feature") %>% 
  mutate(feature = ifelse(feature %in% taxa_vec, feature, "Other")) %>% group_by(feature) %>% 
  summarise_all(sum) %>% ungroup()
otu2 <- merge(dat, otu, by = "feature") %>% column_to_rownames(var = "feature")

# 相关性分析
library(psych)
corr_res <- corr.test(otu2, method = "spearman")
annotation_col <- data.frame(row.names = colnames(corr_res$r)) %>% 
  mutate(group = group$group[match(rownames(.), group$sample)])
group_color <- structure(c("#f47474","#79b1d3","#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f","#e5c494","#b3b3b3"),
                         names = c("embryo","d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42"))
color_list <- list(group = group_color)

pdf("tmp_corr_heatmap_g_lvl.pdf", width = 6, height = 6)
pheatmap(corr_res$r, cellwidth = 2, cellheight = 2, show_colnames = F, show_rownames = F,
         annotation_col = annotation_col, annotation_row = annotation_col, 
         treeheight_row = 10, treeheight_col = 10, annotation_names_row = F,
         annotation_colors = color_list,
         color = colorRampPalette(rev(c("#d7191c","#fdae61","#ffffbf","#abd9e9","#2c7bb6")))(100)) 
dev.off()

