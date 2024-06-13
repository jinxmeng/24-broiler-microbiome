#### info ####
# Jinxin Meng, 2023-10-09, 2023-11-11
pacman::p_load(tibble, dplyr, tidyr, ggplot2, vegan, ggpmisc, ggpubr)
setwd("F:/Microbiome/proj_2023/02.Broiler_fecal_resistome_20230302/taxa/")
source("F:/Code/R_func/taxa.R")

#### phylum composition ####
group <- read.delim("../data/sample_group")
otu <- read.delim("../data/species_tpm", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% 
  data.frame()
taxonomy <- read.delim("../data/species_taxonomy_rename.tsv")

dat <- taxa_trans(otu, taxonomy, group, to = "phylum", out_all = T) %>% 
  rownames_to_column(var = "feature")

taxa_order <- c("p__Bacillota", "p__Bacteroidota", "p__Pseudomonadota", "p__Actinomycetota", "p__Cyanobacteriota", "p__Desulfobacterota", "p__Fusobacteriota")
taxa_color <- c("#80b1d3", "#b3de69", "#fdb462", "#8dd3c7", "#bc80bd", "#fb8072", "#ffed6f")
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2")

plotdat <- dat %>% 
  gather(., key = "group", val = "val", -feature) %>% 
  mutate(group = factor(group, group_order),
         feature = factor(feature, rev(taxa_order)))

# phylum composition xlim 90-100 | d3-d42
p1 <- ggplot(plotdat %>% filter(group != "d1"), aes(val, group, fill = feature)) +
  geom_bar(stat = "identity", position = position_stack(), color = "#000000", size = .2, width = .8) +
  scale_fill_manual(values = rev(taxa_color)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.05), labels = seq(0, 100, by = 5),
                     position = "top") +
  labs(x = "Relative Abundance(%)", y = "", fill = "Phylum taxa") +
  coord_cartesian(x = c(0.9, 1)) +
  theme_classic() +
  theme(axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        axis.text = element_text(size = 8, color = "#000000"),
        axis.title = element_text(size = 8, color = "#000000"),
        plot.title = element_text(size = 10, color = "#000000"),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        legend.title = element_text(size = 10, color = "#000000"),
        panel.grid = element_blank())
# phylum composition xlim 50-100 | d1
p2 <- ggplot(plotdat %>% filter(group == "d1"), aes(val, group, fill = feature)) +
  geom_bar(stat = "identity", position = position_stack(), color = "#000000", size = .2, width = .8) +
  scale_fill_manual(values = rev(taxa_color)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.1), labels = seq(0, 100, by = 10),
                     position = "bottom") +
  labs(x = "Relative Abundance(%)", y = "", fill = "Phylum taxa") +
  coord_cartesian(x = c(0.5, 1)) +
  theme_classic() +
  theme(axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        axis.text = element_text(size = 8, color = "#000000"),
        axis.title = element_text(size = 8, color = "#000000"),
        plot.title = element_text(size = 10, color = "#000000"),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        legend.title = element_text(size = 10, color = "#000000"),
        panel.grid = element_blank())
p1 / p2 + patchwork::plot_layout(heights = c(8, 1))
ggsave("compos_phylum_barplot.pdf", height = 4, width = 6)

# 保存数据
write.table(dat, "compos_phylum_dat.txt", sep = "\t", quote = F, row.names = F)
write_lines("\n", "compos_phylum_dat.txt", append = T)
dat <- taxa_trans(otu, taxonomy, to = "phylum", out_all = T, smp2grp = F)
otu_diff_test(dat, group, group_order, method = "wilcox") %>% 
  mutate(plab = ifelse(pval < 0.0001, "****", 
                       ifelse(pval < 0.001, "***", 
                              ifelse(pval < 0.01, "**", 
                                     ifelse(pval < 0.05, "*", "!"))))) %>% 
  write.table(., "compos_phylum_dat.txt", sep = "\t", quote = F, row.names = F, append = T)

#### F/B ####
plotdat <- taxa_trans(otu, taxonomy, to = "phylum", out_all = T, smp2grp = F) %>% 
  t() %>% 
  data.frame(check.names = F) %>% 
  rownames_to_column(var = "sample") %>% 
  select(sample, p__Bacillota, p__Bacteroidota) %>% 
  mutate(val = p__Bacillota / (p__Bacteroidota + 1e-5)) %>% 
  merge(., group, by = "sample") %>% 
  select(sample, group, val) %>% 
  mutate(group = factor(group, group_order),
         val = log10(val + 1))

comparison <- diff_test(select(plotdat, sample, val), group, method = "wilcox", 
                        group_order = group_order) %>% 
  filter(pval < 0.05) %>% 
  select(group_pair) %>% 
  unlist() %>% 
  as.character() %>% 
  strsplit(x = ., split = "_vs_")

plotdat_point <- plotdat %>% 
  group_by(group) %>% 
  summarise(mean = mean(val), sd = sd(val))

ggplot(plotdat, aes(group, val)) +
  geom_hline(yintercept = 1, lwd = .4, lty = "dashed", color = "#000000") +
  geom_errorbar(data = plotdat_point, aes(group, ymin = mean - sd*.8, ymax = mean + sd*.8),
                width = .3, lwd = .4, inherit.aes = F, show.legend = F) +
  geom_jitter(aes(color = group), width = .2, size = 1.5, show.legend = F, height = .12) +
  stat_smooth(aes(group = 1), level = 0.95, method = "lm",
              color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 2, label.x.npc = .5, label.y.npc = 0, color = "#8da0cb") +
  stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5","d7")), aes(group = 1), 
              level = 0.95, method = "lm", color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5","d7")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 2, label.x.npc = .5, label.y.npc = .03, color = "#fc8d62") +
  stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5","d7")), aes(group = 1), 
              level = 0.95, method = "lm", color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5","d7")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 2, label.x.npc = .5, label.y.npc = .06, color = "#66c2a5") +
  geom_point(data = plotdat_point, aes(group, mean), size = 2.8, shape = 21, stroke = .4,
             fill = "white", color = "black", inherit.aes = F) +
  stat_compare_means(comparisons = comparison, method = "wilcox", method.args = list(exact = F), size = 3,
                     label = "p.signif", tip.length = .01, step.increase = .03, vjust = .9) +
  scale_color_manual(values = group_color) +
  scale_y_continuous(expand = c(.03, .03)) +
  labs(x = "", y = "log10(Overall_Bacillota / Bacteroidota ratio)") +
  theme_classic() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 10, color = "black"),
        aspect.ratio = 4/3)
ggsave("compos_FvsB_boxplot.pdf", width = 3, height = 4)


#### phylum change with lm fitted ####
dat <- taxa_trans(otu, taxonomy, to = "phylum", out_all = T, smp2grp = F) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "sample")

plotlist <- list()
for (i in colnames(dat)[-1]) {
  plotdat <- select(dat, sample, val = all_of(i)) %>% 
    merge(., group, by = "sample") %>% 
    mutate(group = factor(group, group_order),
           val = log10(100 * val))
  
  comparison <- diff_test(select(plotdat, -group), group, method = "wilcox", 
                          group_order = group_order, paired = F) %>% 
    filter(pval < 0.05) %>% 
    select(group_pair) %>% 
    unlist() %>% 
    as.character() %>% 
    strsplit(x = ., split = "_vs_")
  
  p <- ggplot(plotdat, aes(group, val)) +
    stat_compare_means(comparisons = comparison, method = "wilcox", method.args = list(exact = F), lwd = .4,
                       label = "p.signif", tip.length = .01, step.increase = .03, vjust = .95, size = 2.5) +
    geom_boxplot(width = .5, fill = "transparent", lwd = .4, outlier.shape = NA, show.legend = F) +
    geom_jitter(aes(x = group, y = val, color = group), width = .3, size = .8, inherit.aes = F, show.legend = F) +
    stat_smooth(aes(group = 1), level = 0.95, method = "lm",
                color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
    stat_poly_eq(aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
                 size = 1, label.x.npc = 0, label.y.npc = 1, color = "#8da0cb") +
    stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5","d7")), aes(group = 1), 
                level = 0.95, method = "lm", color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
    stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5","d7")), 
                 aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
                 size = 1, label.x.npc = .5, label.y.npc = 1, color = "#fc8d62") +
    stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5","d7")), aes(group = 1), 
                level = 0.95, method = "lm", color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
    stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5","d7")), 
                 aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
                 size = 1, label.x.npc = 1, label.y.npc = 1, color = "#66c2a5") +
    scale_color_manual(values = group_color) +
    labs(x = "", title = i, y = "log10(Relative Abundance %") +
    theme_classic() +
    theme(axis.line = element_line(linewidth = .4, color = "#000000"),
          axis.ticks = element_line(linewidth = .4, color = "#000000"),
          axis.text = element_text(size = 8, color = "#000000"),
          axis.title = element_text(size = 8, color = "#000000"),
          plot.title = element_text(size = 10, face = "italic", color = "#000000"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          aspect.ratio = 1)
    plotlist[[i]] <- p
}
cowplot::plot_grid(plotlist = plotlist, nrow = 2)
ggsave("compos_phylum_boxplot.pdf", width = 12, height = 6)


#### family composition ####
dat <- taxa_trans(otu, taxonomy, group, to = "family", top_n = 20, other_name = "f__Other")

taxa_order <- data.frame(x = rowSums(dat)) %>% arrange(desc(x)) %>% rownames()
taxa_color <- paletteer::paletteer_d("ggthemes::Tableau_20")

plotdat <- dat %>%
  rownames_to_column(var = "feature") %>% 
  gather(., key = "group", val = "val", -feature) %>% 
  mutate(group = factor(group, group_order),
         feature = factor(feature, rev(taxa_order)))

# phylum composition xlim 70-100 | d5-d42
p1 <- ggplot(plotdat %>% filter(!group %in% c("d1", "d3")), aes(val, group, fill = feature)) +
  geom_bar(stat = "identity", position = position_stack(), color = "#000000", size = .2, width = .8) +
  scale_fill_manual(values = rev(taxa_color)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.05), labels = seq(0, 100, by = 5),
                     position = "top") +
  labs(x = "Relative Abundance(%)", y = "", fill = "Phylum taxa") +
  coord_cartesian(x = c(.7, 1)) +
  theme_classic() +
  theme(axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        axis.text = element_text(size = 8, color = "#000000"),
        axis.title = element_text(size = 8, color = "#000000"),
        plot.title = element_text(size = 10, color = "#000000"),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        legend.title = element_text(size = 10, color = "#000000"),
        panel.grid = element_blank())
# phylum composition xlim 0-100 | d1-d3
p2 <- ggplot(plotdat %>% filter(group %in% c("d1", "d3")), aes(val, group, fill = feature)) +
  geom_bar(stat = "identity", position = position_stack(), color = "#000000", size = .2, width = .8) +
  scale_fill_manual(values = rev(taxa_color)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0, 1, by = 0.1), labels = seq(0, 100, by = 10),
                     position = "bottom") +
  labs(x = "Relative Abundance(%)", y = "", fill = "Phylum taxa") +
  theme_classic() +
  theme(axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        axis.text = element_text(size = 8, color = "#000000"),
        axis.title = element_text(size = 8, color = "#000000"),
        plot.title = element_text(size = 10, color = "#000000"),
        legend.text = element_text(size = 8, color = "#000000", face = "italic"),
        legend.title = element_text(size = 10, color = "#000000"),
        panel.grid = element_blank())
p1 / p2 + patchwork::plot_layout(heights = c(7, 2))
ggsave("compos_family_barplot.pdf", height = 4, width = 6)

# 保存数据
write.table(dat, "compos_family_dat.txt", sep = "\t", quote = F, row.names = F)
write_lines("\n", "compos_family_dat.txt", append = T)
dat <- taxa_trans(otu, taxonomy, to = "family", out_all = T, smp2grp = F)
otu_diff_test(dat, group, group_order, method = "wilcox") %>% 
  mutate(plab = ifelse(pval < 0.0001, "****", 
                       ifelse(pval < 0.001, "***", 
                              ifelse(pval < 0.01, "**", 
                                     ifelse(pval < 0.05, "*", "!"))))) %>% 
  write.table(., "compos_family_dat.txt", sep = "\t", quote = F, row.names = F, append = T)

#### family change ####
dat <- taxa_trans(otu, taxonomy, to = "family", out_all = T, smp2grp = F) %>% 
  t() %>% 
  data.frame() %>% 
  rownames_to_column(var = "sample")

plotlist <- list()
for (i in colnames(dat)[-1]) {
  plotdat <- select(dat, sample, val = all_of(i)) %>% 
    merge(., group, by = "sample") %>% 
    mutate(group = factor(group, group_order),
           val = log10(100 * val))
  
  comparison <- diff_test(select(plotdat, -group), group, method = "wilcox", 
                          group_order = group_order, paired = F) %>% 
    filter(pval < 0.05) %>% 
    select(group_pair) %>% 
    unlist() %>% 
    as.character() %>% 
    strsplit(x = ., split = "_vs_")
  
  p <- ggplot(plotdat, aes(group, val)) +
    stat_compare_means(comparisons = comparison, method = "wilcox", method.args = list(exact = F), lwd = .4,
                       label = "p.signif", tip.length = .01, step.increase = .03, vjust = .95, size = 2.5) +
    geom_boxplot(width = .5, fill = "transparent", lwd = .4, outlier.shape = NA, show.legend = F) +
    geom_jitter(aes(x = group, y = val, color = group), width = .3,
                size = .8, inherit.aes = F, show.legend = F) +
    scale_color_manual(values = group_color) +
    labs(x = "", title = i, y = "log10(Relative Abundance %") +
    theme_classic() +
    theme(axis.line = element_line(linewidth = .4, color = "#000000"),
          axis.ticks = element_line(linewidth = .4, color = "#000000"),
          axis.text = element_text(size = 8, color = "#000000"),
          axis.title = element_text(size = 8, color = "#000000"),
          plot.title = element_text(size = 10, face = "italic", color = "#000000"),
          panel.background = element_blank(),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          aspect.ratio = 1)
  plotlist[[i]] <- p
}
cowplot::plot_grid(plotlist = plotlist, nrow = 8)
ggsave("compos_family_boxplot.pdf", width = 24, height = 24, limitsize = F)

dat <- taxa.compos.tib(otu,  group = group, taxonomy = taxonomy,
                       taxa_level = "family", display = "sample", show_taxa = 100) %>%
  select(-group)

# 初步查看每个物种的情况
plotlist <- list()
for (i in colnames(dat)[-1]) {
  plotdat <- dat %>% 
    select(sample, all_of(i)) %>% 
    merge(., group, by = "sample") %>% 
    rename(val = all_of(setdiff(colnames(.), c("sample", "group")))) %>% 
    select(-sample) %>% 
    mutate(group = factor(group, levels = levels),
           val = log10(val*100 + 1))
  
  p <- ggplot(plotdat, aes(group, val)) +
    geom_boxplot(fill = "transparent", width = .5, outlier.shape = NA) +
    geom_jitter(aes(color = group), width = .3, size = 1.8, show.legend = F) +
    stat_smooth(aes(group = 1), span = 2, fill = "grey75", color = "grey50", method = "lm", lwd = .8) +
    scale_color_manual(values = colors) +
    labs(x = "", title = i, y = expression(log[10](Relative~Abundance~'%'))) +
    theme_classic() +
    theme(axis.line = element_line(linewidth = .4, color = "#000000"),
          axis.ticks = element_line(linewidth = .4, color = "#000000"),
          axis.text = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 8, color = "black"),
          plot.title = element_text(size = 10, color = "black", face = "italic"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          aspect.ratio = 1)
  plotlist[[i]] <- p
}
cowplot::plot_grid(plotlist = plotlist, nrow = 10)
ggsave("compos_family_boxplot.pdf", width = 21, height = 30)


#### genus composition ####
dat <- taxa_trans(otu, taxonomy, group, to = "genus", out_all = T) %>% 
  rownames_to_column(var = "feature")

write.table(dat, "compos_genus_dat.txt", sep = "\t", quote = F, row.names = F)
write_lines("\n", "compos_genus_dat.txt", append = T)
dat <- taxa_trans(otu, taxonomy, to = "genus", out_all = T, smp2grp = F)
otu_diff_test(dat, group, group_order, method = "wilcox") %>% 
  mutate(plab = ifelse(pval < 0.0001, "****", 
                       ifelse(pval < 0.001, "***", 
                              ifelse(pval < 0.01, "**", 
                                     ifelse(pval < 0.05, "*", "!"))))) %>% 
  write.table(., "compos_genus_dat.txt", sep = "\t", quote = F, row.names = F, append = T)

#### species composition ####
dat <- taxa_trans(otu, taxonomy, group, to = "species", out_all = T) %>% 
  rownames_to_column(var = "feature")

write.table(dat, "compos_species_dat.txt", sep = "\t", quote = F, row.names = F)
write_lines("\n", "compos_species_dat.txt", append = T)
dat <- taxa_trans(otu, taxonomy, to = "species", out_all = T, smp2grp = F)
otu_diff_test(dat, group, group_order, method = "wilcox") %>% 
  mutate(plab = ifelse(pval < 0.0001, "****", 
                       ifelse(pval < 0.001, "***", 
                              ifelse(pval < 0.01, "**", 
                                     ifelse(pval < 0.05, "*", "!"))))) %>% 
  write.table(., "compos_species_dat.txt", sep = "\t", quote = F, row.names = F, append = T)

#### species corr ####
source("F:/Code/R_func/taxa.R")
group <- read.delim("../data/sample_group.txt")
taxonomy <- read.delim("species_taxonomy_rename.tsv") 
profile <- read.delim("species_tpm_filter", row.names = 1)%>% taxa_trans(taxonomy, to = "species", out_all = T) %>% t() %>% data.frame()
x = cor(profile, method = "spearman")

pdf("tmp_20231226_species_corr_heatmap.pdf", width = 25, height = 25)
pheatmap(x, cellwidth = 4, cellheight = 4, show_colnames = F, fontsize_row = 3, treeheight_row = 150, 
         treeheight_col = 150, cutree_rows = 10, cutree_cols = 10, border_gp = gpar(col = "black"))
dev.off()

#### survey species with family ####
library(ComplexHeatmap)
source("F:/Code/R_func/taxa.R")
group <- read.delim("../data/sample_group.txt")
taxonomy <- read.delim("species_taxonomy_rename.tsv") 
profile <- read.delim("species_tpm_filter", row.names = 1)%>% taxa_trans(taxonomy, to = "species", out_all = T)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

annotation_row <- data.frame(name = rownames(profile)) %>% 
  mutate(phylum = taxonomy$phylum[match(name, taxonomy$species)],
         phylum = forcats::fct_lump_n(phylum, 10),
         family = taxonomy$family[match(name, taxonomy$species)],
         family = forcats::fct_lump_n(family, 19, ties.method = "last")) %>% 
  column_to_rownames(var = "name")
phylum_color <- structure(c("#88cbc0","#7dabc9","#add26a","#b27cb3","#ee7c6e","#f5e571","#f4af63"), names = levels(annotation_row$phylum))
family_color <- structure(paletteer::paletteer_d("ggthemes::Tableau_20") %>% as.character(), names = levels(annotation_row$family))
annotation_col <- data.frame(name = colnames(profile)) %>% 
  mutate(group = group$group[match(name, group$sample)]) %>% 
  column_to_rownames(var = "name")
group_color <- structure(group_color, names = group_order)
annotation_colors <- list(phylum = phylum_color, family = family_color, group = group_color)
col_split <- factor(annotation_col$group, group_order)
row_split <- factor(annotation_row$family)

pdf("tmp_20231226_species_with_family_change_pheatmap_absence_presence.pdf", width = 12, height = 20)
pheatmap(profile, color = circlize::colorRamp2(c(0, 1), c("white","#74add1")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, row_split = row_split, 
         column_gap = unit(1,'mm'), border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 3, fontsize_col = 0,
         cellwidth = 4, cellheight = 4, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

pdf("tmp_20231226_species_with_family_change_pheatmap_abund.pdf", width = 12, height = 20)
pheatmap(log10(profile+1e-3), color = circlize::colorRamp2(c(0, 2), c("white","#74add1")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, row_split = row_split,
         column_gap = unit(1,'mm'), border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 2, fontsize_col = 0,
         cellwidth = 4, cellheight = 4, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()
#### survey species with genus ####
library(ComplexHeatmap)
source("F:/Code/R_func/taxa.R")
group <- read.delim("../data/sample_group.txt")
taxonomy <- read.delim("species_taxonomy_rename.tsv") 
profile <- read.delim("species_tpm_filter", row.names = 1)%>% taxa_trans(taxonomy, to = "species", out_all = T)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

annotation_row <- data.frame(name = rownames(profile)) %>% 
  mutate(phylum = taxonomy$phylum[match(name, taxonomy$species)],
         phylum = forcats::fct_lump_n(phylum, 10),
         genus = taxonomy$genus[match(name, taxonomy$species)]) %>% 
  column_to_rownames(var = "name")
phylum_color <- structure(c("#88cbc0","#7dabc9","#add26a","#b27cb3","#ee7c6e","#f5e571","#f4af63"), names = levels(annotation_row$phylum))
annotation_col <- data.frame(name = colnames(profile)) %>% 
  mutate(group = group$group[match(name, group$sample)]) %>% 
  column_to_rownames(var = "name")
group_color <- structure(group_color, names = group_order)
annotation_colors <- list(phylum = phylum_color, group = group_color)
col_split <- factor(annotation_col$group, group_order)
row_split <- factor(annotation_row$genus)

pdf("tmp_20231226_species_with_genus_change_pheatmap_absence_presence.pdf", width = 12, height = 32)
pheatmap(profile, color = circlize::colorRamp2(c(0, 1), c("white","#74add1")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, row_split = row_split, row_title_rot = 0,
         column_gap = unit(1,'mm'), border_gp = gpar(col = "black"), row_title_gp = gpar(fontsize = 7),
         show_rownames = T, show_colnames = F, fontsize_row = 3, fontsize_col = 0,
         cellwidth = 6, cellheight = 6, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

pdf("tmp_20231226_species_with_genus_change_pheatmap_abund.pdf", width = 18, height = 32)
pheatmap(log10(profile+1e-3), color = circlize::colorRamp2(c(0, 2), c("white","#74add1")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, row_split = row_split, row_title_rot = 0,
         column_gap = unit(1,'mm'), border_gp = gpar(col = "black"), row_title_gp = gpar(fontsize = 7),
         show_rownames = T, show_colnames = F, fontsize_row = 5, fontsize_col = 0,
         cellwidth = 6, cellheight = 6, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

#### survey genus ####
source("F:/Code/R_func/taxa.R")
group <- read.delim("../data/sample_group.txt")
taxonomy <- read.delim("species_taxonomy_rename.tsv")
profile <- read.delim("species_tpm_filter", row.names = 1) %>% taxa_trans(taxonomy, to = "genus", out_all = T)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

annotation_row <- data.frame(name = rownames(profile)) %>% 
  mutate(phylum = taxonomy$phylum[match(name, taxonomy$genus)],
         phylum = forcats::fct_lump_n(phylum, 10),
         family = taxonomy$family[match(name, taxonomy$genus)],
         family = forcats::fct_lump_n(family, 19, ties.method = "last")) %>% 
  column_to_rownames(var = "name")
phylum_color <- structure(c("#88cbc0","#7dabc9","#add26a","#b27cb3","#ee7c6e","#f5e571","#f4af63"), names = levels(annotation_row$phylum))
family_color <- structure(paletteer::paletteer_d("ggthemes::Tableau_20") %>% as.character(), names = levels(annotation_row$family))
annotation_col <- data.frame(name = colnames(profile)) %>% 
  mutate(group = group$group[match(name, group$sample)]) %>% 
  column_to_rownames(var = "name")
group_color <- structure(group_color, names = group_order)
annotation_colors <- list(phylum = phylum_color, family = family_color, group = group_color)
col_split <- factor(annotation_col$group, group_order)
row_split <- factor(annotation_row$family)

pdf("tmp_20231226_genus_change_pheatmap_absence_presence.pdf", width = 12, height = 12)
pheatmap(profile, color = circlize::colorRamp2(c(0, 1), c("white","#74add1")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, row_split = row_split, 
         column_gap = unit(1,'mm'), border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 3, fontsize_col = 0,
         cellwidth = 4, cellheight = 4, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

pdf("tmp_20231226_genus_change_pheatmap_abund.pdf", width = 12, height = 12)
pheatmap(log10(profile+1e-3), color = circlize::colorRamp2(c(log10(1e-3), 3), c("white","#74add1")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, row_split = row_split, 
         column_gap = unit(1,'mm'), border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 3, fontsize_col = 0,
         cellwidth = 5, cellheight = 4, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

#### survey family ####
source("F:/Code/R_func/taxa.R")
group <- read.delim("../data/sample_group.txt")
taxonomy <- read.delim("species_taxonomy_rename.tsv")
profile <- read.delim("species_tpm_filter", row.names = 1) %>% taxa_trans(taxonomy, to = "family", out_all = T)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

annotation_row <- data.frame(name = rownames(profile)) %>% 
  mutate(phylum = taxonomy$phylum[match(name, taxonomy$family)],
         phylum = forcats::fct_lump_n(phylum, 10)) %>% 
  column_to_rownames(var = "name")
phylum_color <- structure(c("#88cbc0","#7dabc9","#add26a","#b27cb3","#ee7c6e","#f5e571","#f4af63"), names = levels(annotation_row$phylum))
annotation_col <- data.frame(name = colnames(profile)) %>% 
  mutate(group = group$group[match(name, group$sample)]) %>% 
  column_to_rownames(var = "name")
group_color <- structure(group_color, names = group_order)
annotation_colors <- list(phylum = phylum_color, group = group_color)
col_split <- factor(annotation_col$group, group_order)

pdf("tmp_20231226_family_change_pheatmap_absence_presence.pdf", width = 12, height = 8)
pheatmap(profile, color = circlize::colorRamp2(c(0, 1), c("white","#74add1")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, column_gap = unit(1,'mm'), border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 4, fontsize_col = 0,
         cellwidth = 5, cellheight = 5, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

pdf("tmp_20231226_family_change_pheatmap_abund.pdf", width = 12, height = 8)
pheatmap(log10(profile+1e-3), color = circlize::colorRamp2(c(log10(1e-3), 3), c("white","#74add1")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, column_gap = unit(1,'mm'), border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 4, fontsize_col = 0,
         cellwidth = 5, cellheight = 5, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()


#### mffuzz 时间动态变化分型 ####
library(Mfuzz)
source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/profile_process.R")
group <- read.delim("../data/sample_group.txt")
taxonomy <- read.delim("species_taxonomy_rename.tsv") %>% 
  mutate(species =  paste0(species," ",feature))
profile <- read.delim("species_tpm_filter", row.names = 1) %>% taxa_trans(taxonomy, to = "species", out_all = T, transRA = T)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

# 构造数据
dat <- profile_smp2grp(profile, group, method = "mean") %>% 
  relocate(all_of(group$group)) %>% as.matrix()

# 构造mfuzz对象
mfuzz <- new("ExpressionSet", exprs = dat)

# 预处理缺失值或者异常值
mfuzz <- filter.NA(mfuzz, thres = 0.25)
mfuzz <- fill.NA(mfuzz, mode = "mean")
mfuzz <- filter.std(mfuzz, min.std = 0)

# 标准化数据
mfuzz <- standardise(mfuzz)

# Mfuzz基于fuzzy c-means的算法进行聚类
# 需手动定义目标聚类群的个数，例如这里我们为了重现原作者的结果，设定为 10，即期望获得 10 组聚类群
# 需要设定随机数种子，以避免再次运行时获得不同的结果
set.seed(2024)
n = 12
cluster <- mfuzz(mfuzz, c = n, m = mestimate(mfuzz))

# 提取所有蛋白所属的聚类群，并和它们的原始表达值整合在一起
cluster_dat <- data.frame(cluster =  cluster$cluster) %>% rownames_to_column(var = "taxa")

# 作图，详情 ?mfuzz.plot2
# time.labels 参数设置时间轴，需要和原基因表达数据集中的列对应
# 颜色、线宽、坐标轴、字体等细节也可以添加其他参数调整，此处略，详见函数帮助
# mfuzz.plot2(mfuzz, cl = cluster, mfrow = c(3, 5), time.labels = colnames(dat),
#             xlab = "", ylab = "Scale(Relative Abundance %)")
# 查看每个聚类群中各自包含的蛋白数量
# data.frame(cluster = paste0("Cluster", 1:n), cluster_n = cluster$size)
# Mfuzz通过计算一个叫membership的统计量判断蛋白质所属的聚类群，以最大的membership值为准
# 查看各蛋白的membership值
# head(cluster$membership)

# 合并profile表
# dat <- merge(data.frame(dat) %>% rownames_to_column(var = "taxa"), cluster_dat, by = "taxa")
# write.table(dat, "marker_cluster.txt", sep = "\t", row.names = F, quote = F)

# 如果您想提取数据分析过程中，标准化后的表达值（绘制曲线图用的那个值，而非原始蛋白表达值）
dat <- merge(mfuzz@assayData$exprs %>% data.frame(check.names = F) %>% rownames_to_column(var = "taxa"), cluster_dat, by = "taxa")
# write.table(dat2, "marker_cluster_scale.txt", sep = "\t", row.names = F, quote = F)

p_lst <- list()
for (i in 1:n) {
  plotdat <- filter(dat, cluster == i) %>% select(-cluster) %>% 
    gather(., key = "group", value = "val", -taxa) %>% 
    mutate(group = factor(group, group_order))
  p <- ggplot(plotdat, aes(x = group, y = val, color = taxa, group = taxa)) +
    geom_smooth(aes(group = taxa), method = "loess", se = F, linewidth = .8, show.legend = F) +
    scale_x_discrete(expand = c(.01, .01)) +
    labs(x = "", y = "Scale (Relative Abundance %)", title = paste0("Cluster", i)) +
    theme_classic() + 
    theme(axis.text = element_text(color = "black", size = 10),
          axis.title = element_text(size = 12),
          axis.ticks = element_line(color = "black", size = .5),
          plot.title = element_text(hjust = .5),
          panel.grid = element_blank(), aspect.ratio = 1)
  p_lst[[i]] <- p
}
cowplot::plot_grid(plotlist = p_lst, nrow = 2, ncol = 6)
ggsave("tmp_20231226_mufzz_lineplot.pdf", width = 18, height = 6)

# 每个cls的物种数量
plotdat <- merge(select(dat, taxa, cluster),
                 select(taxonomy, species, family, genus) , by.x = "taxa", by.y = "species")
family_color <- paletteer::paletteer_d("ggthemes::Tableau_20")
plotdat %>% mutate(family = forcats::fct_lump_n(family, 19, ties.method = "last", other_level = "f__Other")) %>% 
  group_by(cluster, family) %>% summarise(n = n()) %>% ungroup %>% mutate(cluster = as.character(cluster)) %>% 
  ggbarplot(x = "cluster", y = "n", fill = "family", ylab = "Number of species", lwd = .4) +
  scale_fill_manual(values = family_color) +
  scale_y_continuous(expand = c(0, 0))
ggsave("tmp_20231226_mufzz_taxa_family.pdf", width = 6, height = 6)


# 毛螺菌
plotdat <- merge(select(dat, taxa, cluster),
                 select(taxonomy, species, family, genus) , by.x = "taxa", by.y = "species") %>% 
  filter(family == "f__Ruminococcaceae")
family_color <- paletteer::paletteer_d("ggthemes::Tableau_20")
plotdat %>% mutate(genus = forcats::fct_lump_n(genus, 19, ties.method = "last", other_level = "f__Other")) %>% 
  group_by(cluster, genus) %>% summarise(n = n()) %>% ungroup %>% mutate(cluster = as.character(cluster)) %>% 
  ggbarplot(x = "cluster", y = "n", fill = "genus", ylab = "Number of species", lwd = .4) +
  scale_fill_manual(values = family_color) +
  scale_y_continuous(expand = c(0, 0))
ggsave("tmp_20231226_mufzz_taxa_f__Ruminococcaceae.pdf", width = 3, height = 5)



#### WGCNA ####
install.packages("WGCNA") 
source("F:/Code/R_func/taxa.R")
group <- read.delim("../data/sample_group.txt")
taxonomy <- read.delim("species_taxonomy_rename.tsv") 
profile <- read.delim("species_tpm_filter", row.names = 1)%>% taxa_trans(taxonomy, to = "species", out_all = T) %>% t() %>% data.frame()
x = cor(profile, method = "spearman")

pdf("tmp_20231226_species_corr_heatmap.pdf", width = 25, height = 25)
pheatmap(x, cellwidth = 4, cellheight = 4, show_colnames = F, fontsize_row = 3, treeheight_row = 150, 
         treeheight_col = 150, cutree_rows = 10, cutree_cols = 10, border_gp = gpar(col = "black"))
dev.off()