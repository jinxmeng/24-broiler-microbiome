#### info ####
# Jinxin Meng, 20230703, 20231207
pacman::p_load(tidyr, dplyr, tibble, purrr)
pacman::p_load(ggplot2, ggrepel, ggforce, ggpubr, ggpmisc)
setwd("F:/Microbiome/proj_2023/02.Broiler_fecal_resistome_20230302/ARGs/")

#### ARGs drug statistics ####
source("F:/Code/R_func/plot_pie.R")
dat <- read.delim("ARGs_drug_manual_classify.txt")
plot_pie(dat, dat_colnames = c(name = "drug"), top_n = 12, fill = "auto")
ggsave("ARGs_drug_manual_classify_pie.pdf", width = 3, height = 3)

dat <- read.delim("ARGs_drug_shared_split.txt")
plotdat <- data.frame(var = dat$drug) %>% add_column(perc = map(dat$drug, \(x) 1/length(unlist(strsplit(x, split = ";")))) %>% unlist()) %>% 
  t() %>% data.frame() %>% map_dfr(., \(x) map2_dfr(x[1], as.numeric(x[2]), \(i, j) data.frame(var = unlist(strsplit(i, split = ";")), perc = j))) %>% 
  group_by(var) %>% summarise(n = sum(perc))
write.table(plotdat, "ARGs_drug_shared_split2.txt", sep = "\t", quote = F, row.names = F)
plot_pie(plotdat, dat_colnames = c(name = "var"), top = 12, fill = "auto")
ggsave("ARGs_drug_shared_split_pie.pdf", width = 3, height = 3)

#### ARGs mechanism statistics ####
dat <- read.delim("ARGs_mechanism_manual_classify.txt")
plot_pie(dat, dat_colnames = c(name = "mechanism"), top_n = 12, fill = "auto")
ggsave("ARGs_mechanism_manual_classify_pie.pdf", width = 3, height = 3)

dat <- read.delim("ARGs_mechanism_shared_split.txt")
plotdat <- data.frame(var = dat$mechanism) %>% add_column(perc = map(dat$mechanism, \(x) 1/length(unlist(strsplit(x, split = ";")))) %>% unlist()) %>% 
  t() %>% data.frame() %>% map_dfr(., \(x) map2_dfr(x[1], as.numeric(x[2]), \(i, j) data.frame(var = unlist(strsplit(i, split = ";")), perc = j))) %>% 
  group_by(var) %>% summarise(n = sum(perc))
write.table(plotdat, "ARGs_mechanism_shared_split2.txt", sep = "\t", quote = F, row.names = F)
plot_pie(plotdat, dat_colnames = c(name = "var"), top = 12, fill = "auto")
ggsave("ARGs_mechanism_shared_split_pie.pdf", width = 3, height = 3)

#### prevalence and abundance of ARGs ####
# 计算流行率，tpm<10的不纳入计算，计算丰度就不进行过滤了
source("F:/Code/R_func/profile_process.R")
metadata <- read.delim("ARGs_metadata.txt")
profile <- read.delim("ARGs_tpm", row.names = 1)
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

# 流行率和丰度
prevalence <- profile_replace(profile, 10) %>% profile_prevalence()
abundance <- data.frame(avg_tpm = rowMeans(profile)) %>% rownames_to_column(var = "name")
dat <- merge(prevalence, abundance, by = "name") %>% 
  merge(., y = select(metadata, ARO_accession, ARO_short_name, drug = drug3, mechanism = mechanism3) %>% unique.data.frame(), 
        by.x = "name", by.y = "ARO_accession")
write.table(dat, "ARGs_prevalence.txt", sep = "\t", quote = F, row.names = F)

plotdat <- mutate(dat, drug = forcats::fct_lump_n(drug, 9, ties.method = "last"), avg_tpm = log10(avg_tpm+1e-3))
drug_color <- c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f","#fb8072","#80b1d3","#bc80bd","#e5c494","#b3b3b3")
ggplot(plotdat, aes(x = prevalence, y = avg_tpm)) +
  geom_vline(xintercept = 80, lty = 2, lwd = .5) +
  geom_hline(yintercept = log10(10), lty = 2, lwd = .5) +
  geom_point(aes(color = drug), size = 1.2) +
  scale_color_manual(values = drug_color) +
  labs(x = "Prevalence (%)", y = "Average TPM") +
  geom_text(data = filter(plotdat, prevalence > 80 | avg_tpm > log10(10)), aes(label = ARO_short_name), size = 1) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(size = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)
ggsave("ARGs_prevalence_scatterplot.pdf", height = 4, width = 5)

#### ARGs presence/absence heatmap ####
# overall
library(ComplexHeatmap)
metadata <- read.delim("ARGs_metadata.txt")
profile <- read.delim("ARGs_tpm", row.names = 1)
dat <- profile_replace(profile, 10) %>% profile_adjacency()
rownames(dat) <- metadata$ARO_short_name[match(rownames(dat), metadata$ARO_accession)]
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

annotation_row <- data.frame(name = rownames(dat)) %>% 
  mutate(drug = metadata$drug3[match(name, metadata$ARO_short_name)],
         mechanism = metadata$mechanism3[match(name, metadata$ARO_short_name)],
         drug = forcats::fct_lump_n(drug, 9, ties.method = "last")) %>% column_to_rownames(var = "name")
drug_color <- structure(c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f","#fb8072","#80b1d3","#bc80bd","#e5c494","#b3b3b3"), names = sort(unique(as.character(annotation_row$drug))))
mechanism_color <- structure(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f","#bf5b17"), 
                             names = unique(annotation_row$mechanism))
annotation_col <- data.frame(name = colnames(dat)) %>% 
  mutate(group = group$group[match(name, group$sample)]) %>% 
  column_to_rownames(var = "name")
group_color <- structure(group_color, names = group_order)
annotation_colors <- list(drug = drug_color, mechanism = mechanism_color, group = group_color)
col_split <- factor(annotation_col$group, group_order)

pdf("ARGs_prevalence_overall_heatmap.pdf", width = 10, height = 16)
pheatmap(dat, color = circlize::colorRamp2(c(0, 1), c("white","#74add1")),
         cluster_rows = T, cluster_cols = T, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, column_gap = unit(1,'mm'), border  = "black", border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 3, fontsize_col = 0,
         cellwidth = 4, cellheight = 4, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

# partial
metadata <- read.delim("ARGs_metadata.txt")
profile <- read.delim("ARGs_tpm", row.names = 1)
dat <- profile_replace(profile, 10) %>% profile_adjacency()
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"),names = group_order)

prevalence <- read.delim("ARGs_prevalence.txt") %>% filter(prevalence > 10)
dat <- data.frame(dat[as.character(prevalence$name),])
rownames(dat) <- metadata$ARO_short_name[match(rownames(dat), metadata$ARO_accession)]

annotation_row <- data.frame(name = rownames(dat)) %>% 
  mutate(drug = metadata$drug3[match(name, metadata$ARO_short_name)],
         mechanism = metadata$mechanism3[match(name, metadata$ARO_short_name)],
         drug = forcats::fct_lump_n(drug, 9, ties.method = "last")) %>% 
  column_to_rownames(var = "name")
drug_color <- structure(c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f","#fb8072","#80b1d3","#bc80bd","#cdcc5d","#a2a2a2"), names = sort(unique(as.character(annotation_row$drug))))
mechanism_color <- structure(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"), names = unique(annotation_row$mechanism))
annotation_col <- data.frame(name = colnames(dat)) %>% 
  mutate(group = group$group[match(name, group$sample)]) %>% 
  column_to_rownames(var = "name")
group_color <- structure(group_color, names = group_order)
annotation_colors <- list(drug = drug_color, mechanism = mechanism_color, group = group_color)
col_split <- factor(annotation_col$group, group_order)

pdf("ARGs_prevalence_partial_heatmap.pdf", width = 11, height = 7)
pheatmap(dat, color = circlize::colorRamp2(c(0, 1), c("white","#74add1")),
         cluster_rows = T, cluster_cols = T, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, column_gap = unit(1,'mm'), border  = "black", border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 5, fontsize_col = 0,
         cellwidth = 5, cellheight = 5, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

#### upset plot ####
library(UpSetR)
source("F:/Code/R_func/profile_process.R")
profile <- read.delim("ARGs_tpm", row.names = 1) %>% profile_replace(10)
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

plotdat <- profile_smp2grp(profile, group) %>% profile_adjacency()
upset_metadata <- data.frame(id = group_order, set = group_order)

pdf("ARGs_shared_upset.pdf", width = 6, height = 4)
upset(plotdat, 
      nsets = 9, # 传入的组的数目
      # intersections = group_order, # 指定交集
      # 右上，条形图参数
      number.angles = 0, # 上方柱子文字的倾斜程度
      mainbar.y.label = "ARGs number",  # 上方柱状图的注释
      # 左下，条形图的参数
      sets.x.label = "ARGs number",  # 条形图x轴标签
      sets.bar.color = group_color,  # 左下条形图颜色
      sets = group_order, # 指定特殊的集合, 有排序的作用，但实现排序，需要使用keep.order参数
      keep.order = T, # 排序，按照sets参数设置的顺序
      # 右下，交集散点图的参数
      nintersects = 50, # 右下 交集的数量
      point.size = 1.2,  # 右下点大小
      line.size = 0.3, # 右下线的粗细
      order.by = c("freq"), # 对交集进行排序，默认为升序，freq指定的是按照交集数量的多少排序，degree指的是按照交集集合数量的进行排序。这两个参数是先后顺序
      decreasing = c(T, F), # 调整交集排序的方向，和order.by一一对应，决定是否排序。
      set.metadata = list( # 交集的背景颜色
        data = upset_metadata, 
        plots = list(list(type = "matrix_rows", column = "set", alpha = .4, colors = group_color)) ),
      # 总体参数
      text.scale = c(1.2, 1, 1.1, 1, 1.3, 0.75), # text.scale 参数值的顺序为:柱状图的轴标签和刻度,条形图的轴标签和刻度,集合名称,柱子上方表示交集大小的数值
      mb.ratio = c(.5, .5),  # 控制上下部分图形所占的比例
      )
dev.off()
rownames_to_column(plotdat, var = "name") %>% write.table("ARGs_shared_upset_dat.txt", sep = "\t", quote = F, row.names = F)

#### taxa contribute number of ARGs stacked bar plot ####
source("F:/Code/R_func/profile_process.R")
source("F:/Code/R_func/taxa.R")
metadata <- read.delim("ARGs_metadata.txt")
taxonomy <- read.delim("geneset_CARD_res_taxonomy", header = F, col.names = c("gene", "feature", "taxonomy")) %>% 
  taxa_split(., sep = ";", na_fill = "Unknown", taxa_level = "k:s")

# genus
dat <- select(taxonomy, ARO_accession = feature, genus) %>% unique() %>% 
  merge(select(metadata, ARO_accession, drug = drug3) %>% unique, by = "ARO_accession", all.x = T) %>% 
  mutate(genus = forcats::fct_lump_n(genus, 20, ties.method = "last"), drug = forcats::fct_lump_n(drug, 9, ties.method = "last")) %>% 
  group_by(genus, drug) %>% summarise(n = n()) %>% ungroup() %>% filter(!genus %in% "Other")

genus_order <- group_by(dat, genus) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% select(genus) %>% unlist() %>% as.character()
drug_order <- group_by(dat, drug) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% select(drug) %>% unlist %>% as.character()
plotdat <- mutate(dat, genus = factor(genus, genus_order), drug = factor(drug, drug_order))

drug_color <- c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f","#fb8072","#80b1d3","#bc80bd","#e5c494","#b3b3b3")
ggplot(plotdat, aes(genus, n, fill = drug)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .2, width = .8) +
  scale_y_continuous(expand = c(.02, .02)) +
  scale_fill_manual(values = drug_color) +
  labs(x = "", y = "Number of ARGs", fill = "Drug class") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(linewidth = .4),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 1/2)
ggsave("ARGs_number_with_g_barplot.pdf", width = 8, height = 4)

select(taxonomy, ARO_accession = feature, genus) %>%
  merge(select(metadata, ARO_accession, drug = drug3) %>% unique, by = "ARO_accession", all.x = T) %>% 
  group_by(genus, drug) %>% summarise(n = n()) %>% ungroup() %>%
  write.table("ARGs_number_with_g_dat.txt", sep = "\t", row.names = F, quote = F)
write.table(dat, "ARGs_number_with_g_dat2.txt", sep = "\t", row.names = F, quote = F)

# species
dat <- select(taxonomy, ARO_accession = feature,species) %>% unique() %>% 
  merge(select(metadata, ARO_accession, drug = drug3) %>% unique, by = "ARO_accession", all.x = T) %>% 
  mutate(species = forcats::fct_lump_n(species, 20, ties.method = "last"), drug = forcats::fct_lump_n(drug, 9, ties.method = "last")) %>% 
  group_by(species, drug) %>% summarise(n = n()) %>% ungroup() %>% filter(!species %in% "Other")

species_order <- group_by(dat, species) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% select(species) %>% unlist() %>% as.character()
drug_order <- group_by(dat, drug) %>% summarise(n = sum(n)) %>% arrange(desc(n)) %>% select(drug) %>% unlist %>% as.character()
plotdat <- mutate(dat, species = factor(species, species_order), drug = factor(drug, drug_order))

drug_color <- c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f","#fb8072","#80b1d3","#bc80bd","#e5c494","#b3b3b3")
ggplot(plotdat, aes(species, n, fill = drug)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .2, width = .8) +
  scale_y_continuous(expand = c(.02, .02)) +
  scale_fill_manual(values = drug_color) +
  labs(x = "", y = "Number of ARGs", fill = "Drug class") +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid.major = element_line(linewidth = .4),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8, color = "black"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 1/2)
ggsave("ARGs_number_with_s_barplot.pdf", width = 8, height = 4)

select(taxonomy, ARO_accession = feature, species) %>%
  merge(select(metadata, ARO_accession, drug = drug3) %>% unique, by = "ARO_accession", all.x = T) %>% 
  group_by(species, drug) %>% summarise(n = n()) %>% ungroup() %>%
  write.table("ARGs_number_with_s_dat.txt", sep = "\t", row.names = F, quote = F)
write.table(dat, "ARGs_number_with_s_dat2.txt", sep = "\t", row.names = F, quote = F)


#### sankey ####
library(ggsankey)
source("F:/Code/R_func/profile_process.R")
source("F:/Code/R_func/taxa.R")
metadata <- read.delim("ARGs_metadata.txt")
taxonomy <- read.delim("geneset_CARD_res_taxonomy", header = F, col.names = c("gene", "feature", "taxonomy")) %>% 
  taxa_split(., sep = ";", na_fill = "Unknown", taxa_level = "k:s")

dat <- select(taxonomy, ARO_accession = feature, phylum, family, genus, species) %>% unique %>%  
  merge(select(metadata, ARO_accession, drug = drug3, mechanism = mechanism3) %>% unique, by = "ARO_accession", all.x = T) %>% 
  mutate(family = forcats::fct_lump_n(family, 19, other_level = "f__Other", ties.method = "last"), 
         genus = forcats::fct_lump_n(genus, 19, other_level = "g__Other", ties.method = "last"), 
         species = forcats::fct_lump_n(species, 19, other_level = "s__Other", ties.method = "last"), 
         drug = forcats::fct_lump_n(drug, 19, other_level = "Other", ties.method = "last"))

phylum_order <- group_by(dat, phylum) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(phylum) %>% unlist() %>% as.character()
family_order <- group_by(dat, family) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(family) %>% unlist() %>% as.character()
genus_order <- group_by(dat, genus) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(genus) %>% unlist() %>% as.character()
species_order <- group_by(dat, species) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(species) %>% unlist() %>% as.character()
drug_order <- group_by(dat, drug) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(drug) %>% unlist %>% as.character()
mechanism_order <- group_by(dat, mechanism) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(mechanism) %>% unlist %>% as.character()

plotdat <- make_long(select(dat, -1), phylum, family, genus, species, drug, mechanism) %>% 
  mutate(node = as.factor(node), 
         node = forcats::fct_relevel(node, rev(c(phylum_order, family_order, genus_order, species_order, drug_order, mechanism_order))))

ggplot(plotdat, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey(node.color = "black", width = .1, flow.fill = "grey85", flow.alpha = .5, lwd = .2,
              space = 50, smooth = 8, show.legend = F) +
  geom_sankey_text(size = 2, color = "black", width = .1, space = 50, hjust = 0) +
  scale_x_discrete(label = c("Phylum", "Family", "Genus", "Species", "Drug", "Mechanism")) +
  labs(x = "") +
  theme_sankey(base_size = 10) +
  theme(axis.text.x = element_text(color = "black", size = 8))
ggsave("ARGs_number_with_taxa_sankey.pdf", width = 12, height = 8)

dat <- select(taxonomy, ARO_accession = feature, phylum, family, genus, species) %>%
  merge(select(metadata, ARO_accession, drug = drug3, mechanism = mechanism3) %>% unique, by = "ARO_accession", all.x = T)
select(dat, -1) %>% group_by(phylum, family, genus, species, drug, mechanism) %>% 
  summarise(n = n()) %>% arrange(desc(n)) %>% 
  write.table("ARGs_number_with_taxa_sankey_dat2.txt", sep = "\t", quote = F, row.names = F)

# species
dat <- select(taxonomy, ARO_accession = feature, phylum, species) %>%
  merge(select(metadata, ARO_accession, drug = drug3, mechanism = mechanism3) %>% unique, by = "ARO_accession", all.x = T) %>% 
  mutate(species = forcats::fct_lump_n(species, 19, other_level = "s__Other"), drug = forcats::fct_lump_n(drug, 10))
phylum_order <- group_by(dat, phylum) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(phylum) %>% unlist() %>% as.character()
species_order <- group_by(dat, species) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(species) %>% unlist() %>% as.character()
drug_order <- group_by(dat, drug) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(drug) %>% unlist %>% as.character()
mechanism_order <- group_by(dat, mechanism) %>% summarise(n = n()) %>% arrange(desc(n)) %>% select(mechanism) %>% unlist %>% as.character()
plotdat <- make_long(select(dat, -1), phylum, species, drug, mechanism) %>% 
  mutate(node = as.factor(node), 
         node = forcats::fct_relevel(node, rev(c(phylum_order, species_order, drug_order, mechanism_order))))
ggplot(plotdat, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = node, label = node)) +
  geom_sankey(node.color = "black", width = .1, flow.fill = "grey70", flow.alpha = .5, lwd = .2,
              space = 50, smooth = 8, show.legend = F) +
  geom_sankey_text(size = 2, color = "black", width = .1, space = 50, hjust = 0) +
  scale_x_discrete(label = c("Phylum", "Species", "Drug", "Mechanism")) +
  labs(x = "") +
  theme_sankey(base_size = 10) +
  theme(axis.text.x = element_text(color = "black", size = 8))
ggsave("ARGs_number_with_s_sankey.pdf", width = 8, height = 6)
select(dat, -1) %>% group_by(phylum, species, drug, mechanism) %>% 
  summarise(n = n()) %>% arrange(desc(n)) %>% 
  write.table("ARGs_number_with_s_sankey_dat.txt", sep = "\t", quote = F, row.names = F)

#### overall abundance ####
source("F:/Code/R_func/diversity.R")
profile <- read.delim("ARGs_tpm", row.names = 1) %>% apply(., 2, \(x) x/1e6*100) %>% data.frame()
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

plotdat <- data.frame(val = colSums(profile)) %>% rownames_to_column(var = "sample") %>% 
  merge(group, by = "sample") %>% mutate(group = factor(group, group_order)) 
write.table(plotdat, "div_overall_abundance_dat.txt", sep = "\t", quote = F, row.names = F)
comparison <- diff_test(select(plotdat, -group), group, group_order, method = "wilcox") %>% 
  filter(pval < 0.05) %>% select(group_pair) %>% unlist() %>% 
  as.character() %>% strsplit(x = ., split = "_vs_")

ggplot(plotdat, aes(group, val)) +
  stat_compare_means(comparisons = comparison, method = "wilcox", method.args = list(exact = F),
                     label = "p.signif", tip.length = .01, step.increase = .03, vjust = .95, size = 3) +
  geom_boxplot(width = .5, fill = "transparent", lwd = .4, outlier.shape = NA, show.legend = F) +
  geom_jitter(data = plotdat, aes(x = group, y = val, color = group), width = .3,
              size = 1, inherit.aes = F, show.legend = F) +
  stat_smooth(data = plotdat, aes(group = 1), level = 0.95, method = "lm",
              color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0, label.y.npc = 1, color = "#8da0cb") +
  stat_smooth(data = plotdat %>% filter(group %in% c("d1","d3")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = plotdat %>% filter(group %in% c("d1","d3")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0.5, label.y.npc = 1, color = "#fc8d62") +
  stat_smooth(data = plotdat %>% filter(!group %in% c("d1","d3")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = plotdat %>% filter(!group %in% c("d1","d3")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 1, label.y.npc = 1, color = "#66c2a5") +
  scale_color_manual(values = group_color) +
  labs(x = "", title = "Overall abundance (%)", y = "") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "#000000"),
        axis.ticks = element_line(linewidth = .4, color = "#000000"),
        axis.text = element_text(size = 8, color = "#000000"),
        axis.title = element_text(size = 8, color = "#000000"),
        plot.title = element_text(size = 10, color = "#000000"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 4/3)
ggsave("div_overall_abundance_boxplot2.pdf", height = 4, width = 3)

#### α-diversity ####
source("F:/Code/R_func/diversity.R")
source("F:/Code/R_func/profile_process.R")
profile <- read.delim("ARGs_tpm", row.names = 1) %>% profile_replace(10)
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

method = "richness"
plotdat <- calu_alpha(profile, method = method) %>% 
  merge(., group, by = "sample", all.x = T) %>% 
  mutate(group = factor(group, group_order))

comparison <- diff_test(select(plotdat, -group), group, method = "wilcox", group_order = group_order) %>% 
  filter(pval < 0.05) %>% select(group_pair) %>% unlist() %>% as.character() %>% 
  strsplit(x = ., split = "_vs_")

ggplot(plotdat, aes(group, val)) +
  stat_compare_means(comparisons = comparison, method = "wilcox", method.args = list(exact = F), 
                     label = "p.signif", tip.length = .01, step.increase = .03, vjust = .95, size = 3) +
  geom_boxplot(width = .5, fill = "transparent", lwd = .4, outlier.shape = NA, show.legend = F) +
  geom_jitter(aes(group, val, color = group), width = .3, size = .8, inherit.aes = F, show.legend = F) +
  stat_smooth(aes(group = 1), method = "lm", color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1.5, label.x.npc = 1, label.y.npc = 0, color = "#8da0cb") +
  stat_smooth(data = filter(plotdat, group %in% c("d1","d3")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1.5, label.x.npc = 1, label.y.npc = .03, color = "#fc8d62") +
  stat_smooth(data = filter(plotdat, !group %in% c("d1","d3")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1.5, label.x.npc = 1, label.y.npc = .06, color = "#66c2a5") +
  scale_color_manual(values = group_color) +
  labs(x = "", title = paste0(stringr::str_to_title(method), " Index"), y = "") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 10, color = "black"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 4/3)
ggsave(paste0("div_", method, "_boxplot_latest.pdf"), height = 4, width = 3)

plotdat %>% group_by(group) %>% summarise(mean = mean(val), sd = sd(val), count = n()) %>% 
  write.table(., paste0("div_", method, "_latest.txt"), sep = "\t", row.names = F, quote = F)
write("", paste0("div_", method, "_latest.txt"), append = T)
diff_test(select(plotdat, -group), group, method =  "wilcox", group_order = group_order) %>% 
  mutate(plab = ifelse(pval < 0.0001, "****", ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "!"))))) %>% 
  write.table(., paste0("div_", method, "_latest.txt"), sep = "\t", append = T, row.names = F, quote = F)
write("", paste0("div_", method, "_latest.txt"), append = T)
write.table(plotdat, paste0("div_", method, "_latest.txt"), sep = "\t", append = T, row.names = F, quote = F)

#### PCoA ####
source("F:/Code/R_func/plot_PCoA.R")
p <- plot_PCoA(profile, group, group_order = group_order, group_color = group_color)
ggsave("div_PCoA_scatterplot.pdf", height = 4.5, width = 6)
write.table(p$data, "div_PCoA_dat.txt", sep = "\t", quote = F, row.names = F)

#### ARGs abundance heatmap change ####
source("F:/Code/R_func/profile_process.R")
metadata <- read.delim("ARGs_metadata.txt")
profile <- read.delim("ARGs_tpm", row.names = 1)
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

dat <- read.delim("ARGs_prevalence.txt") %>% filter(prevalence > 20)
plotdat <- data.frame(profile[as.character(dat$name),]) %>% profile_smp2grp(., group) %>% relocate(all_of(group_order))
rownames(plotdat) <- dat$ARO_short_name[match(rownames(plotdat), dat$name)]

# diff
diff <- diff_test_profile(data.frame(profile[as.character(dat$name),]), group, group_order, add_plab = T)
diff <- filter(diff, group_pair %in% c("d1_vs_d3", "d3_vs_d5", "d5_vs_d7", "d7_vs_d14", "d14_vs_d21", "d21_vs_d28", "d28_vs_d35", "d35_vs_d42")) %>% 
  mutate(plab = ifelse(plab != "", "", "*"), 
         gp1 = lapply(strsplit(group_pair, "_vs_"), "[[", 1),
         gp2 = lapply(strsplit(group_pair, "_vs_"), "[[", 2)) %>% select(-pval, -method, -group_pair, -gp1) %>% 
  spread(key = "gp2", value = "plab") %>% column_to_rownames(var = "feature") %>% add_column(d1 = "") %>% relocate(d1)
rownames(diff) <- dat$ARO_short_name[match(rownames(diff), dat$name)]

# annotation
annotation_row <- data.frame(name = rownames(plotdat)) %>% 
  mutate(drug = metadata$drug3[match(name, metadata$ARO_short_name)],
         mechanism = metadata$mechanism3[match(name, metadata$ARO_short_name)]) %>% 
  column_to_rownames(var = "name")
drug_color <- structure(c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f","#fb8072","#80b1d3","#bc80bd","#e5c494","#e5c449","#b3b3b3"),
                        names = sort(unique(annotation_row$drug)))
mechanism_color <- structure(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"), names = unique(annotation_row$mechanism))
annotation_colors <- list(drug = drug_color, mechanism = mechanism_color)

plotdat2 <- mutate(plotdat, d3 = d3 - d1, d5 = d5 - d3, d7 = d7 - d5, d14 = d14 - d7, d21 = d21 - d14, 
                   d28 = d28 - d21, d35 = d35 - d28, d42 = d42 - d35, d1 = 0) %>%
  apply(., 2, \(x) ifelse(x == 0, "·", ifelse(x > 0, "+", "-"))) %>% data.frame() %>% 
  map2_df(., diff, \(x, y) paste0(x, "", y)) %>% as.matrix()
rownames(plotdat2) <- rownames(diff)

pdf("ARGs_abundance_partial_heatmap.pdf", width = 6, height = 6)
pheatmap(log10(plotdat+1e-3), color = circlize::colorRamp2(c(0, 3), c("white","#d53e4f")),
         cluster_rows = T, cluster_cols = F, treeheight_col = 20, treeheight_row = 20,
         border = "black", border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = T, fontsize_row = 5, fontsize_col = 5,
         cellwidth = 15, cellheight = 8, display_numbers = plotdat2, number_color = "black",
         annotation_row = annotation_row, annotation_colors = annotation_colors,
         heatmap_legend_param = list(title = "avg.TPM (Log10)", title_gp = grid::gpar(fontsize = 10), border = "black"))
dev.off()

#### resistome stacked bar flow ####
source("F:/Code/R_func/profile_process.R")
metadata <- read.delim("ARGs_metadata.txt")
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
profile <- read.delim("ARGs_tpm", row.names = 1) %>% profile_smp2grp(., group) %>% 
  mutate(drug = metadata$drug3[match(rownames(.), metadata$ARO_accession)]) %>% 
  group_by(drug) %>% summarise_all(sum) %>% ungroup() %>% column_to_rownames(var = "drug")

# dat <- map2_dfr(profile$drug, data.frame(t(profile[,-1])), \(x, y)
#            cbind(data.frame(drug = unlist(strsplit(x, split = ";"))), 
#                  data.frame(t(y/(stringr::str_count(x, ";") + 1)))[rep(1, (stringr::str_count(x, ";") + 1)),])) %>% 
#   group_by(drug) %>% summarise_all(sum) %>% ungroup() %>% column_to_rownames(var = "drug")
# colnames(dat) <- colnames(profile)[-1]

tmp_vec <- sort(rowSums(profile), decreasing = T) %>% head(9) %>% names()
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
plotdat <- relocate(profile, all_of(group_order)) %>% rownames_to_column(var = "drug") %>% 
  mutate(drug = ifelse(drug %in% tmp_vec, drug, "Other")) %>% 
  group_by(drug) %>% summarise_all(sum) %>% column_to_rownames(var = "drug") %>% 
  apply(., 2, \(x) x/sum(x)*100) %>% data.frame() %>% rownames_to_column(var = "drug") %>% 
  gather(key = "group", value = "val", -drug) %>% mutate(group = factor(group, group_order))
drug_color <- structure(c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f","#fb8072","#80b1d3","#bc80bd","#e5c449","#b3b3b3"), names = sort(unique(plotdat$drug)))

library(ggalluvial)
ggplot(plotdat, aes(x = group, y = val, alluvium = drug, stratum = drug)) +
  geom_stratum(aes(fill = drug), width = .4, lwd = .5, color = "transparent") +
  geom_alluvium(aes(fill = drug), alpha = 1, width = .2) +
  scale_fill_manual(values = drug_color) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "", y = "Relative abundance (%)") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = .4, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)
ggsave("ARGs_abundance_drug_dynamic_change.pdf", width = 6, height = 4.5)

profile <- read.delim("ARGs_tpm_m10", row.names = 1) %>%
  mutate(drug = metadata$drug3[match(rownames(.), metadata$ARO_accession)]) %>% 
  group_by(drug) %>% summarise_all(sum) %>% ungroup() %>% column_to_rownames(var = "drug")
source("F:/Code/R_func/diff_test.R")
diff <- diff_test_profile(profile, group, group_order, method = "wilcox", add_plab = T)
write.table(diff, "ARGs_abundance_drug_dynamic_change_diff_dat.txt", sep = "\t", row.names = F, quote = F)

#### resistome fit ####
source("F:/Code/R_func/profile_process.R")
metadata <- read.delim("ARGs_metadata.txt")
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
profile <- read.delim("ARGs_tpm", row.names = 1) %>% 
  mutate(drug = metadata$drug3[match(rownames(.), metadata$ARO_accession)]) %>% 
  group_by(drug) %>% summarise_all(sum) %>% ungroup() %>% column_to_rownames(var = "drug")

plotdat <- profile_filter(profile, group, by_group = T, min_prevalence = 0.5) %>% 
  rownames_to_column(var = "drug") %>% gather(key = "sample", value = "val", -drug) %>% 
  filter(val != 0) %>% 
  mutate(group = group$group[match(sample, group$sample)],
         group = factor(group, group_order), val = log10(val)) %>% 
  filter(!drug %in% c("MLS","aminoglycoside","carbapenem","cephalosporin","elfamycin", "diaminopyrimidine",
                     "glycopeptide","macrolide", "macrolide","multi-drug","phenicol","pleuromutilin","tetracycline"))
  # filter(drug != "oxazolidinone") %>% 
  # filter(!group %in% c("d1", "d3", "d5"))

# drug_color <- c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f","#fb8072","#80b1d3","#bc80bd","#e5c449","#7fc97f",
#                 "#fdc086","#f0027f","#bf5b17","#ff7f00","#984ea3","#4daf4a","#377eb8","#e41a1c","#b3b3b3")

# lm
ggplot(plotdat, aes(x = group, y = val, color = drug)) +
  geom_jitter(size = 1, width = .4) +
  geom_smooth(aes(group = drug), method = "lm", se = F, linewidth = 1) +
  stat_poly_eq(aes(group = drug, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')), size = 1) +
  scale_color_manual(values = drug_color) +
  labs(x = "", y = "TPM (log10)") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = .4, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)

# lm
ggplot(plotdat, aes(x = group, y = val, color = drug)) +
  facet_wrap(~factor(drug), nrow = 3, scales = "free") +
  geom_jitter(size = 1, width = .4, show.legend = F) +
  stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5")), aes(group = 1), 
              color = "grey65", linewidth = .5, method = "lm", se = F, lty = "dashed", show.legend = F) +
  stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5")),
               aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = .5, label.y.npc = .04) +
  stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5")), aes(group = 1), 
              color = "grey65", linewidth = .5, method = "lm", se = F, lty = "dashed", show.legend = F) +
  stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5")),
               aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = .5, label.y.npc = .07) +
  geom_smooth(aes(group = drug), method = "lm", se = F, linewidth = .5, color = "grey65", show.legend = F) +
  stat_poly_eq(aes(group = drug, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')), 
               size = 1, label.x.npc = .5, label.y.npc = .01) +
  labs(x = "", y = "TPM (log10)") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = .4, color = "black"),
        panel.grid = element_blank(),
        strip.text = element_text(size = 8, color = "black", face = "italic", vjust = .5),
        strip.background = element_rect(linewidth = .4, color = "black"))
ggsave("ARGs_abundance_change_lm_satterplot_latest_other.pdf", width = 8, height = 6)

# loess
ggplot(plotdat, aes(x = group, y = val, color = drug)) +
  geom_jitter(size = 1, width = .4) +
  geom_smooth(aes(group = drug), method = "loess", se = F, linewidth = 1) +
  scale_color_manual(values = drug_color) +
  labs(x = "", y = "TPM (log10)") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.border = element_rect(linewidth = .4, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)
ggsave("ARGs_abundance_change_loess_satterplot.pdf", width = 8, height = 6)

#### corr ####
source("F:/Code/R_func/profile_process.R")
source("F:/Code/R_func/corr_process.R")
source("F:/Code/R_func/taxa.R")
metadata <- read.delim("ARGs_metadata.txt")
taxonomy <- read.delim("../taxa/species_taxonomy_rename.tsv")
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

prevalence <- read.delim("ARGs_prevalence.txt") %>% filter(prevalence > 20)
  
profile_args <- read.delim("../ARGs/ARGs_tpm", row.names = 1)[as.character(prevalence$name),] %>% 
  relocate(group$sample) %>% t() %>% data.frame(check.names = F)
profile_bac <- read.delim("../taxa/species_tpm", row.names = 1) %>% 
  taxa_trans(taxonomy, to = "species", out_all = T) %>% 
  profile_filter(min_prevalence = 0.05, by_group = F) %>% 
  relocate(group$sample) %>% t() %>% data.frame(check.names = F)

# 相关性分析
corr <- psych::corr.test(x = profile_args, y = profile_bac, method = "spearman", adjust = "BH")
corr_res <- corr_process(corr, rho = 0.5, padj = 0.05, out_df = T) %>% 
  mutate(direct = ifelse(r > 0, "pos", "neg"))
corr_res$name_x <- metadata$ARO_short_name[match(corr_res$name_x, metadata$ARO_accession)]



