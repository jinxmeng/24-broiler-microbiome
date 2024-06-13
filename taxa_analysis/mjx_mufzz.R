#### info ####
# Jinxin Meng, 20231009, 20231228
pacman::p_load(tibble, dplyr, tidyr, ggplot2)
setwd("F:/Microbiome/proj_2023/02.Broiler_fecal_resistome_20230302/taxa_explore_to_fluctations/")

#### survey species with family ####
library(ComplexHeatmap)
source("F:/Code/R_func/taxa.R")
group <- read.delim("../data/sample_group.txt")
taxonomy <- read.delim("../taxa/species_taxonomy_rename.tsv") 
profile <- read.delim("../taxa/species_tpm_filter", row.names = 1)%>% taxa_trans(taxonomy, to = "species", out_all = T)
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
taxonomy <- read.delim("../taxa/species_taxonomy_rename.tsv") 
profile <- read.delim("../taxa/species_tpm_filter", row.names = 1)%>% taxa_trans(taxonomy, to = "species", out_all = T)
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
taxonomy <- read.delim("../taxa/species_taxonomy_rename.tsv")
profile <- read.delim("../taxa/species_tpm_filter", row.names = 1) %>% taxa_trans(taxonomy, to = "genus", out_all = T)
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
taxonomy <- read.delim("../taxa/species_taxonomy_rename.tsv")
profile <- read.delim("../taxa/species_tpm_filter", row.names = 1) %>% taxa_trans(taxonomy, to = "family", out_all = T)
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
taxonomy <- read.delim("../taxa/species_taxonomy_rename.tsv") %>% 
  mutate(species =  paste0(species," ",feature))
profile <- read.delim("../taxa/species_tpm_filter", row.names = 1) %>% taxa_trans(taxonomy, to = "species", out_all = T, transRA = T)
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

# 保存数据
# merge(cluster_dat, taxonomy, by.x = "taxa", by.y = "species") %>% relocate(cluster) %>% 
#   mutate(species2 = strsplit(taxa, " ") %>% purrr::map_vec(\(x) x[-length(x)] %>% paste(collapse = " "))) %>% 
#   relocate(species2, .after = "taxa") %>% 
#   write.table("tmp_20231226_mufzz_cluster_dat.tsv", sep = "\t", row.names = F, quote = F)


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
          axis.ticks = element_line(color = "black", linewidth = .5),
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
write.table(plotdat, "tmp_20231226_mufzz_taxa_f__Ruminococcaceae.tsv", sep = "\t", row.names = F, quote = F)
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