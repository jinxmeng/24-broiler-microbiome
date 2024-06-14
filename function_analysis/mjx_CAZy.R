# Jinxin Meng, 20231010, 20240611 ------------------------------------------------

pacman::p_load(tidyr, dplyr, tibble, purrr, microeco, ggpubr)
pacman::p_load(ComplexHeatmap)
setwd("F:/Proj/proj_2023/01.Broiler_fecal_resistome_20230302/CAZy/")

# lefse -------------------------------------------------------------------------

source("F:/Code/R_func/profile_process.R")
map <- read.delim("../data/sample_group.txt")

# 数据过滤，至少在一组有三个样本有值
dat <- read.delim("CAZy_tpm_filter_m10", row.names = 1) %>% 
  profile_filter(group = map, by_group = T, min_prevalence = 0.35)
  
# lefse
dataset <- microtable$new(sample_table = map %>% column_to_rownames(var = "sample"),
                          otu_table = dat, tax_table = data.frame(row.names = rownames(dat), name = rownames(dat)))
lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = "group", alpha = 0.05, p_adjust_method = "BH")

lefse_res <- lefse$res_diff %>% 
  data.frame(row.names = NULL) %>% 
  select(name = Taxa, enriched = Group, LDA, pval = P.unadj, padj = P.adj, plab = Significance)
write.table(lefse_res, "CAZy_lefse_dat.txt", sep = "\t", row.names = F, quote = F)

# 热图
lefse_res <- read.delim("CAZy_lefse_dat.txt") %>% 
  filter(LDA > 2 & !grepl("\\.", x = name))
dat <- read.delim("CAZy_tpm_filter_m10", row.names = 1) %>% 
  filter(rownames(.) %in% lefse_res$name)

# 颜色
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)
class_order <- c("AA","CBM","CE","GH","GT","PL")
class_color <- structure(c("#6EE2FF","#F7C530","#95CC5E","#D0DFE6","#F79D1E","#748AA6"), names = class_order)

# 注释数据
col_annotation_dat <- data.frame(name = rownames(dat)) %>% 
  mutate(class = gsub("\\d+", "", x = name),
         class_color = class_color[match(class, names(class_color))]) %>% 
  merge(., select(lefse_res, name, LDA, enriched, plab), by = "name", all.x = T) %>% 
  mutate(LDA = ifelse(is.na(LDA), 0, LDA),
         plab = ifelse(is.na(plab), "", plab),
         enriched_color = group_color[match(enriched, names(group_color))])
row_annotation_dat <- data.frame(name = colnames(dat)) %>% 
  mutate(group = map$group[match(name, map$sample)],
         group_color = group_color[match(group, names(group_color))]) 

# 行列注释
annotation_col <- col_annotation_dat %>% 
  select(name, class) %>% 
  column_to_rownames(var = "name")
annotation_row <- row_annotation_dat %>% 
  select(name, group) %>% 
  column_to_rownames(var = "name")
annotation_colors <- list(group = select(row_annotation_dat, group, group_color) %>% 
                            unique() %>% 
                            pull(name = group), 
                          class = select(col_annotation_dat, class, class_color) %>% 
                            unique() %>% 
                            pull(name = class))

# 分块
col_split <- factor(col_annotation_dat$class)
row_split <- factor(row_annotation_dat$group, group_order)

# 柱状图
col_bar <- HeatmapAnnotation(LDA = anno_barplot(col_annotation_dat$LDA, gp = gpar(fill = col_annotation_dat$enriched_color, col = NA)))
col_sig <- HeatmapAnnotation(pval = anno_text(col_annotation_dat$plab, rot = 90, just = c(0.5, 1), gp = gpar(col = "red", fontsize = 5)))

# 小格子颜色
color <- c("#ffffff","#fff5f0","#fee0d2","#fcbba1","#fc9272","#fb6a4a","#ef3b2c")

# 主图
# pdf(file = "CAZy_heatmap_LDA.pdf", height = 8, width = 18)
pheatmap(t(log10(dat+1)),
        scale = "none", heatmap_legend_param = list(title = "TPM (Log10)"),
        cluster_rows = F, cluster_cols = T,
        treeheight_row = 15, treeheight_col = 15,
        show_rownames = F, show_colnames = F,
        cellwidth = 4, cellheight = 4,
        color = color, 
        border_color = "#ffffff",
        row_split = row_split, column_split = col_split,
        column_gap = unit(0.5, "mm"), row_gap = unit(0.5, "mm"),
        annotation_row = annotation_row, annotation_colors = annotation_colors,
        top_annotation = col_bar, bottom_annotation = col_sig)
# dev.off()

group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

lefse_res <- read.delim("CAZy_lefse_dat.txt") %>% 
  filter(LDA > 2 & !grepl("\\.", x = name))

table(lefse_res$enriched) %>%
  data.frame() %>% 
  select(sample = Var1, val = Freq) %>%
  mutate(sample = factor(sample, group_order)) %>%
  ggpubr::ggbarplot(., x = "sample", y = "val", fill = "sample", xlab = "", ylab = "CAZymes count") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = group_color) +
  theme(legend.position = "none")
ggsave("CAZy_lefse_CAZymes_count_barplot.pdf", width = 4, height = 4)

# diff ----------------------------------------------------------------------------

source("F:/Code/R_func/diff_test.R")
map <- read.delim("../data/sample_group.txt")
dat <- read.delim("CAZy_tpm_filter_m10", row.names = 1) %>% 
  profile_filter(group = map, by_group = T, min_prevalence = 0.35)

dat2 <- data.frame(t(dat), check.names = F) %>% 
  mutate(group = map$group[match(rownames(.), map$sample)]) %>% 
  group_by(group) %>% 
  summarise_all(mean) %>% 
  column_to_rownames(var = "group") %>% 
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "name")
write.table(dat2, "CAZy_abundance_dat.txt", sep = "\t", row.names = F, quote = F)

diff <- diff_test_profile(dat, map, group_order)
write.table(diff, "CAZy_diff_dat.txt", sep = "\t", quote = F, row.names = F)

# lefse substrate --------------------------------------------------------------------

fam_sub <- read.delim("fam-substrate-mapping-08252022.tsv")
CAZy_sub <- read.delim("CAZy_lefse_dat.txt") %>% 
  filter(LDA > 2 & !grepl("\\.", x = name)) %>%
  merge(fam_sub, by.x = "name", by.y = "Family")
# write.table(CAZy_sub, "CAZy_lefse_substrate_dat.txt", sep = "\t", quote = F, row.names = F)

CAZy_sub <- group_by(CAZy_sub, name, enriched, LDA, padj, plab) %>% 
  group_modify(~ { 
    .x %>% 
      head(n = 1) %>% 
      mutate(substrate = paste(unique(.x$Substrate_high_level), collapse = "; "), .after = "pval") 
    }) %>% 
  ungroup() %>% 
  select(name, enriched, LDA, padj, plab, substrate)
# write.table(CAZy_sub, "CAZy_lefse_substrate_dat2.txt", quote = F, sep = "\t", row.names = F)

CAZy_sub <- group_by(CAZy_sub, name, enriched, LDA, padj, plab) %>% 
  group_modify(~ { 
    .x %>% 
      head(n = 1) %>% 
      mutate(substrate_high = paste(unique(.x$Substrate_high_level), collapse = "; "), .after = "pval") %>% 
      mutate(substrate_simp = paste(unique(.x$Substrate_Simple), collapse = "; "), .after = "substrate_high")
  }) %>% 
  ungroup() %>% 
  select(name, enriched, LDA, padj, plab, substrate_high, substrate_simp)
# write.table(CAZy_sub, "CAZy_lefse_substrate_dat3.txt", quote = F, sep = "\t", row.names = F)

# 底物类型统计 high level

source("F:/Code/R_func/plot_pie.R")
plotdat <- CAZy_sub %>% 
  select(enriched, substrate) %>% 
  add_column(perc = map_int(.$substrate, \(x) 
                            strsplit(x, split = "[;|,] ") %>% 
                              unlist() %>% 
                              length())) %>% 
  group_by(enriched) %>% 
  group_modify(~ map_df(data.frame(t(.x)), \(x) 
                        map2_dfr(unlist(strsplit(x[1], split = "[;|,] ")), as.numeric(x[2]), \(i, j) 
                                 data.frame(substrate = i, n = 1/j) ) ) %>% 
                 group_by(substrate) %>% 
                 summarise(n = sum(n))
               )
# write.table(plotdat, "CAZy_lefse_substrate_stat.txt", sep = "\t", quote = F, row.names = F)

group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
p_list <- list()
for (i in group_order) {
  tmp_plotdat <- filter(plotdat, enriched %in% i) %>% mutate(n = round(n, 2))
  p_list[[i]] <- plot_pie(tmp_plotdat, dat_colnames = c(name = "substrate"), top_n = 10, add_n = T, add_perc = F, fill = "auto", title = i)
}
cowplot::plot_grid(plotlist = p_list, nrow = 2)
ggsave("CAZy_lefse_substrate_pie2.pdf", width = 15, height = 6)

plotdat %>% 
  mutate(enriched = factor(enriched, group_order)) %>% 
  ggscatter(x = "enriched", y = "substrate", color = "n", size = "n") +
  scale_size(range = c(1, 4)) +
  scale_colour_viridis_c(begin = .3)

# 底物类型统计 simp level

plotdat <- CAZy_sub %>% 
  select(enriched, substrate_simp) %>% 
  add_column(perc = map_int(.$substrate_simp, \(x) 
                            strsplit(x, split = "[;|,] ") %>% 
                              unlist() %>% 
                              length())) %>% 
  group_by(enriched) %>% 
  group_modify(~ map_df(data.frame(t(.x)), \(x) 
                        map2_dfr(unlist(strsplit(x[1], split = "[;|,] ")) %>% trimws, as.numeric(x[2]) , \(i, j) 
                                 data.frame(substrate_simp = i, n = 1/j) ) ) %>% 
                 group_by(substrate_simp) %>% 
                 summarise(n = sum(n))
  ) %>% 
  mutate(substrate_high = fam_sub$Substrate_high_level[match(substrate_simp, fam_sub$Substrate_Simple)])
# # write.table(plotdat, "CAZy_lefse_substrate_stat.txt", sep = "\t", quote = F, row.names = F)
# 
# group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
# p_list <- list()
# for (i in group_order) {
#   tmp_plotdat <- filter(plotdat, enriched %in% i) %>% mutate(n = round(n, 2))
#   p_list[[i]] <- plot_pie(tmp_plotdat, dat_colnames = c(name = "substrate"), top_n = 10, add_n = T, add_perc = F, fill = "auto", title = i)
# }
# cowplot::plot_grid(plotlist = p_list, nrow = 2)
# ggsave("CAZy_lefse_substrate_pie2.pdf", width = 15, height = 6)

plotdat %>% 
  mutate(enriched = factor(enriched, group_order)) %>% 
  ggscatter(x = "enriched", y = "substrate_simp", color = "n", size = "n") +
  scale_size(range = c(1, 4)) +
  scale_colour_viridis_c(begin = .3)

#### tmp lefse substrate ####
# pul <- read.delim("dbCAN_PUL.txt")
# pul_map <- map2_dfr(strsplit(pul$cazymes, split = ","), pul$id, \(x, y)
#                     strsplit(x, split = "\\|") %>% expand.grid()
#                     %>% t() %>% data.frame() %>%
#                       map(., \(i) paste(unique(i), collapse = " ")) %>%
#                       data.frame() %>% t() %>%
#                       data.frame(CAZy = ., row.names = NULL) %>%
#                       add_column(id = y, .before = 1))
# write.table(pul_map, "dbCAN_PUL_map.tsv", sep = "\t", quote = F, row.names = F)
# 
# lefse_res <- read.delim("CAZy_lefse_dat.txt") %>% filter(LDA > 2)
# pul_map <- read.delim("dbCAN_PUL_map.tsv", sep = "\t") %>% 
#   filter(!grepl("_", x = CAZy))
# 
# group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
# pul_res <- rbind()
# for (i in group_order) {
#   x = filter(lefse_res, enriched == i) %>% select(name) %>% unlist() %>% as.character()
#   res <- data.frame(id = character(), CAZy = character())
#   for (j in 1:nrow(pul_map)) {
#     if (all(unlist(stringr::str_split(pul_map[j,2], " ")) %in% x)) {
#       res <- add_row(res, id = pul_map[j, 1], CAZy = pul_map[j, 2])
#     }
#   }
#   res <- add_column(res, enriched = i)
#   pul_res <- rbind(pul_res, res)
# }
# pul_res <- merge(pul_res, select(pul, id, substrate, direction), by = "id", all.x = T)
# write.table(pul_res, "CAZy_lefse_PUL.txt", quote = F, sep = "\t", row.names = F)
# 
# source("F:/Code/R_func/plot_pie.R")
# lefse_res <- read.delim("CAZy_lefse_dat.txt") %>% filter(LDA > 2) 
# 
# group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
# group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)
# 
# table(lefse_res$enriched) %>% 
#   data.frame() %>% select(sample = Var1, val = Freq) %>% 
#   mutate(sample = factor(sample, group_order)) %>% 
#   ggpubr::ggbarplot(., x = "sample", y = "val", fill = "sample", xlab = "", ylab = "CAZymes count") +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_manual(values = group_color) +
#   theme(legend.position = "none")
# ggsave("CAZy_lefse_CAZymes_count_barplot.pdf", width = 4, height = 4)
# 
# pul_res <- read.delim("CAZy_lefse_PUL.txt")
# pul_res %>% select(id, enriched) %>% unique.data.frame() %>% group_by(enriched) %>% summarise(n = n()) %>% 
#   ungroup() %>% add_row(enriched = "d14", n = 0) %>% mutate(enriched = factor(enriched, group_order)) %>% 
#   ggpubr::ggbarplot(., x = "enriched", y = "n", fill = "enriched", xlab = "", ylab = "PULs count") +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_fill_manual(values = group_color) +
#   theme(legend.position = "none")
# ggsave("CAZy_lefse_PULs_count_barplot.pdf", width = 4, height = 4)
# 
# pul_res %>% select(enriched, id, substrate) %>% unique.data.frame() %>% select(-id) %>% group_by(enriched, substrate) %>% 
#   summarise(n = n()) %>% ungroup() %>% data.frame() %>% filter(enriched == "d5") %>% select(-enriched) %>% 
#   plot_pie(., dat_colnames = c(name = "substrate"), top_n = 10, add_n = T)
# ggsave("CAZy_lefse_PULs_substrate_d5_pie.pdf", width = 3, height = 3)
# 
