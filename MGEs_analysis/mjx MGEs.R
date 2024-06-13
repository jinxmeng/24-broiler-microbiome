#### info ####
# Jinxin Meng, 20230606, 20231223
pacman::p_load(tidyr, dplyr, tibble, purrr)
pacman::p_load(ggplot2, ggrepel, ggforce, ggpubr, ggpmisc)
setwd("F:/Microbiome/proj_2023/02.Broiler_fecal_resistome_20230302/MGEs/")

#### MGEs drug statistics ####
source("F:/Code/R_func/plot_pie.R")
dat <- read.delim("clipboard")
plot_pie(dat, top_n = 12, fill = "auto")
ggsave("MGEs_class_pieplot.pdf", width = 3, height = 3)

#### prevalence and abundance of MGEs ####
# 计算流行率，tpm<10的不纳入计算，计算丰度就不进行过滤了
source("F:/Code/R_func/profile_process.R")
metadata <- read.delim("MGEs_metadata.txt")
profile <- read.delim("MGEs_tpm", row.names = 1)
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

# 流行率和丰度
prevalence <- profile_replace(profile, 10) %>% profile_prevalence()
abundance <- data.frame(avg_tpm = rowMeans(profile)) %>% rownames_to_column(var = "name")
dat <- merge(prevalence, abundance, by = "name") %>% 
  merge(., y = select(metadata, MGE_name, type) %>% unique.data.frame(), 
        by.x = "name", by.y = "MGE_name")
write.table(dat, "MGEs_prevalence.txt", sep = "\t", quote = F, row.names = F)

plotdat <- mutate(dat, avg_tpm = log10(avg_tpm+1e-3))

type_color <- structure(c("#1b9e77","#d95f02","#7570b3","#66a61e","#e6ab02"), names = c("Transposon","Plasmid","Insertion_element","Integron","Transposase"))
ggplot(plotdat, aes(x = prevalence, y = avg_tpm)) +
  geom_vline(xintercept = 80, lty = 2, lwd = .5) +
  geom_hline(yintercept = log10(10), lty = 2, lwd = .5) +
  geom_point(aes(fill = type), size = 2, shape = 21, stroke = .3) +
  scale_fill_manual(values = type_color) +
  labs(x = "Prevalence (%)", y = "Average TPM") +
  geom_text(data = filter(plotdat, prevalence > 30 | avg_tpm > log10(30)), aes(label = name), size = 1) +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(size = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)
ggsave("MGEs_prevalence_scatterplot.pdf", height = 4, width = 5)

#### MGEs presence/absence heatmap ####
# overall
library(ComplexHeatmap)
source("F:/Code/R_func/profile_process.R")
metadata <- read.delim("MGEs_metadata.txt")
profile <- read.delim("MGEs_tpm", row.names = 1)
dat <- profile_replace(profile, 10) %>% profile_adjacency()
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

annotation_row <- data.frame(name = rownames(dat)) %>% 
  mutate(type = metadata$type[match(name, metadata$MGE_name)]) %>% column_to_rownames(var = "name")
type_color <- structure(c("#1b9e77","#d95f02","#7570b3","#66a61e","#e6ab02"), names = sort(unique(as.character(annotation_row$type))))
annotation_col <- data.frame(name = colnames(dat)) %>% 
  mutate(group = group$group[match(name, group$sample)]) %>% 
  column_to_rownames(var = "name")
group_color <- structure(group_color, names = group_order)
annotation_colors <- list(group = group_color, type = type_color)
col_split <- factor(annotation_col$group, group_order)

pdf("MGEs_prevalence_overall_heatmap.pdf", width = 10, height = 12)
pheatmap(dat, color = circlize::colorRamp2(c(0, 1), c("white","#74add1")),
         cluster_rows = T, cluster_cols = T, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, column_gap = unit(1,'mm'), border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 3, fontsize_col = 0,
         cellwidth = 4, cellheight = 4, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

# partial
metadata <- read.delim("MGEs_metadata.txt")
profile <- read.delim("MGEs_tpm", row.names = 1)
dat <- profile_replace(profile, 10) %>% profile_adjacency()
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"),names = group_order)

prevalence <- read.delim("MGEs_prevalence.txt") %>% filter(prevalence > 10)
dat <- data.frame(dat[as.character(prevalence$name),])

annotation_row <- data.frame(name = rownames(dat)) %>% 
  mutate(type = metadata$type[match(name, metadata$MGE_name)]) %>% column_to_rownames(var = "name")
type_color <- structure(c("#1b9e77","#d95f02","#7570b3","#66a61e","#e6ab02"), names = sort(unique(as.character(annotation_row$type))))
annotation_col <- data.frame(name = colnames(dat)) %>% 
  mutate(group = group$group[match(name, group$sample)]) %>% 
  column_to_rownames(var = "name")
group_color <- structure(group_color, names = group_order)
annotation_colors <- list(group = group_color, type = type_color)
col_split <- factor(annotation_col$group, group_order)

pdf("MGEs_prevalence_partial_heatmap.pdf", width = 12, height = 6)
pheatmap(dat, color = circlize::colorRamp2(c(0, 1), c("white","#74add1")),
         cluster_rows = T, cluster_cols = T, treeheight_col = 20, treeheight_row = 20,
         column_split = col_split, column_gap = unit(1,'mm'), border_gp = gpar(col = "black"),
         show_rownames = T, show_colnames = F, fontsize_row = 5, fontsize_col = 0,
         cellwidth = 5, cellheight = 5, use_raster = F,
         annotation_row = annotation_row, annotation_col = annotation_col, annotation_colors = annotation_colors)
dev.off()

#### upset plot ####
library(UpSetR)
source("F:/Code/R_func/profile_process.R")
profile <- read.delim("MGEs_tpm", row.names = 1) %>% profile_replace(10)
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- structure(c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2"), names = group_order)

plotdat <- profile_smp2grp(profile, group) %>% profile_adjacency()
upset_metadata <- data.frame(id = group_order, set = group_order)

pdf("MGEs_shared_upset.pdf", width = 6, height = 4)
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
rownames_to_column(plotdat, var = "name") %>% write.table("MGEs_shared_upset_dat.txt", sep = "\t", quote = F, row.names = F)

#### overall abundance ####
source("F:/Code/R_func/diversity.R")
profile <- read.delim("MGEs_tpm", row.names = 1) %>% apply(., 2, \(x) x/1e6*100) %>% data.frame()
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
ggsave("div_overall_abundance_boxplot.pdf", height = 4, width = 3)

#### α-diversity ####
source("F:/Code/R_func/diversity.R")
source("F:/Code/R_func/profile_process.R")
profile <- read.delim("MGEs_tpm", row.names = 1) %>% profile_replace(10)
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
  stat_smooth(aes(group = 1), method = "loess", color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1.5, label.x.npc = 1, label.y.npc = .03, color = "#fc8d62") +
  stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5")), 
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
ggsave(paste0("div_", method, "_boxplot.pdf"), height = 4, width = 3)

plotdat %>% group_by(group) %>% summarise(mean = mean(val), sd = sd(val), count = n()) %>% 
  write.table(., paste0("div_", method, "_dat.txt"), sep = "\t", row.names = F, quote = F)
write("", paste0("div_", method, "_dat.txt"), append = T)
diff_test(select(plotdat, -group), group, method =  "wilcox", group_order = group_order) %>% 
  mutate(plab = ifelse(pval < 0.0001, "****", ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "!"))))) %>% 
  write.table(., paste0("div_", method, "_dat.txt"), sep = "\t", append = T, row.names = F, quote = F)
write("", paste0("div_", method, "_dat.txt"), append = T)
write.table(plotdat, paste0("div_", method, "_dat.txt"), sep = "\t", append = T, row.names = F, quote = F)

#### PCoA ####
source("F:/Code/R_func/plot_PCoA.R")
p <- plot_PCoA(profile, group, group_order = group_order, group_color = group_color)
ggsave("div_PCoA_scatterplot.pdf", height = 4.5, width = 6)
write.table(p$data, "div_PCoA_dat.txt", sep = "\t", quote = F, row.names = F)

#### corr div ####
div <- read.delim("clipboard")
ggscatter(data = div, x = "arg_richness", y = "mge_richness", rug = T, 
          add = "reg.line", conf.int = T, conf.int.level = .95, 
          xlab = "ARGs richness", ylab = "MGEs richness", size = 2, 
          add.params = list(color = "#f2390c", fill = "lightgray", size = 1)) +
  stat_cor(label.sep = "\n", color = "black", label.x = 100, 
           label.y = 60, method = "spearman", size = 3) +
  theme(aspect.ratio = 1)
ggsave("corr_div_richness_scatterplot.pdf", width = 4, height = 4)

ggscatter(data = div, x = "arg_shannon", y = "mge_shannon", rug = T, 
          add = "reg.line", conf.int = T, conf.int.level = .95, 
          xlab = "ARGs shannon", ylab = "MGEs shannon", size = 2, 
          add.params = list(color = "#f2390c", fill = "lightgray", size = 1)) +
  stat_cor(label.sep = "\n", color = "black", label.x = 4, 
           label.y = 1, method = "spearman", size = 3) +
  theme(aspect.ratio = 1)
ggsave("corr_div_shannon_scatterplot.pdf", width = 4, height = 4)

ggscatter(data = div, x = "arg_abund", y = "mge_abund", rug = T, 
          add = "reg.line", conf.int = T, conf.int.level = .95, 
          xlab = "ARGs overall abundance", ylab = "MGEs overall abundance", size = 2, 
          add.params = list(color = "#f2390c", fill = "lightgray", size = 1)) +
  stat_cor(label.sep = "\n", color = "black", label.x = 0.1, 
           label.y = 2.5, method = "spearman", size = 3) +
  theme(aspect.ratio = 1)
ggsave("corr_div_overall_scatterplot.pdf", width = 4, height = 4)


#### network ####
source("F:/Code/R_func/taxa.R")
source("F:/Code/R_func/profile_process.R")
source("F:/Code/R_func/corr_process.R")
metadata_args <- read.delim("../ARGs/ARGs_metadata.txt")
metadata_mges <- read.delim("MGEs_metadata.txt")

# 位置关系
# 以ARG为中心去推断两侧是否携带MGE
# 所以如果ARG本身就是MGE，那么这个contig就被保留了
# 有可能这条contig上仅有这一个基因，同时作为ARG和MGE
# 统计这种关系的数量，要考虑全，一个基因本身既是ARG也是MGE，就算是一组了，另外，这个ARG周围还有别的ARG的话，那么关系数量+1

gene_info_args <- read.delim("gene_info_args.tsv", header = F, col.names = c("name", "gene", "gene_name"))
gene_info_mges <- read.delim("gene_info_mges.txt", header = F, col.names = c("name", "gene", "gene_name"))

colinearity <- read.delim("total_find.tsv", header = F, col.names = c("name", "gene", "type", "start", "end", "dir")) %>% 
  filter(type != "Other") %>% select(name, gene, type) %>% 
  mutate(gene_name = ifelse(type %in% c("ARGs", "ARGs|MGEs"), 
                            gene_info_args$gene_name[match(gene, gene_info_args$gene)], 
                            gene_info_mges$gene_name[match(gene, gene_info_mges$gene)]),
         gene_name = ifelse(type %in% c("ARGs", "ARGs|MGEs"), 
                            metadata_args$ARO_short_name[match(gene_name, metadata_args$ARO_accession)],
                            gene_name)
         )

tmp_group <- unique(colinearity$name) 
tmp_pair <- data.frame(name_x = character(), name_y = character())
for(i in tmp_group){
  tmp_dat <- dplyr::filter(colinearity, name == i)
  if(sum(tmp_dat$type=="ARGs|MGEs") > 0){
    x <- tmp_dat[tmp_dat$type=="ARGs|MGEs",]
    x[1,3] <- "MGEs"
    y = gene_info[gene_info$gene %in% x$name, 3]
    x[1,4] = y[!grepl("^300", y)]
    tmp_dat[tmp_dat$type=="ARGs|MGEs",3] <- "ARGs"
    tmp_dat <- rbind(tmp_dat, x)
  }
  args <- tmp_dat[tmp_dat$type=="ARGs", 4]
  mges <- tmp_dat[tmp_dat$type=="MGEs", 4]
  for(j in mges){
    tmp_pair <- tmp_pair %>% add_row(name_x = args, name_y = j)
  }
}

plotdat <- group_by(tmp_pair, name_x, name_y) %>% summarise(n = n()) %>% arrange(desc(n))
write.table(plotdat, "total_co-occurence_pair.tsv", sep = "\t", quote = F, row.names = F)

# 网络图
library(igraph)
library(ggraph)
library(tidygraph)
edges <- select(plotdat, name_x, name_y, n) %>%
  rename(from = name_x, to = name_y)

nodes <-rbind(data.frame(name = unique(edges$from)) %>% 
                add_column(type = "ARG") %>% 
                mutate(class = metadata_args$drug3[match(name, metadata_args$ARO_short_name)],
                       class = forcats::fct_lump_n(class, 9, ties.method = "last")),
              data.frame(name = unique(edges$to)) %>% 
                add_column(type = "MGE") %>% 
                mutate(class = metadata_mges$type[match(name, metadata_mges$MGE_name)],
                       class = forcats::fct_lump_n(class, 9, ties.method = "last"))) %>% 
  mutate(type = factor(type, labels = c("ARG", "MGE")))

graph <- tbl_graph(nodes = nodes, edges = edges)
nodes$degree <- degree(graph)                  #边的数量给点赋一个值
graph <- tbl_graph(nodes = nodes, edges = edges)

class_color <- c("#66c2a5","#fc8d62","#8da0cb","#a6d854","#ffd92f",
                 "#fb8072","#80b1d3","#bc80bd","#e5c494","#b3b3b3",
                 "#1b9e77","#d95f02","#7570b3","#66a61e","#e6ab02")

setseed(2023)
ggraph(graph, layout = 'fr') + 
  geom_edge_arc(aes(edge_width = n), color = "grey70", strength = .15, show.legend = T) + 
  scale_edge_width(range = c(.2, 1)) +
  scale_edge_alpha_manual(values = c(1, 2)) +
  geom_node_point(aes(fill = class, size = degree), shape = 21, stroke = .3) +
  scale_size(range = c(2, 4)) +
  scale_fill_manual(values = class_color) +
  geom_node_text(aes(label = name), size = .6, fontface = "italic") +
  theme(panel.background = element_rect(fill = NA, color = "black"),
        legend.key = element_rect(fill = NA),
        aspect.ratio = 1)
ggsave("corr_cooccurence.pdf", width = 7, height = 7)

#### gene arrow plot ####
library(gggenes)
source("F:/Code/R_func/taxa.R")
tmp_vec <- unique(colinearity$name)
tmp_list_gene_name <- list()
tmp_list_gene <- list()
for (i in tmp_vec) tmp_list_gene_name[[i]] <- filter(colinearity, name == i) %>% select(gene_name) %>% unlist() %>% as.character()
for (i in tmp_vec) tmp_list_gene[[i]] <- filter(colinearity, name == i) %>% select(gene) %>% unlist() %>% as.character()

show_line <- plotdat %>% head(21)
show_line <- show_line[-1,]
x = c()
y = c()
for (i in 1:nrow(show_line)) {
  tmp_x <- show_line[i,1:2] %>% as.character()
  for (j in names(tmp_list_gene_name)) {
    if (all(tmp_x %in% tmp_list_gene_name[[j]]) & length(tmp_list_gene_name[[j]]) == 2 & !j %in% x) {
      x <- c(x, j)
      y <- c(y, paste0(tmp_x, collapse = "-"))
      break
    } else if (all(tmp_x %in% tmp_list_gene_name[[j]]) & length(tmp_list_gene_name[[j]]) == 3 & !j %in% x) {
      x <- c(x, j)
      y <- c(y, paste0(tmp_x, collapse = "-"))
      break
    }
  }
}

tmp_dat <- data.frame(name = x, lab = y)

dat <- read.delim("total_find.tsv", header = F, col.names = c("name", "gene", "class", "start", "end", "direct")) %>% 
  mutate(direct = ifelse(direct == "+", TRUE, FALSE),
         contig = gsub("_\\d+$", "", name)) %>% 
  filter(name %in% x) %>% 
  mutate(class = ifelse(class == "ARGs|MGEs", "ARGs", class)) %>% 
  merge(., tmp_dat, by = "name")


ggplot(dat, aes(xmin = start, xmax = end, y = lab, fill = class, forward = direct)) +
  geom_gene_arrow(arrowhead_height = unit(3.5, "mm"), arrowhead_width = unit(2, "mm"),
                  arrow_body_height = unit(2, "mm"), color = "black", size = .2) +
  facet_wrap(~factor(lab), scale = "free_y", ncol = 1) +
  scale_fill_manual(values = c("#f8766d", "#00ba38", "#ffffff")) +
  theme_genes() +
  theme(panel.grid.major.y = element_line(color = "grey70", size = .5),
        plot.caption = element_text(size = 6, color = "black", hjust = 0, lineheight = 1.2),
        plot.caption.position = "panel",
        plot.subtitle = element_text(size = 6, color = "black", lineheight = 1.2),
        axis.line.x = element_line(linewidth = .4, color = "black"),
        axis.ticks.x = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        # axis.text.y = element_blank(),
        axis.title.x = element_text(size = 8, color = "black", hjust = 1))
ggsave("corr_cooccurence_arrowplot_partial2.pdf", width = 6, height = 10)


#### co-abundance ####
# 相关性
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
profile_args <- read.delim("../ARGs/ARGs_tpm", row.names = 1) %>%
  profile_filter(min_abundance = 1, min_prevalence = 0.05, by_group = F) %>%
  relocate(group$sample) %>% t() %>% data.frame(check.names = F)
profile_mges <- read.delim("MGEs_tpm", row.names = 1) %>%
  profile_filter(min_abundance = 1, min_prevalence = 0.05, by_group = F) %>%
  relocate(group$sample) %>% t() %>% data.frame(check.names = F)

# 相关性分析
corr <- psych::corr.test(x = profile_args, y = profile_mges, method = "spearman", adjust = "BH")
corr_res <- corr_process(corr, rho = 0.5, padj = 0.05, out_df = T) %>%
  mutate(direct = ifelse(r > 0, "pos", "neg"))
corr_res$name_x <- metadata_args$ARO_short_name[match(corr_res$name_x, metadata_args$ARO_accession)]
write.table(corr_res, "corr_res.tsv", sep = "\t", row.names = F, quote = F)

# 整合相关性结果和位置关系
index = c()
for(i in 1:nrow(corr_res)){
  tmp_vec <- as.character(corr_res[i,1:2])
  a = 0
  for (j in tmp_list) {
    if(all(tmp_vec %in% j)){
      a = 1
    }
  }
  if (a == 1){
    index = c(index, i)
  }
}
