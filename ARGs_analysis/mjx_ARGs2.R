#### info ####
# Jinxin Meng, 20231225, 20231228
pacman::p_load(tidyr, dplyr, tibble, purrr)
setwd("F:/Proj/proj_2023/02.Broiler_fecal_resistome_20230302/ARGs_related_to_taxa_and_mges/")

#### div cor ####
source("F:/Code/R_func/diversity.R")
profile_args <- read.delim("../ARGs/ARGs_tpm", row.names = 1)
profile_mges <- read.delim("../MGEs/MGEs_tpm", row.names = 1)
profile_taxa <- read.delim("../taxa/species_tpm", row.names = 1)

dat <- list(
  args = calu_alpha(profile_args, method = "shannon", out_colnames = "args"),
  taxa = calu_alpha(profile_taxa, method = "shannon", out_colnames = "taxa"),
  mges = calu_alpha(profile_mges, method = "shannon", out_colnames = "mges")
  ) %>% reduce(\(x, y) merge(x, y, by = "sample"))

grp_pair <- list(c("args", "mges"), c("args", "taxa"))
p_list <- list()
x = 2
p_list[[x]] <- ggscatter(dat, x = grp_pair[[x]][1], y = grp_pair[[x]][2], rug = T, add = "reg.line",
                         conf.int = T, conf.int.level = .95, title = "shannon",
                         xlab = grp_pair[[x]][1], ylab = grp_pair[[x]][2], size = .8,
                         add.params = list(color = "#0ea4b0", fill = "lightgray", size = .5)) +
  stat_cor(label.sep = "\n", color = "black", method = "spearman", size = 3) +
  theme(aspect.ratio = 1)
cowplot::plot_grid(plotlist = p_list)
ggsave("div_cor_shannon.pdf", width = 6, height = 3)


#### spearman ####
source("F:/Code/R_func/corr_process.R")
source("F:/Code/R_func/taxa.R")
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
prevalence_args <- read.delim("../ARGs/ARGs_prevalence.txt") %>% filter(prevalence > 20)
profile_args <- read.delim("../ARGs/ARGs_tpm", row.names = 1)
prevalence_mges <- read.delim("../MGEs/MGEs_prevalence.txt")
profile_mges <- read.delim("../MGEs/MGEs_tpm", row.names = 1)
taxonomy <- read.delim("../taxa/species_taxonomy_rename.tsv")
profile_taxa <- read.delim("../taxa/species_tpm", row.names = 1) %>% taxa_trans(taxonomy, to = "family", out_all = T)

corr <- psych::corr.test(x = prevalence_args, y = profile_taxa, method = "spearman", adjust = "BH")
corr_res <- corr_process(corr, rho = 0.5, padj = 0.05, out_df = T)

#### 普鲁克分析 ####
library(vegan)
profile_args <- read.delim("../ARGs/ARGs_tpm", row.names = 1)
profile_mges <- read.delim("../MGEs/MGEs_tpm", row.names = 1)
profile_taxa <- read.delim("../taxa/species_tpm", row.names = 1)

# 距离
dist_1 <- vegdist(t(profile_args))
dist_2 <- vegdist(t(profile_taxa))

# 降维
PCoA_1 <- cmdscale(dist_1)
PCoA_2 <- cmdscale(dist_2)

# 以对称模式进行普氏分析（symmetric = TRUE）
proc <- procrustes(PCoA_1, PCoA_2, symmetric = T)
summary(proc)

# 评价一致性
# 如果样本中物种与环境一致性（相似性）越近，则对应的残差越小，
# 反之物种与环境的相似性越远，则残差越大（三条辅助线对应的位置分别为残差25%、50%和75%）
plot(proc, kind = 2)
residuals(proc)

# 事后检验999次
# 普氏分析中M2统计量的显著性检验
# 在 Procrustes 分析中，M2 常常指代的是 Procrustes 统计量，它是度量两个形状间差异程度的一个指标。
# 更具体地说，M2 是源数据集的点经过旋转、缩放和/或平移后与目标数据集中对应点的平方距离和。
# 它代表了变换后源数据集的点与目标数据集点之间的不匹配程度。M2 的数值越小，表明两组数据集的形状越相似；反之，则表明它们之间的差异越大。
set.seed(2024)
proc_test <- protest(PCoA_1, PCoA_2, permutations = 999)
proc_test
# Call:
#   protest(X = PCoA_args, Y = PCoA_taxa, permutations = 999) 
# Procrustes Sum of Squares (m12 squared):        0.5268  # M2统计值
# Correlation in a symmetric Procrustes rotation: 0.6879 
# Significance:  0.001  # 结果显示: p = 0.001具有很强显著性
# Permutation: free
# Number of permutations: 999

# 提取结果
#偏差平方和（M2统计量）
proc_test$ss
# 对应p值结果
proc_test$signif

# 基于 ggplot2 包ggplot函数将其结果绘制成图，并利用 export 包导出 ppt 。
# 首先提取降维后的数据轴1和2的坐标，并且提取转换的坐标；然后进行绘制。
# 获得x和y轴的坐标及旋转过的坐标
# Pro_Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
# Pro_X <- data.frame(proc$rotation)
proc_point <- cbind(
  data.frame(proc$Yrot) %>% rename(X1_rotated = X1, X2_rotated = X2),# Y-矩阵旋转后的坐标，就是物种矩阵为了逼近X做了调整，调整后的坐标。
  data.frame(proc$X) %>% rename(X1_target = Dim1, X2_target = Dim2) # X-目标矩阵，就是proc分析输入的第一个矩阵，作为目标矩阵没变化。
  )
proc_coord <- data.frame(proc$rotation) # Y旋转后的坐标轴

# 绘图
ggplot(proc_point) + # 旋转坐标到目的坐标的一半上一个颜色，另一半上不同的颜色
  geom_segment(aes(x = X1_rotated, y = X2_rotated, xend = (X1_rotated + X1_target)/2, yend = (X2_rotated + X2_target)/2), 
               arrow = arrow(length = unit(0, 'cm')), color = "#9BBB59", size = .4) +
  geom_segment(aes(x = (X1_rotated + X1_target)/2, y = (X2_rotated + X2_target)/2, xend = X1_target, yend = X2_target),
               arrow = arrow(length = unit(0.2, 'cm')), color = "#957DB1", size = .4) +
  geom_point(aes(X1_rotated, X2_rotated), color = "#9BBB59", size = 1.6, shape = 16) + # 旋转坐标点
  geom_point(aes(X1_target, X2_target), color = "#957DB1", size = 1.6, shape = 16) + # 目的坐标点
  labs(x = 'Dim 1', y = 'Dim 2',
       subtitle = paste0("coefficients: M2 = ", round(proc_test$ss, 4), ", p = ", proc_test$signif)) +
  labs(title = "Correlation analysis by Procrustes analysis") +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.4) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.4) +
  geom_abline(intercept = 0, slope = proc_coord[1,2]/proc_coord[1,1], size = 0.4) +
  geom_abline(intercept = 0, slope = proc_coord[2,2]/proc_coord[2,1], size = 0.4) +
  theme_bw() +
  theme(axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.line = element_blank(),
        plot.title = element_text(size = 10, color = "black"),
        plot.subtitle = element_text(size = 10, color = "black"),
        panel.border = element_rect(linewidth = .4, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"),
        aspect.ratio = 3/4)
ggsave("procrustes_corr_analysis.pdf", width = 6, height = 4.5)

#### mantel test for taxa and ARGs phenotype ####
library(linkET)
source("F:/Code/R_func/taxa.R")
group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
metadata_args <- read.delim("../ARGs/ARGs_metadata.txt")
profile_args <- read.delim("../ARGs/ARGs_tpm", row.names = 1) %>% 
  mutate(drug = metadata_args$drug3[match(rownames(.), metadata_args$ARO_accession)]) %>% 
  group_by(drug) %>% summarise_all(sum) %>% ungroup() %>% column_to_rownames(var = "drug") %>% 
  select(all_of(group$sample)) %>%  t() %>% data.frame(check.names = F)
profile_taxa <- read.delim("../taxa/species_tpm", row.names = 1) %>% select(all_of(group$sample)) %>%
  t() %>% data.frame(check.names = F)
taxonomy <- read.delim("../taxa/species_taxonomy_rename.tsv")

tmp_df <- data.frame(taxa = taxonomy$family[match(colnames(profile_taxa), taxonomy$feature)],
                     ncol = 1:ncol(profile_taxa)) %>% arrange(taxa)
tmp_list <- list()
for (i in unique(tmp_df$taxa)) {
  tmp_list[[i]] = c(tmp_df$ncol[tmp_df$taxa %in% i])
}

mantel <- mantel_test(profile_taxa, profile_args, spec_select = tmp_list)
# write.table(mantel, "mantel_test_args_family.tsv", sep = "\t", row.names = F, quote = F)
mantel <- read.delim("mantel_test_args_family.tsv")
mantel_df <- mutate(mantel, 
                    rd = cut(r, breaks = c(-Inf, 0.3, 0.5, 0.7, Inf), labels = c("< 0.3", "0.3 - 0.6", "0.6 - 0.9", ">= 0.9")),
                    pd = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"))) %>% 
  filter(r > 0.3 & p < 0.05) %>% mutate(spec = forcats::fct_lump_min(spec, 2)) %>% filter(spec != "Other")

qcorrplot(correlate(profile_args), type = "upper", diag = F) +
  geom_square() +
  geom_couple(data = mantel_df, aes(colour = pd, size = rd), curvature = nice_curvature()) +
  scale_fill_viridis_c() +
  scale_size_manual(values = c(0.6, 0.9, 1.2)) +
  scale_colour_manual(values = c("#b1bb52", "#9794a4", "#aaaba6")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
ggsave("mantel_test_args_family_linkET.pdf", width = 10, height = 8)

# genus
tmp_df <- data.frame(family = taxonomy$family[match(colnames(profile_taxa), taxonomy$feature)],
                     taxa = taxonomy$species[match(colnames(profile_taxa), taxonomy$feature)],
                     ncol = 1:ncol(profile_taxa)) %>%  
  filter(family %in% c("f__Enterobacteriaceae", "f__Enterococcaceae", "f__Pseudomonadaceae", "f__Clostridiaceae"))

tmp_list <- list()
for (i in unique(tmp_df$taxa)) {
  tmp_list[[i]] = c(tmp_df$ncol[tmp_df$taxa %in% i])
}

mantel <- mantel_test(profile_taxa, profile_args, spec_select = tmp_list)
# write.table(mantel, "mantel_test_args_species.tsv", sep = "\t", row.names = F, quote = F)
mantel <- read.delim("mantel_test_args_species.tsv")
mantel_df <- mutate(mantel, 
                    rd = cut(r, breaks = c(-Inf, 0.3, 0.5, 0.7, Inf), labels = c("< 0.3", "0.3 - 0.6", "0.6 - 0.9", ">= 0.9")),
                    pd = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), labels = c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", ">= 0.05"))) %>% 
  filter(r > 0.3 & p < 0.05)

qcorrplot(correlate(profile_args), type = "upper", diag = F) +
  geom_square() +
  geom_couple(data = mantel_df, aes(colour = pd, size = rd), curvature = nice_curvature()) +
  scale_fill_viridis_c() +
  scale_size_manual(values = c(0.6, 0.9, 1.2)) +
  scale_colour_manual(values = c("#b1bb52", "#9794a4", "#aaaba6")) +
  guides(size = guide_legend(title = "Mantel's r", override.aes = list(colour = "grey35"), order = 2),
         colour = guide_legend(title = "Mantel's p", override.aes = list(size = 3), order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))
ggsave("mantel_test_args_species_linkET.pdf", width = 10, height = 8)

