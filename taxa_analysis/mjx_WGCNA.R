#### info ####
# Jinxin Meng, 2023-10-09, 2023-11-11
pacman::p_load(tibble, dplyr, tidyr, ggplot2, vegan, ggpmisc, ggpubr)
setwd("F:/Microbiome/proj_2023/02.Broiler_fecal_resistome_20230302/taxa/")
source("F:/Code/R_func/taxa.R")

#### WGCNA ####
library(WGCNA)
source("F:/Code/R_func/taxa.R")
group <- read.delim("../data/sample_group.txt")
profile <- read.delim("species_tpm", row.names = 1)%>% 
  taxa_trans(taxonomy, to = "species", out_all = T, transRA = T) %>% t()


# 确定软阈值
softPower = 6

# 邻接矩阵
adjacency = adjacency((profile), power = softPower)

# TOM矩阵
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM

# 层次聚类
tree = hclust(as.dist(dissTOM), method = "average")
plot(tree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", 
     labels = FALSE, hang = 0.04) 

# 拆解树，分析模块 层级聚类树展示各个模块
# minModuleSize规定最小的模块不小于30
# dendro = geneTree,  distM = dissTOM 输入的是前面的树和TOM矩阵
# deepSplit 是切割强度，越大得到的模块越多
# pamStage = TRUE, pamRespectsDendro = FALSE,
# 这是对PAM方法的限制，pamStage = TRUE是使用PAM方法，pamRespectsDendro = FALSE，是说PAM不会考虑树的感受
# 翻译就是，如果设置为TRUE，那么PAM在分配模棱两可基因的时候会把他们限定在同一个分支中，这是残血版
# 如果设置为FALSE，PAM方法就是满血版本了，按照自己的方法去分配。
# 结果就是，有很多灰色的不属于任何模块的基因，也找到了归属，有时候会导致结果不好解释。
# 该函数最终返回的结果是各个基因的模块编号
# 总共找到22个模块，0里面的88个基因不属于任何一个模块，就是垃圾箱。
# labels2colors这个函数可以把数字转换为颜色。
minModuleSize = 15
dynamicMods = cutreeDynamic(dendro = tree, distM = dissTOM, deepSplit = 2, pamStage = T, pamRespectsDendro = F,
                            minClusterSize = minModuleSize)
table(dynamicMods)
color <- names(dynamicMods)
color[color==""] = 0
dynamicColors = paletteer::paletteer_d("ggthemes::Tableau_20")[as.numeric(color)+1] %>% as.character()

# Convert numeric lables into colors 
# dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors) 

# Plot the dendrogram and colors underneath 
# pdf(file="5_Dynamic Tree Cut.pdf",width=8,height=6) 
plotDendroAndColors(tree, dynamicColors, "Dynamic Tree Cut", 
                    dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, 
                    main = "Gene dendrogram and module colors") 

color = colorRampPalette(c("#4575b4","#74add1","#e0f3f8"))(100)
pdf("tmp_WCGNA.pdf", width = 7, height = 7)
TOMplot(dissTOM^5, tree, dynamicColors, col = color)
dev.off()

#### 模块重现 #### 
dat <- data.frame(t(profile)) %>% 
  add_column(module = dynamicMods)

group <- read.delim("../Data/sample_group.txt", sep = "\t", header = T, row.names = NULL)
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")

# lm
p_list <- list()
x = 1
for (i in unique(dat$module)){
  plotdat <- rownames_to_column(dat[dat$module==i,], var = "name") %>% select(-module) %>% 
    gather(key = "sample", value = "val", -name) %>%
    mutate(group = group$group[match(sample, group$sample)],
           group = factor(group, group_order), val = log10(val+1e-3))
  p <- ggplot(plotdat, aes(x = group, y = val, color = name)) +
    # geom_jitter(size = 1, width = .4, show.legend = F) +
    geom_smooth(aes(group = name), method = "loess", se = F, linewidth = 1, show.legend = F) +
    labs(x = "", y = "TPM (log10)") +
    theme_bw() +
    theme(axis.line = element_line(linewidth = .4, color = "black"),
          axis.ticks = element_line(linewidth = .4, color = "black"),
          axis.text = element_text(size = 8, color = "black"),
          axis.title = element_text(size = 8, color = "black"),
          panel.border = element_rect(linewidth = .4, color = "black"),
          panel.grid = element_blank(),
          aspect.ratio = 1)
  p_list[[x]] <- p
  x = x + 1
  }
cowplot::plot_grid(plotlist = p_list, nrow = 2)


# 模块间进行矫正合并
# 计算模块合并的高度
# Calculate eigengenes
MEList = moduleEigengenes(profile, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
# pdf(file="6_Clustering of module eigengenes.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

MEDissThres = 0.4 ###剪切高度可修改
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()

# 重新绘制合并后的模块层次聚类图
# # Call an automatic merging function 
# merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
# # The merged module colors 
# mergedColors = merge$colors 
# # Eigengenes of the new merged modules: 
# mergedMEs = merge$newMEs 
# table(mergedColors) 
# 
# #sizeGrWindow(12, 9) 
# pdf(file="7_merged dynamic.pdf", width = 9, height = 6) 
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
#                     c("Dynamic Tree Cut", "Merged dynamic"), 
#                     dendroLabels = FALSE, hang = 0.03, 
#                     addGuide = TRUE, guideHang = 0.05) 
# dev.off()


# 2.8 可视化基因网络 (TOM plot)

# plotDiss = selectTOM^7 
# diag(plotDiss) = NA
# library("gplots") 
# pdf(file="Network heatmap plot_selected genes.pdf",width=9, height=9) 
# mycol = colorpanel(250,'red','orange','lemonchiffon') 
# TOMplot(profile, tree, selectColors, col=mycol ,main = "Network heatmap plot, selected genes") 
# dev.off()
