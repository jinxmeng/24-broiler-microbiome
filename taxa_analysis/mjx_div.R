#### info ####
# Jinxin Meng, 2023-10-10, 2023-11-12
pacman::p_load(tibble, dplyr, tidyr, ggplot2, vegan, ggpmisc, ggpubr)
setwd("/share/data1/mjx/proj/01.broiler_metagenome_20230301/analysis/taxa/")
source("/share/data1/mjx/Rproj/script/diversity.R")
source("/share/data1/mjx/Rproj/script/plot_PCoA.R")

#### α-diversity ####
group <- read.delim("../data/sample_group")
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% data.frame()

group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2")

method = "shannon"
plotdat <- calu_alpha(otu, method = method) %>% 
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
  stat_smooth(aes(group = 1), level = 0.95, method = "lm", color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1.5, label.x.npc = 1, label.y.npc = 0, color = "#8da0cb") +
  stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5","d7")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5","d7")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1.5, label.x.npc = 1, label.y.npc = .03, color = "#fc8d62") +
  stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5","d7")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5","d7")), 
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

plotdat %>% group_by(group) %>% summarise(mean = mean(index), sd = sd(index), count = n()) %>% 
  write.table(., paste0("div_", method, "_dat.txt"), sep = "\t", row.names = F, quote = F)
write("", paste0("div_", method, "_dat.txt"), append = T)
diff_test(select(plotdat, -group), group, method =  "wilcox", group_order = group_order) %>% 
  mutate(plab = ifelse(pval < 0.0001, "****", ifelse(pval < 0.001, "***", ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "!"))))) %>% 
  write.table(., paste0("div_", method, "_dat.txt"), sep = "\t", append = T, row.names = F, quote = F)
write("", paste0("div_", method, "_dat.txt"), append = T)
write.table(plotdat, paste0("div_", method, "_dat.txt"), sep = "\t", append = T, row.names = F, quote = F)

#### PCoA ####
p <- plot_PCoA(select(otu, -d3O), group, group_order = group_order, group_color = group_color, dis_method = "bray")
ggsave(file = "div_PCoA_scatterplot.pdf", width = 6, height = 4.5)
write.table(p$data, "div_PCoA_dat.txt", sep = "\t", row.names = T, quote = F)

#### adonis ####
dat <- calu_adonis_pair(otu, group, group_order = group_order, dis_method = "bray", permutations = 999)
plotdat <- data.frame(x = sapply(strsplit(dat$group_pair, split = "_vs_"), "[[", 1),
                      y = sapply(strsplit(dat$group_pair, split = "_vs_"), "[[", 2),
                      r2 = dat$r2, pval = dat$pval) %>% 
  mutate(x = factor(x, rev(group_order)),
         y = factor(y, group_order),
         plab = ifelse(pval<=0.001, "p<=0.001", 
                       ifelse(pval<0.01, "p<0.01", 
                              ifelse(pval<0.05, "p<0.05", "p>=0.05"))),
         plab = factor(plab, levels = c("p<=0.001", "p<0.01", "p<0.05", "p>=0.05")))
write.table(plotdat, "div_adonis_dat.txt", sep = "\t", row.names = T, quote = F)

ggplot(plotdat, aes(x, y)) +
  geom_tile(fill = "transparent", color = "black", width = 1, height = 1, lwd = .4) +
  # geom_tile(aes(width = r2, height = r2, fill = plab), lwd = .4, color = "#000000") +
  geom_point(aes(size = r2, fill = plab), shape = 21, color = "#000000", stroke = .4) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = c("#fb6a4a", "#fc9272", "#fee0d2", "#ffffff")) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        plot.title = element_text(size = 10, color = "black"),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        aspect.ratio = 1)
ggsave(file = "div_adonis_scatterplot.pdf", width = 4.5, height = 4.5)