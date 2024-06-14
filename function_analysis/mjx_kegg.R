#### info ####
# Jinxin Meng, 20231021, 20231207
pacman::p_load(tidyr, dplyr, tibble, ggplot2, ggpmisc, ggpubr)
setwd("F:/Microbiome/proj_2023/02.Broiler_fecal_resistome_20230302/kegg/")

#### KO overall abundance ####
otu <- read.delim("KO_tpm", row.names = 1)
group <- read.delim("../data/sample_group.txt")
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2")

source("F:/Code/R_func/diff_test.R")
plotdat <- data.frame(val = colSums(otu)) %>% rownames_to_column(var = "sample") %>% merge(group, by = "sample") %>% 
  mutate(group = factor(group, group_order), val = val/1000000*100)

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
  stat_smooth(data = plotdat %>% filter(group %in% c("d1","d3","d5")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = plotdat %>% filter(group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0.5, label.y.npc = 1, color = "#fc8d62") +
  stat_smooth(data = plotdat %>% filter(!group %in% c("d1","d3","d5")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = plotdat %>% filter(!group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 1, label.y.npc = 1, color = "#66c2a5") +
  scale_color_manual(values = group_color) +
  labs(x = "", title = "Overall abundance", y = "") +
  coord_cartesian(ylim = c(5.8, NA)) +
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
ggsave("KO_tpm_overall_abundance.pdf", height = 4, width = 3)

#### KO diversity ####
source("F:/Code/R_func/diversity.R")

otu <- read.delim("KO_tpm_filter_m10", row.names = 1)
group <- read.delim("../data/sample_group.txt")
group_order <- c("d1", "d3", "d5", "d7", "d14", "d21", "d28", "d35", "d42")
group_color <- c("#729ece","#ff9e4a","#67bf5c","#ed665d","#ad8bc9","#a8786e","#ed97ca","#cdcc5d","#a2a2a2")

method = "richness"
plotdat <- calu_alpha(otu, method = method) %>% 
  merge(., group, by = "sample", all.x = T) %>% 
  mutate(group = factor(group, group_order))

comparison <- diff_test(select(plotdat, -group), group, group_order = group_order, method = "wilcox") %>% 
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
  stat_smooth(data = plotdat %>% filter(group %in% c("d1","d3","d5")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = plotdat %>% filter(group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0.5, label.y.npc = 1, color = "#fc8d62") +
  stat_smooth(data = plotdat %>% filter(!group %in% c("d1","d3","d5")), 
              aes(group = 1), level = 0.95, method = "lm",
              color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = plotdat %>% filter(!group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 1, label.y.npc = 1, color = "#66c2a5") +
  scale_color_manual(values = group_color) +
  labs(x = "", title = paste0(stringr::str_to_title(method), " Index"), y = "") +
  coord_cartesian(ylim = c(5.8, NA)) +
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
p <- plot_PCoA(otu, group, group_order = group_order, group_color = group_color, dis_method = "bray", display_type = "line")
ggsave(file = "div_PCoA_scatterplot.pdf", width = 6, height = 4.5)
write.table(p$data, "div_PCoA_dat.txt", sep = "\t", row.names = T, quote = F)

#### adonis ####
dat <- calu_pairwise_adonis(otu, group, group_order = group_order, dis_method = "bray", permutations = 999)
plotdat <- data.frame(x = sapply(strsplit(dat$group_pair, split = "_vs_"), "[[", 1),
                      y = sapply(strsplit(dat$group_pair, split = "_vs_"), "[[", 2),
                      r2 = dat$r2, pval = dat$pval) %>% 
  mutate(x = factor(x, levels = rev(group_order)),
         y = factor(y, levels = (group_order)),
         plab = ifelse(pval<=0.001, "p<=0.001", ifelse(pval<0.01, "p<0.01", ifelse(pval<0.05, "p<0.05", "p>=0.05"))),
         plab = factor(plab, levels = c("p<=0.001", "p<0.01", "p<0.05", "p>=0.05")))
write.table(plotdat, "div_adonis_dat.txt", sep = "\t", row.names = T, quote = F)

ggplot(plotdat, aes(x, y)) +
  geom_tile(fill = "transparent", color = "black", width = 1, height = 1, lwd = .4) +
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
ggsave(file = "div_adnois_scatterplot.pdf", width = 4.5, height = 4.5)

#### kegg level A ####
db <- read.delim("KO_level_A_B_C_D_Description", header = F) %>% 
  select(lv1 = V1, lv1name = V2) %>% unique.data.frame() %>% 
  filter(!lv1 %in% c("A09150","A09160","A09190"))

otu <- read.delim("KO_tpm.levelA", row.names = 1) %>% 
  apply(2, \(x) x/sum(x)) %>% data.frame()

plotdat <- dplyr::filter(otu, rownames(otu) %in% db$lv1) %>%  t() %>% data.frame() %>% 
  mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
  gather(key = "lv1", value = "val", -group) %>% 
  mutate(lv1 = db$lv1name[match(lv1, db$lv1)],
         group = factor(group, group_order))

ggplot(plotdat, aes(group, val)) +
  facet_wrap(~ factor(lv1), scale = "free", nrow = 1) +
  geom_boxplot(position = position_dodge(), fill = "transparent", width = .6, 
               color = "#000000", linewidth = .3, outlier.shape = NA) +
  geom_jitter(aes(color = group), size = 1.5, width = .3, show.legend = F) +
  stat_smooth(aes(group = 1), level = 0.95, method = "lm",
              color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0, label.y.npc = 1, color = "#8da0cb") +
  stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5")), aes(group = 1), 
              level = 0.95, method = "lm", color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0.5, label.y.npc = 1, color = "#fc8d62") +
  stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5")), aes(group = 1), 
              level = 0.95, method = "lm", color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 1, label.y.npc = 1, color = "#66c2a5") +
  scale_color_manual(values = group_color) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 10, color = "black", face = "italic"),
        panel.background = element_rect(linewidth = .4, color = "black", fill = "transparent"),
        panel.grid = element_line(linewidth = .5, color = "grey85"))
ggsave("KO_tpm_not_to_RA_levelA_boxplot.pdf", width = 12, height = 3)

# diff test 
source("F:/Code/R_func/diff_test.R")
diff <- diff_test_profile(otu, group, group_order)
write.table(diff, "KO_tpm_levelA_diff_test.txt", sep = "\t", quote = F, row.names = F)

#### kegg level A phylum ####
source("F:/Code/R_func/profile_process.R")
source("F:/Code/R_func/utilities.R")

otu <- read.delim("kO_phylum_tpm", row.names = NULL)
db <- read.delim("KO_level_A_B_C_D_Description", header = F, quote = "") %>% 
  select(lv1 = V1, lv1name = V2, lv4 = V7) %>% filter(lv4 %in% otu$KO) %>% unique.data.frame()
db2 <- data.frame(unique(db[1:2]), row.names = NULL) %>% filter(!lv1 %in% c("A09150","A09160","A09190"))

# KO-->lv1 一对一
dat_x <- filter(otu, KO %in% get_freq(db$lv4, eq = 1)) %>% 
  mutate(lv1 = db$lv1[match(KO, db$lv4)]) %>% relocate(lv1)

# KO-->lv1 一对多
dat_y <- filter(otu, KO %in% get_freq(db$lv4, ne = 1))
tmp_list <- list()
pb <- txtProgressBar(style = 3)
x = 1
for (i in unique(dat_y$KO)) {
  kegg_vec <- db[db$lv4 %in% i, 1]
  tmp_dat <- filter(dat_y, KO %in% i)
  tmp_dat <- cbind(tmp_dat[,1:2], as.matrix(tmp_dat[,-(1:2)])/length(kegg_vec))
  tmp_dat <- map_dfr(kegg_vec, \(x) add_column(tmp_dat, lv1 = x, .before = "KO"))
  tmp_list[[i]] <- tmp_dat
  setTxtProgressBar(pb, x/length(unique(dat_y$KO))) 
  x = x + 1
}
close(pb)

dat_y = map_dfr(tmp_list, \(x) x)
dat_z = filter(otu, KO == "Unknown") %>% 
  add_column(lv1 = "Unknown", .before = "KO")

dat <- reduce(list(dat_x, dat_y, dat_z), \(x, y) rbind(x, y))
dat <- cbind(dat[1:3], apply(dat[-(1:3)], 2, \(x) x/sum(x))) %>% data.frame() %>% select(-KO)

# select(dat, -lv1) %>% group_by(taxa) %>% summarise_all(sum) %>% column_to_rownames(var = "taxa") %>%
#   rowMeans() %>% data.frame(taxa = .) %>% arrange(desc(taxa))
taxa_vec  <- c("Bacillota", "Unknown", "Bacteroidota", "Pseudomonadota", "Actinomycetota", 
               "Uroviricota", "Campylobacterota", "Fusobacteriota", "Verrucomicrobiota")

plotdat <- mutate(dat, taxa = ifelse(taxa %in% taxa_vec, taxa, "Other phyla")) %>% 
  group_by(lv1, taxa) %>% summarise_all(sum) %>% ungroup() %>% data.frame() %>% 
  filter(!lv1 %in% c("A09150","A09160","A09190")) %>% 
  profile_smp2grp(., group, unsmp_colnames = c("lv1", "taxa")) %>% 
  gather(key = "group", value = "val", -lv1, -taxa) %>% 
  mutate(taxa = factor(taxa, c(taxa_vec, "Other phyla")),
         group = factor(group, group_order),
         lv1 = ifelse(lv1 != "Unknown", db$lv1name[match(lv1, db$lv1)], "Unknown"))

taxa_color <- structure(c("#7dabc9", "#b37cb3", "#add26a", "#f4af63", "#89cbc0", "#ee7c6f", "#c8e3c1", "#f5e672", "#c86e8e", "#b9afab"), 
                        names = c(taxa_vec, "Other phyla"))

ggplot(plotdat, aes(group, val * 100, fill = taxa)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .1, width = .8) +
  facet_wrap(~lv1, scales = "free", nrow = 1) +
  scale_fill_manual(values = taxa_color) +
  scale_y_continuous(expand = c(.01, .01)) +
  labs(x = "", y = "Relative Abundance(%)", fill = "Phylum taxa") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 10, color = "black", face = "italic"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 10, color = "black"),
        panel.grid = element_line(linewidth = .4, color = "grey85"))
ggsave("kO_phylum_tpm_lv1_barplot.pdf", width = 18, height = 4)

#### kegg level A family ####
otu <- read.delim("kO_family_tpm", row.names = NULL)
db <- read.delim("KO_level_A_B_C_D_Description", header = F, quote = "") %>% 
  select(lv1 = V1, lv1name = V2, lv4 = V7) %>% filter(lv4 %in% otu$KO) %>% unique.data.frame()
db2 <- data.frame(unique(db[1:2]), row.names = NULL) %>% filter(!lv1 %in% c("A09150","A09160","A09190"))

# KO-->lv1 一对一
dat_x <- filter(otu, KO %in% get_freq(db$lv4, eq = 1)) %>% 
  mutate(lv1 = db$lv1[match(KO, db$lv4)]) %>% relocate(lv1)

# KO-->lv1 一对多
dat_y <- filter(otu, KO %in% get_freq(db$lv4, ne = 1))
tmp_list <- list()
pb <- txtProgressBar(style = 3)
x = 1
for (i in unique(dat_y$KO)) {
  kegg_vec <- db[db$lv4 %in% i, 1]
  tmp_dat_i <- filter(dat_y, KO %in% i)
  tmp_dat_i <- cbind(tmp_dat_i[,1:2], as.matrix(tmp_dat_i[,-(1:2)])/length(kegg_vec))
  tmp_dat_i <- map_dfr(kegg_vec, \(x) add_column(tmp_dat_i, lv1 = x, .before = "KO"))
  tmp_list[[i]] <- tmp_dat_i
  setTxtProgressBar(pb, x/length(unique(dat_y$KO))) 
  x = x + 1
}
close(pb)
dat_y = map_dfr(tmp_list, \(x) x)
dat_z = filter(otu, KO == "Unknown") %>% add_column(lv1 = "Unknown", .before = "KO")

dat <- reduce(list(dat_x, dat_y, dat_z), \(x, y) rbind(x, y))
dat <- cbind(dat[1:3], apply(dat[-(1:3)], 2, \(x) x/sum(x))) %>% data.frame() %>% select(-KO)

# select(dat, -lv1) %>% group_by(taxa) %>% summarise_all(sum) %>% column_to_rownames(var = "taxa") %>% 
#   rowSums() %>% data.frame(taxa = .) %>% arrange(desc(taxa)) %>% head(n = 30)
taxa_vec  <- c("Lactobacillaceae","Unknown","Enterococcaceae","Planococcaceae", "Lachnospiraceae","Bacteroidaceae",
               "Oscillospiraceae","Enterobacteriaceae", "Bacillaceae","Clostridiaceae","Bifidobacteriaceae")

plotdat <- mutate(dat, taxa = ifelse(taxa %in% taxa_vec, taxa, "Other families")) %>% 
  group_by(lv1, taxa) %>% summarise_all(sum) %>% ungroup() %>% data.frame() %>% 
  filter(!lv1 %in% c("A09150","A09160","A09190")) %>%  
  profile_smp2grp(., group, unsmp_colnames = c("lv1", "taxa")) %>% 
  gather(key = "group", value = "val", -lv1, -taxa) %>% 
  mutate(taxa = factor(taxa, c(taxa_vec, "Other families")),
         group = factor(group, group_order),
         lv1 = ifelse(lv1 != "Unknown", db$lv1name[match(lv1, db$lv1)], "Unknown"))

taxa_color <- structure(as.character(paletteer::paletteer_d("ggthemes::Tableau_20"))[1:12], names = c(taxa_vec, "Other families"))

ggplot(plotdat, aes(group, val * 100, fill = taxa)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .1, width = .8) +
  facet_wrap(~lv1, scales = "free", nrow = 1) +
  scale_fill_manual(values = taxa_color) +
  scale_y_continuous(expand = c(.01, .01)) +
  labs(x = "", y = "Relative Abundance(%)", fill = "Family taxa") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 10, color = "black", face = "italic"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 10, color = "black"),
        panel.grid = element_line(linewidth = .4, color = "grey85"))
ggsave("kO_family_tpm_lv1_barplot.pdf", width = 18, height = 4)

#### kegg level B metabolism ####
db <- read.delim("KO_level_A_B_C_D_Description", header = F) %>% 
  select(lv1 = V1, lv1ame = V2, lv2 = V3, lv2name = V4) %>% unique.data.frame() %>% 
  filter(lv1 == "A09100")

otu <- read.delim("KO_tpm.levelB", row.names = 1) %>% 
  apply(2, \(x) x/sum(x)) %>% data.frame()

plotdat <- dplyr::filter(otu, rownames(otu) %in% db$lv2) %>% t() %>% data.frame() %>% 
  mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
  gather(key = "lv2", value = "val", -group) %>% 
  mutate(lv2 = db$lv2name[match(lv2, db$lv2)],
         group = factor(group, group_order))

ggplot(plotdat, aes(group, val * 100)) +
  facet_wrap(~ factor(lv2), scale = "free") +
  geom_boxplot(position = position_dodge(), fill = "transparent", width = .6, 
               color = "#000000", linewidth = .3, outlier.shape = NA) +
  geom_jitter(aes(color = group), size = 1.5, width = .3, show.legend = F) +
  stat_smooth(aes(group = 1), level = 0.95, method = "lm",
              color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0, label.y.npc = 1, color = "#8da0cb") +
  stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5")), aes(group = 1), 
              level = 0.95, method = "lm", color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0.5, label.y.npc = 1, color = "#fc8d62") +
  stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5")), aes(group = 1), 
              level = 0.95, method = "lm", color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 1, label.y.npc = 1, color = "#66c2a5") +
  scale_color_manual(values = group_color) +
  labs(x = "", y = "") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 10, color = "black", face = "italic"),
        panel.background = element_rect(linewidth = .4, color = "black", fill = "transparent"),
        panel.grid = element_line(linewidth = .5, color = "grey85"))
ggsave("KO_tpm_not_to_RA_levelB_metabolism_boxplot.pdf", width = 11, height = 9)

# diff test 
source("F:/Code/R_func/diff_test.R")
diff <- diff_test_profile(otu, group, group_order)
write.table(diff, "KO_tpm_levelB_metabolism_diff_test.txt", sep = "\t", quote = F, row.names = F)

#### kegg level B overall ####
db <- read.delim("KO_level_A_B_C_D_Description", header = F) %>% 
  select(lv1 = V1, lv1ame = V2, lv2 = V3, lv2name = V4) %>% unique.data.frame()

otu <- read.delim("KO_tpm.levelB", row.names = 1) %>% 
  apply(2, \(x) x/sum(x)) %>% data.frame()

plotdat <- dplyr::filter(otu, rownames(otu) %in% db$lv2) %>% t() %>% data.frame() %>% 
  mutate(group = group$group[match(rownames(.), group$sample)]) %>% 
  gather(key = "lv2", value = "val", -group) %>% 
  mutate(lv2 = db$lv2name[match(lv2, db$lv2)],
         group = factor(group, group_order))

ggplot(plotdat, aes(group, val * 100)) +
  facet_wrap(~ factor(lv2), scale = "free") +
  geom_boxplot(position = position_dodge(), fill = "transparent", width = .6, 
               color = "#000000", linewidth = .3, outlier.shape = NA) +
  geom_jitter(aes(color = group), size = 1.5, width = .3, show.legend = F) +
  stat_smooth(aes(group = 1), level = 0.95, method = "lm",
              color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0, label.y.npc = 1, color = "#8da0cb") +
  stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5")), aes(group = 1), 
              level = 0.95, method = "lm", color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 0.5, label.y.npc = 1, color = "#fc8d62") +
  stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5")), aes(group = 1), 
              level = 0.95, method = "lm", color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
  stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5")), 
               aes(group = 1, label = paste(after_stat(eq.label), after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
               size = 1, label.x.npc = 1, label.y.npc = 1, color = "#66c2a5") +
  scale_color_manual(values = group_color) +
  labs(x = "", y = "TPM") +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 10, color = "black", face = "italic"),
        panel.background = element_rect(linewidth = .4, color = "black", fill = "transparent"),
        panel.grid = element_line(linewidth = .5, color = "grey85"))
ggsave("KO_tpm_not_to_RA_levelB_overall_boxplot.pdf", width = 24, height = 21)

# diff test 
source("F:/Code/R_func/diff_test.R")
diff <- diff_test_profile(otu, group, group_order)
write.table(diff, "KO_tpm_levelB_overall_diff_test.txt", sep = "\t", quote = F, row.names = F)

#### kegg level B phylum ####
source("F:/Code/R_func/profile_process.R")
source("F:/Code/R_func/utilities.R")
otu <- read.delim("kO_phylum_tpm", row.names = NULL)
db <- read.delim("KO_level_A_B_C_D_Description", header = F, quote = "") %>% 
  select(lv1 = V1, lv1name = V2, lv2 = V3, lv2name = V4, lv4 = V7) %>% 
  filter(lv4 %in% otu$KO) %>% unique.data.frame() %>% filter(lv1 == "A09100")
db2 <- data.frame(unique(db[1:4]), row.names = NULL)

# KO 一对一 lv1
dat_x <- filter(otu, KO %in% get_freq(db$lv4, eq = 1)) %>% 
  mutate(lv2 = db$lv2[match(KO, db$lv4)]) %>% relocate(lv2)
# KO 一对多 lv1
dat_y <- filter(otu, KO %in% get_freq(db$lv4, ne = 1))
tmp_list <- list()
pb <- txtProgressBar(style = 3)
x = 1
for (i in unique(dat_y$KO)) {
  kegg_vec <- db[db$lv4 %in% i, 3]
  tmp_dat_i <- filter(dat_y, KO %in% i)
  tmp_dat_i <- cbind(tmp_dat_i[,1:2], as.matrix(tmp_dat_i[,-(1:2)])/length(kegg_vec))
  tmp_dat_i <- map_dfr(kegg_vec, \(x) add_column(tmp_dat_i, lv2 = x, .before = "KO"))
  tmp_list[[i]] <- tmp_dat_i
  setTxtProgressBar(pb, x/length(unique(dat_y$KO))) 
  x = x + 1
}
close(pb)
dat_y = map_dfr(tmp_list, \(x) x)
dat_z = filter(otu, KO == "Unknown") %>% add_column(lv2 = "Unknown", .before = "KO")

dat <- reduce(list(dat_x, dat_y, dat_z), \(x, y) rbind(x, y))
dat <- cbind(dat[1:3], apply(dat[-(1:3)], 2, \(x) x/sum(x))) %>% data.frame() %>% select(-KO)

# select(dat, -lv1) %>% group_by(taxa) %>% summarise_all(sum) %>% column_to_rownames(var = "taxa") %>%
#   rowMeans() %>% data.frame(taxa = .) %>% arrange(desc(taxa))
taxa_vec  <- c("Bacillota", "Unknown", "Bacteroidota", "Pseudomonadota", "Actinomycetota", 
               "Uroviricota", "Campylobacterota", "Fusobacteriota", "Verrucomicrobiota")

plotdat <- mutate(dat, taxa = ifelse(taxa %in% taxa_vec, taxa, "Other phyla")) %>% 
  group_by(lv2, taxa) %>% summarise_all(sum) %>% ungroup() %>% data.frame() %>% 
  profile_smp2grp(., group, unsmp_colnames = c("lv2", "taxa")) %>% 
  gather(key = "group", value = "val", -lv2, -taxa) %>% 
  mutate(taxa = factor(taxa, c(taxa_vec, "Other phyla")),
         group = factor(group, group_order),
         lv2 = ifelse(lv2 != "Unknown", db$lv2name[match(lv2, db$lv2)], "Unknown"))
taxa_color <- structure(c("#7dabc9", "#b37cb3", "#add26a", "#f4af63", "#89cbc0", "#ee7c6f", "#c8e3c1", "#f5e672", "#c86e8e", "#b9afab"), 
                        names = c(taxa_vec, "Other phyla"))

ggplot(plotdat, aes(group, val * 100, fill = taxa)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .1, width = .8) +
  facet_wrap(~lv2, scales = "free", nrow = 3) +
  scale_fill_manual(values = taxa_color) +
  scale_y_continuous(expand = c(.01, .01)) +
  labs(x = "", y = "Relative Abundance(%)", fill = "Phylum taxa") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 10, color = "black", face = "italic"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 10, color = "black"),
        panel.grid = element_line(linewidth = .4, color = "grey85"))
ggsave("KO_phylum_tpm_lv2_barplot.pdf", width = 11, height = 9)


#### kegg level B family ####
source("F:/Code/R_func/profile_process.R")
source("F:/Code/R_func/utilities.R")
otu <- read.delim("KO_family_tpm", row.names = NULL)
db <- read.delim("KO_level_A_B_C_D_Description", header = F, quote = "") %>% 
  select(lv1 = V1, lv1name = V2, lv2 = V3, lv2name = V4, lv4 = V7) %>% 
  filter(lv4 %in% otu$KO) %>% unique.data.frame() %>% filter(lv1 == "A09100")
db2 <- data.frame(unique(db[1:4]), row.names = NULL)

# KO 一对一 lv1
dat_x <- filter(otu, KO %in% get_freq(db$lv4, eq = 1)) %>% 
  mutate(lv2 = db$lv2[match(KO, db$lv4)]) %>% relocate(lv2)
# KO 一对多 lv1
dat_y <- filter(otu, KO %in% get_freq(db$lv4, ne = 1))
tmp_list <- list()
pb <- txtProgressBar(style = 3)
x = 1
for (i in unique(dat_y$KO)) {
  kegg_vec <- db[db$lv4 %in% i, 3]
  tmp_dat_i <- filter(dat_y, KO %in% i)
  tmp_dat_i <- cbind(tmp_dat_i[,1:2], as.matrix(tmp_dat_i[,-(1:2)])/length(kegg_vec))
  tmp_dat_i <- map_dfr(kegg_vec, \(x) add_column(tmp_dat_i, lv2 = x, .before = "KO"))
  tmp_list[[i]] <- tmp_dat_i
  setTxtProgressBar(pb, x/length(unique(dat_y$KO))) 
  x = x + 1
}
close(pb)
dat_y = map_dfr(tmp_list, \(x) x)
dat_z = filter(otu, KO == "Unknown") %>% add_column(lv2 = "Unknown", .before = "KO")

dat <- reduce(list(dat_x, dat_y, dat_z), \(x, y) rbind(x, y))
dat <- cbind(dat[1:3], apply(dat[-(1:3)], 2, \(x) x/sum(x))) %>% data.frame() %>% select(-KO)

# select(dat, -lv1) %>% group_by(taxa) %>% summarise_all(sum) %>% column_to_rownames(var = "taxa") %>% 
#   rowSums() %>% data.frame(taxa = .) %>% arrange(desc(taxa)) %>% head(n = 30)
taxa_vec  <- c("Lactobacillaceae","Unknown","Enterococcaceae","Planococcaceae", "Lachnospiraceae","Bacteroidaceae",
               "Oscillospiraceae","Enterobacteriaceae", "Bacillaceae","Clostridiaceae","Bifidobacteriaceae")

plotdat <- mutate(dat, taxa = ifelse(taxa %in% taxa_vec, taxa, "Other families")) %>% 
  group_by(lv2, taxa) %>% summarise_all(sum) %>% ungroup() %>% data.frame() %>% 
  profile_smp2grp(., group, unsmp_colnames = c("lv2", "taxa")) %>% 
  gather(key = "group", value = "val", -lv2, -taxa) %>% 
  mutate(taxa = factor(taxa, c(taxa_vec, "Other families")),
         group = factor(group, group_order),
         lv2 = ifelse(lv2 != "Unknown", db$lv2name[match(lv2, db$lv2)], "Unknown"))
taxa_color <- structure(as.character(paletteer::paletteer_d("ggthemes::Tableau_20"))[1:12], names = c(taxa_vec, "Other families"))

ggplot(plotdat, aes(group, val * 100, fill = taxa)) +
  geom_bar(stat = "identity", position = position_stack(), color = "black", lwd = .1, width = .8) +
  facet_wrap(~lv2, scales = "free", nrow = 3) +
  scale_fill_manual(values = taxa_color) +
  scale_y_continuous(expand = c(.01, .01)) +
  labs(x = "", y = "Relative Abundance(%)", fill = "Family taxa") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        strip.text = element_text(size = 10, color = "black", face = "italic"),
        legend.text = element_text(size = 8, color = "black", face = "italic"),
        legend.title = element_text(size = 10, color = "black"),
        panel.grid = element_line(linewidth = .4, color = "grey85"))
ggsave("KO_family_tpm_lv2_barplot.pdf", width = 11, height = 9)

#### kegg level C metabolism ####
db <- read.delim("KO_level_A_B_C_D_Description", header = F) %>% 
  select(lv1 = V1, lv1ame = V2, lv2 = V3, lv2name = V4, lv3 = V5, lv3name = V6) %>% unique.data.frame() %>% 
  filter(lv1 == "A09100")

otu <- read.delim("KO_tpm.levelC", row.names = 1) %>% 
  apply(2, \(x) x/sum(x)) %>% data.frame()

# tmp_sample <- filter(group, group %in% c("d1", "d3", "d5")) %>% select(sample) %>% unlist() %>% as.character()
tmp_sample <- filter(group, group %in% group_order) %>% select(sample) %>% unlist() %>% as.character()
plotdat <- dplyr::filter(otu, rownames(otu) %in% db$lv3) %>% select(all_of(tmp_sample)) %>% t() %>% data.frame()
lm_res <- data.frame(id = character(), coef = numeric(), slope = numeric(), rr = numeric(), pval = numeric())
for (i in 1:ncol(plotdat)) {
  tmp_var <- colnames(plotdat)[i]
  dat <- select(plotdat, val = all_of(i)) %>% rownames_to_column(var = "sample") %>% merge(group, by = "sample") %>% 
    mutate(group = factor(group, group_order, labels = seq_len(9)) %>% as.numeric(), val = val * 100)
  lm <- lm(val ~ group, dat)
  lm_sum <- summary(lm)
  lm_res <- add_row(lm_res, 
                    id = tmp_var,
                    coef = lm$coefficients[1], 
                    slope = lm$coefficients[2], 
                    rr = lm_sum$r.squared, 
                    pval = pf(lm_sum$fstatistic[1],lm_sum$fstatistic[2],lm_sum$fstatistic[3], lower.tail = F))
}
lm_res <- mutate(lm_res, name = db$lv3name[match(id, db$lv3)], .after = id)
write.table(lm_res, "KO_tpm_levelC_lm_d1-42.txt", sep = "\t", row.names = F, quote = F)

# plotdat <- dplyr::filter(otu, rownames(otu) %in% db$lv3) %>% t() %>% data.frame() %>% 
#   mutate(group = group$group[match(rownames(.), group$sample)]) %>%
#   gather(key = "lv3", value = "val", -group) %>%
#   mutate(lv3 = db$lv3name[match(lv3, db$lv3)],
#          group = factor(group, group_order))

# ggplot(plotdat, aes(group, val * 100)) +
#   facet_wrap(~ factor(lv3), scale = "free") +
#   geom_boxplot(position = position_dodge(), fill = "transparent", width = .6, 
#                color = "#000000", linewidth = .3, outlier.shape = NA) +
#   geom_jitter(aes(color = group), size = 1.5, width = .3, show.legend = F) +
#   stat_smooth(aes(group = 1), level = 0.95, method = "lm",
#               color = "#8da0cb", lwd = .5, fill = "grey85", lty = "solid") +
#   stat_poly_eq(aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
#                size = 1, label.x.npc = 0, label.y.npc = 1, color = "#8da0cb") +
#   stat_smooth(data = filter(plotdat, group %in% c("d1","d3","d5")), aes(group = 1), 
#               level = 0.95, method = "lm", color = "#fc8d62", lwd = .5, fill = "grey85", lty = "solid") +
#   stat_poly_eq(data = filter(plotdat, group %in% c("d1","d3","d5")), 
#                aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
#                size = 1, label.x.npc = 0.5, label.y.npc = 1, color = "#fc8d62") +
#   stat_smooth(data = filter(plotdat, !group %in% c("d1","d3","d5")), aes(group = 1), 
#               level = 0.95, method = "lm", color = "#66c2a5", lwd = .5, fill = "grey85", lty = "solid") +
#   stat_poly_eq(data = filter(plotdat, !group %in% c("d1","d3","d5")), 
#                aes(group = 1, label = paste(after_stat(rr.label), after_stat(p.value.label), sep='~`,`~')),
#                size = 1, label.x.npc = 1, label.y.npc = 1, color = "#66c2a5") +
#   scale_color_manual(values = group_color) +
#   labs(x = "", y = "Relative Abundance (%)") +
#   theme_bw() +
#   theme(axis.line = element_blank(),
#         axis.ticks = element_line(linewidth = .4, color = "black"),
#         axis.text = element_text(size = 8, color = "black"),
#         axis.title = element_text(size = 8, color = "black"),
#         strip.text = element_text(size = 10, color = "black", face = "italic"),
#         panel.background = element_rect(linewidth = .4, color = "black", fill = "transparent"),
#         panel.grid = element_line(linewidth = .5, color = "grey85"))
# ggsave("KO_tpm_levelC_metabolism_boxplot.pdf", width = 40, height = 40, limitsize = F)

# diff test 
source("F:/Code/R_func/diff_test.R")
diff <- diff_test_profile(otu, group, group_order)
write.table(diff, "KO_tpm_levelB_diff_test.txt", sep = "\t", quote = F, row.names = F)
