# Jinxin Meng, 20231226, 20240608 ------------------------------------------------

pacman::p_load(tibble, dplyr, tidyr, ggplot2)
setwd("F:/Proj/proj_2023/01.Broiler_fecal_resistome_20230302/bacterial_genome/")

# ggvenn -------------------------------------------------------------------------

library(ggvenn)

dat <- read.delim("strains_fna_cls_95_adj")

# adj2list

plotdat <- dat %>%
  mutate(pattern_d = ifelse(pattern_d > 0, cluster, NA),
         pattern_HC = ifelse(pattern_HC > 0, cluster, NA)) %>% 
  column_to_rownames(var = "cluster") %>% 
  as.list() %>% 
  lapply(., function(x) as.character(na.omit(x)))

names(plotdat) <- c("Isolate \n(42 species)", 
                    "Assembly \n(326 species)")

ggvenn(plotdat, fill_color = c("#e18d73", "#91c6b6", "#60aac8"), fill_alpha = .6,
       stroke_color = "black", stroke_size = .6, stroke_linetype = "dashed",
       set_name_color = "black", set_name_size = 3, text_size = 3)

ggsave("vennplot.pdf", height = 3, width = 3)

# ggvenn taxa stat ---------------------------------------------------------------

dat %>% 
  filter(pattern_HC > 0 & pattern_d == 0) %>% 
  nrow()


#### genome_info_stat ####
# genome size and N50
dat <- read.delim("strains_fna_statistics.tsv")
plotdat <- dat %>% 
  select(total_len,  GC) %>% 
  mutate(GenomeSize = total_len / 1e6)

p <- ggplot(plotdat, aes(x = GenomeSize, y = GC)) +
  geom_point(color = "#63a5c8", size = 1) +
  geom_vline(xintercept = 2.038571, lty = 2, linewidth = .4) +
  geom_hline(yintercept = 45.30158, lty = 2, linewidth = .4) +
  scale_x_continuous(expand = c(.03, .03)) +
  scale_y_continuous(expand = c(.03, .03)) +
  labs(x = "Genome Size (Mbp)", y = "GC content(%)") +
  theme_bw() +
  theme(axis.line = element_line(linewidth = .4, color = "black"),
        axis.ticks = element_line(linewidth = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)
p <- ggExtra::ggMarginal(p, type = "histogram", fill = "#63a5c8", size = 6, lwd = .3)
ggsave(p, file = "genome_statistics_GC_len_scatterplot.pdf", width = 5, height = 5)

# completeness and contamination
plotdat <- read.delim("strains_fna_qs", header = F) %>% 
  select(feature = V1, completeness = V2, contamination = V3)

p <- ggplot(plotdat, aes(x = completeness, y = contamination)) +
  geom_point(color = "#f5b75d", size = 1) +
  geom_vline(xintercept = 87.58035, lty = 2, linewidth = .4) +
  geom_hline(yintercept = 1.330865, lty = 2, linewidth = .4) +
  scale_x_continuous(expand = c(.03, .03)) +
  scale_y_continuous(expand = c(.03, .03)) +
  labs(x = "Completeness (%)", y = "Contamination (%)") +
  theme_bw() +
  theme(axis.line = element_line(size = .4, color = "black"),
        axis.ticks = element_line(size = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 10, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)
p <- ggExtra::ggMarginal(p, type = "histogram", fill = "#f5b75d", size = 6, lwd = .3)
ggsave(p, file = "genome_statistics_conpleteness_and_contamination_scatterplot.pdf", width = 5, height = 5)
