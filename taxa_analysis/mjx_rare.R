#### info ####
# Jinxin Meng, 2023-10-09, 2023-11-12
pacman::p_load(tibble, dplyr, ggpubr)
setwd("/share/data1/mjx/proj/01.broiler_metagenome_20230301/analysis/taxa")
source("/share/data1/mjx/Rproj/script/rare_curve.R")

#### rarefaction ####
group <- read.delim("../data/sample_group")
otu <- read.delim("../data/species_tpm_filter", row.names = 1) %>% 
  apply(., 2, function(x) x/sum(x)) %>% 
  data.frame()
taxonomy <- read.delim("../data/species_taxonomy_rename.tsv")

dat <- rare_specaccum(otu)
write.table(dat, "rare_sample.txt", sep = "\t", row.names = F, quote = F)

plotdat <- dat %>%
  add_row(sample = 0, obs = 0, sd = 0) %>% 
  filter(sample%%3==1)

ggplot(plotdat, aes(x = sample, y = obs)) +
  geom_ribbon(aes(ymin = obs-sd, ymax = obs+sd), fill = "#fcbba1") +
  geom_line(lwd = .4, color = "#000000") +
  # geom_errorbar(aes(ymin = obs - sd, ymax = obs + sd), width = 1, lwd = .4, color = "red") +
  labs(x = "Number of samples", y = "Richness", color = NULL) +
  scale_color_manual(values = colors) +
  scale_x_continuous(expand = c(.02, .02)) +
  scale_y_continuous(expand = c(.02, .02)) +
  theme_classic() +
  theme(axis.line = element_line(size = .4, color = "black"),
        axis.ticks = element_line(size = .4, color = "black"),
        axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8, color = "black"),
        panel.grid = element_blank(),
        aspect.ratio = 1)
ggsave("rare_sample.pdf", width = 4, height = 4)
dev.off()
