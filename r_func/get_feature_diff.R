# Encoding: utf-8
# Author: Jinxin Meng
# Email: mengjx855@163.com
# Created Data：2022-9-6
# Modified Data: 2023-11-01
# Version: 1.3
# 2023-10-25: add xlim parameter in plot_roc.
# 2023-10-25: add trans parameter in get_feature_diff.
# 2023-10-25: add map_name parameter in the get_feature_diff.
# 2023-10-25: add names in the plot_volcano.
# 2023-11-01: add prevalence statistics in output of get_feature_diff, min_abundance as a threshold of presence.
# 2023-11-02: modify group_color in plot_volcano using function structure

library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

#### Function1 ####
# 对feature进行组间富集分析，绘制火山图并保存数据
# otu输入正常的otu表
# group表中分组列必须为group
# group_pair传入一个向量；例如c("A","B")，之后的分析即为AvsB，上调在A中富集，下调在B中富集，初步版本只能是"Disease", "Control"
# log2fc为认定富集的差异倍数阈值
# padj为认定富集的显著性阈值;
# group[data.frame]: metadata contain sample and group column, also specify by map_names parameter.
# trans[log]: used in difference analysis to transfrom profile, only including LOG in current version. [NULL, LOG]
# min_abundace[num]: the threshold of a feature presencing when evaluate prevalence rate (%).
get_feature_diff <- function(otu, group, group_pair = NULL, trans = NULL, min_abundance = 0,
                             group_colnames = NULL, diff_method = "wilcox", qval = 0.05, log2fc = 1) {
  if (!all(colnames(group) %in% c("sample", "group")) & is.null(group_colnames)) stop("group field (sample|group)")
  if (!is.null(group_colnames)) group <- data.frame(group, check.names = F) %>% dplyr::rename(all_of(group_colnames))

  # difference analysis order is define.
  if (is.null(group_pair)) { group_pair <- unique(group$group)
  } else { group_pair <- group_pair }
  
  # is the profile file transformation in DA?
  profile <- data.frame(otu, check.names = F)
  if(is.null(trans)) {
    dat_DA <- data.frame(t(profile), check.names = F) %>% mutate(group = group$group[match(rownames(.), group$sample)])
  } else if (trans == "LOG") {
    source("F:/Code/R_Func/utilities.R")
    dat_DA <- data.frame(t(trans_LOG(profile)), check.names = F) %>% 
      mutate(group = group$group[match(rownames(.), group$sample)])
  } else {
    stop("in profile transformation, parameter only including LOG in current version.")
  }
  
  # difference analysis using diff_method.
  feature <- setdiff(colnames(dat_DA), "group")
  cat(paste0("[", Sys.time(), "] Hypothesis testing using ", diff_method, " method.\n"))
  pb <- txtProgressBar(style = 3)
  pval_vec <- c()
  for (i in 1:length(feature)) {
    dat_i <- select(dat_DA, feature[i], group)
    if (diff_method == "wilcox") { # wilcox
      test <- wilcox.test(unlist(subset(dat_i, group%in%group_pair[1])[1]),
                          unlist(subset(dat_i, group%in%group_pair[2])[1]),
                          paired = F, exact = F)
    } else if (diff_method == "t") { # t
      test <- stats::t.test(unlist(subset(dat_i, group%in%group_pair[1])[1]),
                            unlist(subset(dat_i, group%in%group_pair[2])[1]), 
                            paired = F)
    }
    pval_vec[i] <- test$p.value
    setTxtProgressBar(pb, i/length(feature)) 
  }
  close(pb)
  diff <- data.frame(feature = feature, pval = pval_vec) %>% 
    mutate(padj = p.adjust(pval, method = "BH"))
  
  # calu feature abundance and foldchange
  cat(paste0("[", Sys.time(), "] Abundance analysis.\n"))
  dat <- data.frame(t(profile), check.names = F) %>% mutate(group = group$group[match(rownames(.), group$sample)])
  abund <- dat %>%
    group_by(group) %>%
    summarise_all(mean) %>%
    column_to_rownames(var = "group") %>% 
    t(.) %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column(var = "feature") %>% 
    dplyr::select(feature, index1_abund = all_of(group_pair[1]), index2_abund = all_of(group_pair[2])) %>% 
    mutate(comparsion = paste(group_pair, collapse = "_vs_"),
           FC = index1_abund / index2_abund,
           log2FC = log2(FC))
  
  # calu feature prevalence
  cat(paste0("[", Sys.time(), "] Prevalence statistic.\n"))
  prevalence <- dat %>% 
    group_by(group) %>% 
    group_modify(~ purrr::map_df(.x,\(x) sum(x > min_abundance)/length(x) *100)) %>% 
    column_to_rownames(var = "group") %>% t() %>% 
    data.frame(check.names = F) %>% 
    rownames_to_column(var = "feature") %>% 
    dplyr::select(feature, index1_prevalence = all_of(group_pair[1]), index2_prevalence = all_of(group_pair[2]))
  
  # reshape the output
  cat(paste0("[", Sys.time(), "] Reshape result.\n"))  
  res <- merge(abund, prevalence, by = "feature", all = T) %>% 
    merge(., diff, by = "feature", all = T) %>% 
    relocate(index1_prevalence, .after = index1_abund) %>% 
    relocate(index2_prevalence, .after = index2_abund)
  res <- res %>% subset(!is.nan(FC)) %>% 
    within(., {
           FC[FC == Inf] <- max(FC[!is.infinite(FC)])
           log2FC[log2FC == Inf] <- max(log2FC[!is.infinite(log2FC)])
           log2FC[log2FC == -Inf] <- min(log2FC[!is.infinite(log2FC)])
           }) %>% 
    mutate(enriched = ifelse(log2FC > log2fc & padj < qval, group_pair[1],
                             ifelse(log2FC < -log2fc & padj < qval, group_pair[2], "None")))
  
  # rename the output
  names_vec <- c("index1_abund", "index2_abund", "index1_prevalence", "index2_prevalence")
  names(names_vec) <- c(paste0(group_pair[1], "_abund"), paste0(group_pair[2], "_abund"), 
                        paste0(group_pair[1], "_prevalence"), paste0(group_pair[2], "_prevalence"))
  res <- res %>% rename(all_of(names_vec))
  cat(paste0("[", Sys.time(), "] Program end.\n"))
  return(res)
}


