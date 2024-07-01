##############################################################################
#
# eFC and agitation
# 
# For this pipeline, we explored the DMN-SMN interaction with respect
# to the agitation severity.
#
# Functional connectivity was obtained from the ETS estimation and divided into
# Peak FC and Trough FC
#
# Liang Qunjun 2023/11/13

library(tidyverse)
library(bruceR)
library(ggstatsplot)
library(ggridges)
library(psych)
library(RColorBrewer)
library(emmeans)
library(ggeasy)
library(ggsci)
library(cowplot)
library(scales)
library(ggsignif)
library(patchwork)
library(sjPlot)
source("scripts/function_ShowCorMat.R")
source("scripts/function_ExtractConnectivity.R")

# load the data
dat_fc <- rio::import("inputs/Analysis3_FC_extraction.xlsx")
dat_hama <- dat_fc %>% select(starts_with("HAMA"))
dat_fc <- dat_fc %>% mutate(HAMA_wave1_total = rowSums(dat_hama))

#################################################
#
# statistics - group difference
#
##################################################
dat_fc %>% 
  filter(FC == 'FC') %>% # static FC
  filter(subtype_psycho != "NA-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "SMN_DMN", 
                 between = "subtype_psycho", 
                 covariate = c("gender","age","education")) %>% 
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# Pairwise Comparisons of "group":
# ───────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────
# MDD - HC   -0.028 (0.014) 133 -2.102  .037 *   -0.425 [-0.824, -0.025]
# ───────────────────────────────────────────────────────────────────────
dat_fc %>% 
  filter(FC == 'peak') %>% # high-amplitude FC
  filter(subtype_psycho != "NA-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "SMN_DMN", 
                 between = "subtype_psycho", 
                 covariate = c("gender","age","education")) %>% 
  bruceR::EMMEANS(effect = "subtype_psycho")
# ───────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────
# MDD - HC   -0.048 (0.029) 133 -1.671  .097 .    -0.338 [-0.737, 0.062]
# ───────────────────────────────────────────────────────────────────────
dat_fc %>% 
  filter(FC == 'trough') %>% # low-amplitude FC
  filter(subtype_psycho != "NA-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "SMN_DMN", 
                 between = "group", 
                 covariate = c("gender","age","education")) %>% 
  bruceR::EMMEANS(effect = "group")
# ───────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────
# MDD - HC   -0.019 (0.009) 133 -2.139  .034 *   -0.432 [-0.832, -0.032]
# ───────────────────────────────────────────────────────────────────────

p.adjust(c(0.037, 0.034, 0.097), method = "fdr")
# 0.0555 0.0555 0.0970

# NA-MDD vs. HC -----------------------------------------------------
dat_fc %>% 
  filter(FC == 'FC') %>%
  filter(subtype_psycho != "A-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "SMN_DMN", 
                 between = "subtype_psycho", 
                 covariate = c("gender","age","education")) %>% 
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# Pairwise Comparisons of "subtype_psycho":
# ────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - HC   -0.017 (0.012) 163 -1.454  .148      -0.259 [-0.610, 0.093]
# ────────────────────────────────────────────────────────────────────────────

dat_fc %>% 
  filter(FC == 'peak') %>%
  filter(subtype_psycho != "A-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "SMN_DMN", 
                 between = "group", 
                 covariate = c("gender","age","education")) %>% 
  bruceR::EMMEANS(effect = "group")
# Pairwise Comparisons of "group":
# ───────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────
# MDD - HC   -0.026 (0.026) 163 -0.985  .326      -0.175 [-0.527, 0.176]
# ───────────────────────────────────────────────────────────────────────

dat_fc %>% 
  filter(FC == 'trough') %>%
  filter(subtype_psycho != "A-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "SMN_DMN", 
                 between = "group", 
                 covariate = c("gender","age","education")) %>% 
  bruceR::EMMEANS(effect = "group")
# Pairwise Comparisons of "group":
# ───────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────
# MDD - HC   -0.013 (0.008) 163 -1.620  .107      -0.288 [-0.639, 0.063]
# ───────────────────────────────────────────────────────────────────────

p.adjust(c(0.148, 0.107, .326), method = "fdr")
#  0.222 0.222 0.326

#############################################################################
#
# Obtain the node and edge parameter. Visualization should be proceed in MATLAB
# using BrainNetViewer
#
############################################################################

## reshape the network annotation 
net_anna <- rio::import("inputs/Power264_Yeo7.xlsx")
net_anna["ROI.Name"] <- paste0("R",formatC(1:214, digits = 2, flag = "0"))
net_anna <- net_anna %>%
  rename(x.mni = "x", y.mni="y", z.mni="z")
net_anna$x.mni <- as.integer(net_anna$x.mni)
net_anna$y.mni <- as.integer(net_anna$y.mni)
net_anna$z.mni <- as.integer(net_anna$z.mni)


x <- example_weighted_undirected
brainconn(atlas ="schaefer300_n7", conmat=x, node.size = 3, view="ortho", edge.color.weighted = T)

net_anna_plot <- net_anna %>% mutate(net_color = ifelse(network_new == "somMot", "cyan",
                                                   ifelse(network_new == "ATN","green",
                                                          ifelse(network_new == "DMN", "red",
                                                                 ifelse(network_new == "salience", "purple",
                                                                        ifelse(network_new == "visual", "blue",
                                                                               ifelse(network_new == "FPN", "orange", "black"))))))) %>%
  rename(network = "network_new") %>%
  filter(network == "somMot" | network == "DMN")
check_atlas(net_anna_plot)

net_anna_plot[,c(1,2,3,8)] %>%
  mutate(cluster = ifelse(network == "somMot", 1, 2)) %>%
  mutate(size = 1) %>%
  select(1:3,5,6,4)  %>%
readr::write_delim(file = "outputs/FC_network_plto.node", col_names = F, delim = "\t")

## obtain A-MDD FC 
sbj_list <- dat_use %>% filter(subtype_psycho == "A-MDD") %>% .$participant_id
file_top <- file_topFC_names[grep(paste(sbj_list, collapse = "|"), file_topFC_names)]
file_bottom <- file_bottomFC_names[grep(paste(sbj_list, collapse = "|"), file_bottomFC_names)]
file_FC <- file_FC_names[grep(paste(sbj_list, collapse = "|"), file_FC_names)]

fc_amdd <- map(file_FC, ~ ExtractSMNDMN(.x, net_annotation = net_anna)) %>%
  reduce(`+`) 
fc_amdd_plot <- fc_amdd/length(sbj_list)
diag(fc_amdd_plot) <- 0  
readr::write_delim(data.frame(fc_amdd_plot), file = "outputs/FC_network_plto_AMDD.edge", col_names = F, delim = "\t")

## obtain NA-MDD FC 
sbj_list <- dat_use %>% filter(subtype_psycho == "NA-MDD") %>% .$participant_id
file_top <- file_topFC_names[grep(paste(sbj_list, collapse = "|"), file_topFC_names)]
file_bottom <- file_bottomFC_names[grep(paste(sbj_list, collapse = "|"), file_bottomFC_names)]
file_FC <- file_FC_names[grep(paste(sbj_list, collapse = "|"), file_FC_names)]

fc_amdd <- map(file_FC, ~ ExtractSMNDMN(.x, net_annotation = net_anna)) %>%
  reduce(`+`) 
fc_amdd_plot <- fc_amdd/length(sbj_list)
diag(fc_amdd_plot) <- 0  
readr::write_delim(data.frame(fc_amdd_plot), file = "outputs/FC_network_plto_NAMDD.edge", col_names = F, delim = "\t")

## obtain HC FC 
sbj_list <- dat_use %>% filter(subtype_psycho == "HC") %>% .$participant_id
file_top <- file_topFC_names[grep(paste(sbj_list, collapse = "|"), file_topFC_names)]
file_bottom <- file_bottomFC_names[grep(paste(sbj_list, collapse = "|"), file_bottomFC_names)]
file_FC <- file_FC_names[grep(paste(sbj_list, collapse = "|"), file_FC_names)]

fc_amdd <- map(file_FC, ~ ExtractSMNDMN(.x, net_annotation = net_anna)) %>%
  reduce(`+`) 
fc_amdd_plot <- fc_amdd/length(sbj_list)
diag(fc_amdd_plot) <- 0  
readr::write_delim(data.frame(fc_amdd_plot), file = "outputs/FC_network_plto_HC.edge", col_names = F, delim = "\t")

#################################################
#
# Plot the correlation matrix by amplitude and trough FCs
#
##################################################
## group average
# indicate the path to the inputs and participant information
path_ets <- "inputs/timeseries_Power/Edge_time_series_estimation/"
file_topFC_names <- list.files(path_ets, pattern = 'sub-[0-9]*_topFC.csv', full.names = T)
file_bottomFC_names <- list.files(path_ets, pattern = 'sub-[0-9]*_bottomFC.csv', full.names = T)

# obtain FC in health control --------------------------------------------------
sbj_list_hc <- dat_use %>% filter(group == "HC") %>% .$participant_id
file_top <- file_topFC_names[grep(paste(sbj_list_hc, collapse = "|"), file_topFC_names)]
file_bottom <- file_bottomFC_names[grep(paste(sbj_list_hc, collapse = "|"), file_bottomFC_names)]

fc_hc_top <- matrix(0, ncol = 214, nrow = 214)
for (i in file_top) {
  mat_tmp <- readr::read_csv(i, col_names = F, show_col_types = F) %>% as.matrix()
  fc_hc_top <- (fc_hc_top + mat_tmp)/2
}

p_mat_hc_top <- ShowCorMat(result_mat = fc_hc_top, title_text = "Peak FC in HC")
p_mat_hc_top

fc_hc_bottom <- matrix(0, ncol = 214, nrow = 214)
for (i in file_bottom) {
  mat_tmp <- readr::read_csv(i, col_names = F, show_col_types = F) %>% as.matrix()
  fc_hc_bottom <- (fc_hc_bottom + mat_tmp)/2
}

# obtain agited MDD ------------------------------------------------------------
sbj_list_agit <- dat_use %>% filter(subtype_psycho == "A-MDD") %>% .$participant_id
file_top <- file_topFC_names[grep(paste(sbj_list_agit, collapse = "|"), file_topFC_names)]
file_bottom <- file_bottomFC_names[grep(paste(sbj_list_agit, collapse = "|"), file_bottomFC_names)]

fc_agit_top <- matrix(0, ncol = 214, nrow = 214)
for (i in file_top) {
  mat_tmp <- readr::read_csv(i, col_names = F, show_col_types = F) %>% as.matrix()
  fc_agit_top <- (fc_agit_top + mat_tmp)/2
}

p_mat_agit_top <- ShowCorMat(result_mat = fc_agit_top, title_text = "Peak FC in Agited MDD")
p_mat_agit_top

fc_agit_bottom <- matrix(0, ncol = 214, nrow = 214)
for (i in file_bottom) {
  mat_tmp <- readr::read_csv(i, col_names = F, show_col_types = F) %>% as.matrix()
  fc_agit_bottom <- (fc_agit_bottom + mat_tmp)/2
}

# obtain non-agited MDD ------------------------------------------------------------
sbj_list_Nagit <- dat_use %>% filter(subtype_psycho == "NA-MDD") %>% .$participant_id
file_top <- file_topFC_names[grep(paste(sbj_list_Nagit, collapse = "|"), file_topFC_names)]
file_bottom <- file_bottomFC_names[grep(paste(sbj_list_Nagit, collapse = "|"), file_bottomFC_names)]

fc_Nagit_top <- matrix(0, ncol = 214, nrow = 214)
for (i in file_top) {
  mat_tmp <- readr::read_csv(i, col_names = F, show_col_types = F) %>% as.matrix()
  fc_Nagit_top <- (fc_Nagit_top + mat_tmp)/2
}

p_mat_Nagit_top <- ShowCorMat(result_mat = fc_Nagit_top, title_text = "Peak FC in non-agited MDD")
p_mat_Nagit_top

fc_Nagit_bottom <- matrix(0, ncol = 214, nrow = 214)
for (i in file_bottom) {
  mat_tmp <- readr::read_csv(i, col_names = F, show_col_types = F) %>% as.matrix()
  fc_Nagit_bottom  <- (fc_Nagit_bottom  + mat_tmp)/2
}

# combine the figures ------------------------------------------------------------

# combine the matrix for HC
merged_matrix <- matrix(NA, nrow = nrow(fc_hc_bottom), ncol = ncol(fc_hc_bottom))
merged_matrix[lower.tri(fc_hc_top)] <- fc_hc_top[lower.tri(fc_hc_top)]
merged_matrix[upper.tri(fc_hc_bottom)] <- fc_hc_bottom[upper.tri(fc_hc_bottom)]
p1 <- ShowCorMat_whole(merged_matrix, title_text = "Health control")
p1

# combine the matrix for A-MDD
merged_matrix <- matrix(NA, nrow = nrow(fc_agit_bottom), ncol = ncol(fc_agit_bottom))
merged_matrix[lower.tri(fc_agit_top)] <- fc_agit_top[lower.tri(fc_agit_top)]
merged_matrix[upper.tri(fc_agit_bottom)] <- fc_agit_bottom[upper.tri(fc_agit_bottom)]
p2 <- ShowCorMat_whole(merged_matrix, title_text = "Agitated MDD")
p2

# combine the matrix for NA-MDD 
merged_matrix <- matrix(NA, nrow = nrow(fc_Nagit_bottom), ncol = ncol(fc_Nagit_bottom))
merged_matrix[lower.tri(fc_agit_top)] <- fc_Nagit_top[lower.tri(fc_Nagit_top)]
merged_matrix[upper.tri(fc_agit_bottom)] <- fc_Nagit_bottom[upper.tri(fc_Nagit_bottom)]
p3 <- ShowCorMat_whole(merged_matrix, title_text = "Non-agitated MDD")
p3 + easy_remove_x_axis(what = "text")

p_fc_merge <- p1 + p2 + p3 # merge FC matrices 
p_fc_merge
ggsave(plot = p_fc_merge, filename = "outputs/Analysis4_eFC.png", dpi = 300, width = 20, height = 6)

##################################################################
#
# Plot the group difference by boxplot
#
###################################################################

f1 <- dat_fc %>% filter(subtype_psycho != "NA-MDD") %>% filter(FC == "peak") %>% .$subtype_psycho
f1 <- c(.7, 1.7)[as.integer(factor(f1))]

p_diff_fc <- dat_fc %>% filter(subtype_psycho != "NA-MDD") %>%
  filter(FC == "FC") %>%
  ggplot(aes(x = subtype_psycho, y = SMN_DMN, fill = subtype_psycho)) +
  geom_boxplot(size = 1, width = .4, outlier.colour = NA, alpha = .6) +
  geom_point(aes(x = f1), size = 4, alpha = .6, position = position_jitter(.03)) +
  geom_signif(annotations = rep("*", times = 1),
              textsize = 7, vjust = .7,
              y_position = c(0.28),
              xmin = c(1),
              xmax = c(2),
              tip_length = 0) +
  ylim(c(-.15,.3)) + ylab("Connectivity (a.u.)") +
  ggtitle(label = "static FC") +
  scale_fill_lancet() +
  theme_classic() +
  easy_remove_legend() +
  easy_text_size(13) +
  easy_center_title() +
  easy_remove_x_axis(what = "title")
p_diff_fc

p_diff_peak <- dat_fc %>% filter(subtype_psycho != "NA-MDD") %>%
  filter(FC == "peak") %>%
  ggplot(aes(x = subtype_psycho, y = SMN_DMN, fill = subtype_psycho)) +
  geom_boxplot(size = 1, width = .4, outlier.colour = NA, alpha = .6) +
  geom_point(aes(x = f1), size = 4, alpha = .6, position = position_jitter(.03)) +
  ylim(c(-.25,.3)) + ylab("Connectivity (a.u.)") +
  ggtitle(label = "high-amplitude FC") +
  scale_fill_lancet() +
  theme_classic() +
  easy_remove_legend() +
  easy_text_size(13) +
  easy_center_title() +
  easy_remove_x_axis(what = "title")
p_diff_peak

p_diff_trough <- dat_fc %>% filter(subtype_psycho != "NA-MDD") %>%
  filter(FC == "trough") %>%
  ggplot(aes(x = subtype_psycho, y = SMN_DMN, fill = subtype_psycho)) +
  geom_boxplot(size = 1, width = .4, outlier.colour = NA, alpha = .6) +
  geom_point(aes(x = f1), size = 4, alpha = .6, position = position_jitter(.03)) +
  geom_signif(annotations = rep("*", times = 1),
              textsize = 7, vjust = .7,
              y_position = c(0.15),
              xmin = c(1),
              xmax = c(2),
              tip_length = 0) +
  ylim(c(-.12,.16)) + ylab("Connectivity (a.u.)") +
  ggtitle(label = "low-amplitude FC") +
  scale_fill_lancet() +
  theme_classic() +
  easy_remove_legend() +
  easy_text_size(13) +
  easy_center_title() +
  easy_remove_x_axis(what = "title")
p_diff_trough

p_diff_boxplot <- p_diff_fc + p_diff_peak + p_diff_trough
p_diff_boxplot

ggsave(plot = p_diff_boxplot, filename = "outputs/Fig3C.png", dpi = 300,
       width = 14, height = 6)
