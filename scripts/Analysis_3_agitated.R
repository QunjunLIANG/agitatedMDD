#############################################################################
#
# Temporal feature abnormality for Agitated MDD 
#
# For this pipeline, we tested the temporal structure among agitated MDD, HC 
# and non-agitated MDD. 
#
# The alternative features include:
#   1. Peak duration (ETS)
#   2. Number of Peaks (ETS)
#   3. Time delay in SMN (TD)
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
library(patchwork)
library(cowplot)
library(scales)
library(ggsignif)
library(patchwork)
library(sjPlot)
source("scripts/function_PvalueForTable1.R")
# load the data
dat_use <- rio::import("inputs/Analysis2_collect_fMRI.xlsx")

#################################################
#
# coupling between peak duration and peak height
#
#################################################

model_lm <- dat_use %>% filter(subtype_psycho != "NA-MDD") %>%
  lm(peak_mean ~ duration_mean*subtype_psycho + 
       age + gender + education, data = .)
bruceR::GLM_summary(model_lm) # check the results
# ───────────────────────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p        [95% CI of β] r(partial) r(part)
# ───────────────────────────────────────────────────────────────────────────────────────────────────────
# duration_mean                        0.648 (0.142)  4.545 <.001 *** [ 0.367,  0.929]      0.301   0.283
# subtype_psychoHC                     2.889 (1.275)  2.265  .025 *   [ 0.375,  5.403]      0.156   0.141
# subtype_psychoNA-MDD                 1.161 (1.297)  0.895  .372     [-1.396,  3.719]      0.062   0.056
# duration_mean:subtype_psychoHC      -2.933 (1.305) -2.248  .026 *   [-5.505, -0.361]     -0.154  -0.140
# duration_mean:subtype_psychoNA-MDD  -1.189 (1.322) -0.899  .370     [-3.795,  1.417]     -0.062  -0.056
# ───────────────────────────────────────────────────────────────────────────────────────────────────────

### between-group contrast for somMot
emtrends(model_lm, pairwise ~ subtype_psycho, var="duration_mean")
# subtype_psycho duration_mean.trend   SE  df lower.CL upper.CL
# A-MDD          11.9 2.67 131     6.56    17.14
# HC             4.7 1.75 131     1.24     8.17
#
# contrast     estimate   SE  df t.ratio p.value
# (A-MDD) - HC     7.15 3.16 131   2.261  0.0254

interactions::sim_slopes(model_lm, johnson_neyman = F,
                         pred = "duration_mean", modx = "subtype_psycho")
# SIMPLE SLOPES ANALYSIS 
# 
# Slope of duration_mean when subtype_psycho = HC: 
#   
#   Est.   S.E.   t val.      p
# ------ ------ -------- ------
#   4.70   1.75     2.68   0.01
# 
# Slope of duration_mean when subtype_psycho = A-MDD: 
#   
#   Est.   S.E.   t val.      p
# ------- ------ -------- ------
#   11.85   2.67     4.43   0.00

###########################################
#
# difference between A-MDD and NA-MDD 
#
###########################################

## between non-agit and agit ------------------------------------------
### trough duration 
dat_use %>% filter(group != "HC") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "subtype_psycho", 
                 covariate = c("HAMD_wave1_A")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho")
# ────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD)    0.140 (0.060) 116 2.309  .023 *      0.458 [0.065, 0.852]
# ────────────────────────────────────────────────────────────────────────────────

### Number of peak
dat_use %>% filter(group != "HC") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "peak_n", 
                 between = "subtype_psycho", 
                 covariate = c("HAMD_wave1_A")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho")
# ─────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD)   -1.485 (0.610) 115 -2.433  .017 *   -0.483 [-0.877, -0.090]
# ─────────────────────────────────────────────────────────────────────────────────

### time delay in SMN
dat_use %>% filter(group != "HC") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "subtype_psycho", 
                 covariate = c("HAMD_wave1_A")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho")
# ─────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD)   -1.485 (0.607) 116 -2.445  .016 *   -0.486 [-0.879, -0.092]
# ─────────────────────────────────────────────────────────────────────────────────

###########################################
#
# difference between A-MDD and HC
#
###########################################
### trough duration
dat_use %>% filter(subtype_psycho != "NA-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "subtype_psycho", 
                 covariate = c("age","gender","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "none")
# ──────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ──────────────────────────────────────────────────────────────────────────
# HC - (A-MDD)    0.140 (0.062) 133 2.255  .026 *      0.456 [0.056, 0.855]
# ──────────────────────────────────────────────────────────────────────────

### number of peak
dat_use %>% filter(subtype_psycho != "NA-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "peak_n", 
                 between = "subtype_psycho", 
                 covariate = c("age","gender","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "none")
# ───────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────────
# HC - (A-MDD)   -1.441 (0.615) 133 -2.341  .021 *   -0.473 [-0.873, -0.073]
# ───────────────────────────────────────────────────────────────────────────

### Time delay in SMN
dat_use %>% filter(subtype_psycho != "NA-MDD") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "subtype_psycho", 
                 covariate = c("age","gender","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "none")
# ──────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ──────────────────────────────────────────────────────────────────────────
# HC - (A-MDD)    0.003 (0.011) 133 0.263  .793       0.053 [-0.346, 0.453]
# ──────────────────────────────────────────────────────────────────────────

######################################################
#
# correlation to agitation 
#
#####################################################

## time delay in somMot and agitation
model_lm <- dat_use %>% filter(subtype_psycho != "HC") %>%
  mutate_at(vars(td_SMN), ~ scale(., scale = F)) %>% # demean the time delay
  lm(HAMD_wave1_item9 ~ td_SMN + 
       age + gender + education + QC_bold_fd_mean, data = .)
bruceR::GLM_summary(model_lm) # check the results
# ────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p        [95% CI of β] r(partial) r(part)
# ────────────────────────────────────────────────────────────────────────────────────
# td_somMot        -0.200 (0.099) -2.030  .045 *   [-0.396, -0.005]     -0.188  -0.188
# age               0.020 (0.104)  0.187  .852     [-0.187,  0.226]      0.018   0.017
# gendermale        0.051 (0.093)  0.542  .589     [-0.134,  0.235]      0.051   0.050
# education        -0.040 (0.103) -0.386  .700     [-0.244,  0.164]     -0.036  -0.036
# QC_bold_fd_mean  -0.052 (0.098) -0.535  .594     [-0.246,  0.141]     -0.051  -0.050
# ────────────────────────────────────────────────────────────────────────────────────

#########################################################
#
# Result visualization
#
#########################################################
dat_use$subtype_psycho <- factor(dat_use$subtype_psycho, 
                                 levels = c("A-MDD","HC","NA-MDD"))

# ETS correlation ------------------------------------------------------
p_cor_peak <- dat_use %>% 
  ggplot(aes(x = duration_mean, y  = peak_mean)) +
  geom_point(size = 3, alpha = .6) +
  geom_smooth(method = "lm", color = "darkred", size = .9) +
  ylab("Peak height") + xlab("Peak duration") +
  facet_grid(.~subtype_psycho, scales = "free") +
  theme_ggstatsplot() +
  easy_text_size(13)
p_cor_peak
ggsave(filename = "outputs/Analysis3_corDot.png", dpi = 300, height = 5, width = 11)

# group difference in duration -----------------------------------------
f1 <- dat_use %>% .$subtype_psycho
f1 <- c(.7, 1.7, 2.7)[as.integer(factor(f1))]

p_diff_duration <- dat_use %>%
  ggplot(aes(y = duration_mean, x = subtype_psycho, fill = subtype_psycho)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(aes(x = f1), size =3, alpha = .6, position = position_jitter(.03)) +
  geom_signif(annotations = rep("*", times = 2), 
              textsize = 7, vjust = .7,
              y_position = c(5.8, 5.9), 
              xmin = c(1, 1), 
              xmax = c(2, 3),
              tip_length = 0) +
  scale_fill_lancet() +
  ylab("Trough duration") +
  theme_classic() +
  easy_text_size(13) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_diff_duration
ggsave(filename = "outputs/Analysis3_boxplot_duration.png", dpi = 300, 
       height = 5, width = 8)

# group difference in N peak -------------------------------------------
p_diff_npeak <- dat_use %>%
  ggplot(aes(y = peak_n, x = subtype_psycho, fill = subtype_psycho)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(aes(x = f1), size =3, alpha = .6, position = position_jitter(.03)) +
  geom_signif(annotations = rep("*", times = 2),
              textsize = 7, vjust = .7,
              y_position = c(55.8, 56.8),
              xmin = c(1, 1),
              xmax = c(2, 3),
              tip_length = 0) +
  scale_fill_lancet() +
  ylab("Number of Peaks") +
  theme_classic() +
  easy_text_size(13) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_diff_npeak
ggsave(filename = "outputs/Analysis3_boxplot_Npeak.png", dpi = 300, 
       height = 5, width = 8)

# group difference in Time delay in SMN -------------------------------------
p_diff_SMN <- dat_use %>%
  ggplot(aes(y = td_SMN, x = subtype_psycho, fill = subtype_psycho)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(aes(x = f1), size =3, alpha = .6, position = position_jitter(.03)) +
  geom_signif(annotations = rep("*", times = 1),
              textsize = 7, vjust = .7,
              y_position = c(0.12),
              xmin = c(1),
              xmax = c(3),
              tip_length = 0) +
  scale_fill_lancet() +
  # coord_flip() +
  ylim(c(-.18, 0.13)) +
  ylab("Time delay in SMN") +
  theme_classic() +
  easy_text_size(13) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_diff_SMN
ggsave(filename = "outputs/Analysis3_boxplot_tdSMN.png", dpi = 300, 
       height = 5, width = 8)

# correlation between TD in SMN and agitation ---------------------------------
p_cor_SMN <- dat_use %>%
  ggplot(aes(x = HAMD_wave1_item9, y = td_SMN)) +
  geom_point(size =3 , alpha = .6, position = position_jitter(.03)) +
  geom_smooth(method = "lm", color = "darkred") +
  xlab("Agitation score") + ylab("Time delay in SMN") +
  ylim(c(-0.14, 0.14)) +
  theme_classic() +
  easy_text_size(13)
p_cor_SMN
ggsave(filename = "outputs/Analysis3_corDot_SMNagitation.png", dpi = 300, 
       height = 5, width = 6)

#########################################################
#
# Preparing the materials for visualizing time delay in
# SMN with MATLAB, BrainNetViewer
#
#########################################################

## reshape the network annotation 
net_anna <- rio::import("inputs/Power264_Yeo7.xlsx")
net_anna["ROI.Name"] <- paste0("R",formatC(1:214, digits = 2, flag = "0"))
net_anna <- net_anna %>%
  rename(x.mni = "x", y.mni="y", z.mni="z")
net_anna$x.mni <- as.integer(net_anna$x.mni)
net_anna$y.mni <- as.integer(net_anna$y.mni)
net_anna$z.mni <- as.integer(net_anna$z.mni)

net_anna_plot <- net_anna %>% mutate(net_color = ifelse(network_new == "somMot", "cyan",
                                                        ifelse(network_new == "ATN","green",
                                                               ifelse(network_new == "DMN", "red",
                                                                      ifelse(network_new == "salience", "purple",
                                                                             ifelse(network_new == "visual", "blue",
                                                                                    ifelse(network_new == "FPN", "orange", "black"))))))) %>%
  rename(network = "network_new") 

## add time delay 
path_td <- "inputs/timeseries_Power/time_lag_estimation/"
sbj_info <- rio::import('inputs/Analysis1_subject_table.xlsx')
file_names <- list.files(path_td, pattern = 'sub-[0-9]*_projection_map_weighted', full.names = T)
# for HC
sbj_hc <- sbj_info %>% filter(group == "HC") %>% .$participant_id 
file_use <- file_names[grep(paste(sbj_hc, collapse = "|"), file_names)]
dat_td_raw <- map_dfr(data.frame(file_use), ObtainBrainData_individual, .progress = T)
dat_td_SMN <- dat_td_raw[,-1] %>% colMeans()

net_anna_plot_this <- net_anna_plot %>% mutate(td = dat_td_SMN) %>%
  filter(network == "somMot") %>%
  mutate(cluster = 1) %>%
  select(1:3,11,12)
readr::write_delim(net_anna_plot_this, 
                   file = "outputs/SMN_Timedelay_HC.node", col_names = F, delim = "\t")

# for A-MDD
sbj_hc <- sbj_info %>% filter(subtype_psycho == "A-MDD") %>% .$participant_id 
file_use <- file_names[grep(paste(sbj_hc, collapse = "|"), file_names)]
dat_td_raw <- map_dfr(data.frame(file_use), ObtainBrainData_individual, .progress = T)
dat_td_SMN <- dat_td_raw[,-1] %>% colMeans()

net_anna_plot_this <- net_anna_plot %>% mutate(td = dat_td_SMN) %>%
  filter(network == "somMot") %>%
  mutate(cluster = 1) %>%
  select(1:3,11,12)
readr::write_delim(net_anna_plot_this, 
                   file = "outputs/SMN_Timedelay_AMDD.node", col_names = F, delim = "\t")

# for NA-MDD
sbj_hc <- sbj_info %>% filter(subtype_psycho == "NA-MDD") %>% .$participant_id 
file_use <- file_names[grep(paste(sbj_hc, collapse = "|"), file_names)]
dat_td_raw <- map_dfr(data.frame(file_use), ObtainBrainData_individual, .progress = T)
dat_td_SMN <- dat_td_raw[,-1] %>% colMeans()

net_anna_plot_this <- net_anna_plot %>% mutate(td = dat_td_SMN) %>%
  filter(network == "somMot") %>%
  mutate(cluster = 1) %>%
  select(1:3,11,12)
readr::write_delim(net_anna_plot_this, 
                   file = "outputs/SMN_Timedelay_NAMDD.node", col_names = F, delim = "\t")
