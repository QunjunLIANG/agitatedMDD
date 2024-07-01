#############################################################################
#
# Control analysis for Temporal feature abnormality for Agitated MDD 
#
# For this pipeline, we repeated the analysis of ETS and TD, but use difference
# group criteria. 
#
# The control analysis include:
#.  1. Peak height difference
#   2. MDD as a whole group
#   3. subgroup with MDD severity
#
# Liang Qunjun 2024/06/13

library(tidyverse)
library(ggstatsplot)
library(ggeasy)
library(ggsci)
library(patchwork)
library(ggsignif)
library(sjPlot)
# load the data
dat_use <- rio::import("inputs/Analysis2_collect_fMRI.xlsx")
dat_hama <- dat_use %>% select(starts_with("HAMA"))
dat_use <- dat_use %>% mutate(HAMA_wave1_total = rowSums(dat_hama))

#################################################
#
# 1 difference of peak height
#
#################################################

dat_use %>%  # difference between HC and MDD whole group
  bruceR::MANOVA(subID = "participant_id", dv = "peak_mean", 
                 between = "group", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "group", p.adjust = "fdr")
# ───────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────
# MDD - HC   -0.769 (0.923) 207 -0.833  .406      -0.130 [-0.439, 0.178]
# ───────────────────────────────────────────────────────────────────────

dat_use %>% filter(group != "HC") %>% # between agitation subtypes
  bruceR::MANOVA(subID = "participant_id", dv = "peak_mean", 
                 between = "subtype_psycho", covariate = c("gender","age","education","HAMA_wave1_total")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD)    1.236 (1.323) 112 0.934  .352       0.190 [-0.213, 0.594]
# ────────────────────────────────────────────────────────────────────────────────

dat_use %>% filter(subtype_psycho != "NA-MDD") %>% # between AMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "peak_mean", 
                 between = "subtype_psycho", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ──────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ──────────────────────────────────────────────────────────────────────────
# HC - (A-MDD)    1.937 (1.150) 133 1.683  .095 .     0.340 [-0.059, 0.740]
# ─────────────────────────────────────────────────────────────────────────

dat_use %>% filter(subtype_psycho != "A-MDD") %>% # between NAMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "peak_mean", 
                 between = "subtype_psycho", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - HC   -0.000 (1.006) 163 -0.000 1.000      -0.000 [-0.351, 0.351]
# ────────────────────────────────────────────────────────────────────────────

p.adjust(c(.352, .095, 1), method = "fdr") # p-value adjust
# 0.528 0.285 1.000

#################################################
#
# 1 control analysis: anxiety subgroup
#
#################################################

dat_anxiety <- dat_use %>% 
  mutate(anxiety_group = ifelse(group == "HC", "HC",
                                ifelse(HAMA_wave1_total > 21, "severe", "moderated")))

dat_anxiety %>%  # difference between HC and MDD whole group
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "anxiety_group", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "anxiety_group", p.adjust = "fdr")
# ─────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────
# moderated - HC       -0.057 (0.055) 206 -1.037  .452      -0.188 [-0.625, 0.249]
# severe - HC          -0.068 (0.056) 206 -1.201  .452      -0.222 [-0.667, 0.224]
# severe - moderated   -0.010 (0.058) 206 -0.179  .858      -0.034 [-0.489, 0.422]
# ─────────────────────────────────────────────────────────────────────────────────

dat_anxiety %>% # SMN td for HC and MDD whole group 
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "anxiety_group", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "anxiety_group", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────
# moderated - HC        0.003 (0.010) 206 0.267  .790       0.048 [-0.389, 0.485]
# severe - HC           0.014 (0.010) 206 1.399  .402       0.258 [-0.187, 0.703]
# severe - moderated    0.011 (0.010) 206 1.111  .402       0.210 [-0.246, 0.665]
# ────────────────────────────────────────────────────────────────────────────────

#################################################
#
# control analysis: MDD severity subgroup
#
#################################################

dat_use <- dat_use %>%
  mutate(severity = ifelse(group=="HC","HC",
                     ifelse(HAMD_wave1_total >= 24,"severe","moderated")))

# For trough duration ------------------------------------------------------
dat_use %>% filter(group != "HC") %>% # between agitation subtypes
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "severity", covariate = c("gender","age","education","HAMA_wave1_total")) %>%
  bruceR::EMMEANS(effect = "severity", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────
# severe - moderated    0.049 (0.068) 112 0.717  .475       0.158 [-0.279, 0.596]
# ────────────────────────────────────────────────────────────────────────────────

dat_use %>% filter(severity != "moderated") %>% # between AMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "severity", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "severity", p.adjust = "fdr")
# ──────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ──────────────────────────────────────────────────────────────────────────
# severe - HC   -0.063 (0.053) 173 -1.201  .231      -0.206 [-0.544, 0.132]
# ──────────────────────────────────────────────────────────────────────────#

dat_use %>% filter(severity != "severe") %>% # between NAMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "severity", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "severity", p.adjust = "fdr")
# ─────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────
# moderated - HC   -0.076 (0.068) 123 -1.110  .269      -0.246 [-0.685, 0.193]
# ─────────────────────────────────────────────────────────────────────────────

p.adjust(c(.475, .231, .269), method = "fdr") # p-value adjust
# 0.726 0.726 0.726

# For time delay in SMN ------------------------------------------------------
dat_use %>% filter(group != "HC") %>% # between agitation subtypes
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "severity", covariate = c("gender","age","education","HAMA_wave1_total")) %>%
  bruceR::EMMEANS(effect = "severity", p.adjust = "fdr")
# ─────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────
# severe - moderated   -0.020 (0.011) 112 -1.755  .082 .    -0.388 [-0.825, 0.050]
# ─────────────────────────────────────────────────────────────────────────────────

dat_use %>% filter(severity != "moderated") %>% # between AMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "severity", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "severity", p.adjust = "fdr")
# ─────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────
# severe - HC    0.003 (0.009) 173 0.293  .770       0.050 [-0.288, 0.388]
# ─────────────────────────────────────────────────────────────────────────

dat_use %>% filter(severity != "severe") %>% # between NAMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "severity", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "severity", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────
# moderated - HC    0.020 (0.012) 123 1.750  .083 .     0.388 [-0.051, 0.827]
# ────────────────────────────────────────────────────────────────────────────

p.adjust(c(.082, .771, .083), method = "fdr") # p-value adjust
# 0.1245 0.7710 0.1245

#######################################
#
#
# Brain-sympton Correlation of other aspects
#
#####################################

## correlation of anxiety dimension
model_lm <- dat_use %>% filter(subtype_psycho != "HC") %>%
  mutate_at(vars(td_SMN), ~ scale(., scale = F)) %>% # demean the time delay
  lm(HAMA_wave1_total ~ td_SMN + 
       age + gender  + education, data = .)
bruceR::GLM_summary(model_lm) # check the results
# ──────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p       [95% CI of β] r(partial) r(part)
# ──────────────────────────────────────────────────────────────────────────────
# td_SMN       0.117 (0.089)  1.325  .188     [-0.058, 0.293]      0.124   0.115
# age          0.434 (0.098)  4.435 <.001 *** [ 0.240, 0.628]      0.385   0.385
# gendermale  -0.030 (0.087) -0.345  .731     [-0.202, 0.142]     -0.032  -0.030
# education    0.169 (0.096)  1.760  .081 .   [-0.021, 0.360]      0.163   0.153
# ──────────────────────────────────────────────────────────────────────────────