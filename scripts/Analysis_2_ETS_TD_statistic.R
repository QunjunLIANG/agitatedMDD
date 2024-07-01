#############################################################################
#
# Temporal feature abnormality for Agitated MDD 
#
# For this pipeline, we tested the temporal structure among agitated MDD, HC 
# and non-agitated MDD. 
#
# The alternative features include:
#   1. Peak duration (ETS)
#   2. Time delay in SMN (ETS)
#   3. brain-symptom correlation (TD)
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
# 1 difference of peak duration
#
#################################################

dat_use %>% filter(group != "HC") %>% # between agitation subtypes
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "subtype_psycho", covariate = c("gender","age","education","HAMA_wave1_total")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD)    0.131 (0.061) 112 2.134  .035 *      0.435 [0.031, 0.839]
# ────────────────────────────────────────────────────────────────────────────────

dat_use %>% filter(subtype_psycho != "NA-MDD") %>% # between AMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "subtype_psycho", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ──────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ──────────────────────────────────────────────────────────────────────────
# HC - (A-MDD)    0.140 (0.062) 133 2.255  .026 *      0.456 [0.056, 0.855]
# ──────────────────────────────────────────────────────────────────────────

dat_use %>% filter(subtype_psycho != "A-MDD") %>% # between NAMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "duration_mean", 
                 between = "subtype_psycho", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - HC   -0.017 (0.054) 163 -0.311  .756      -0.055 [-0.407, 0.296]
# ────────────────────────────────────────────────────────────────────────────

p.adjust(c(.035, .026, .756), method = "fdr") # p-value adjust
# 0.0525 0.0525 0.7560

###########################################
#
# 2 difference in SMN time delay
#
###########################################

dat_use %>% filter(group != "HC") %>% # between agitat ion subtypes
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "subtype_psycho", covariate = c("gender","age","education","HAMA_wave1_total")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD)    0.027 (0.010) 112 2.689  .008 **     0.548 [0.144, 0.952]
# ────────────────────────────────────────────────────────────────────────────────

dat_use %>% filter(subtype_psycho != "NA-MDD") %>% # between AMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "subtype_psycho", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ──────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ──────────────────────────────────────────────────────────────────────────
# HC - (A-MDD)    0.003 (0.011) 133 0.263  .793       0.053 [-0.346, 0.453]
# ──────────────────────────────────────────────────────────────────────────

dat_use %>% filter(subtype_psycho != "A-MDD") %>% # between NAMDD and HC
  bruceR::MANOVA(subID = "participant_id", dv = "td_SMN", 
                 between = "subtype_psycho", covariate = c("gender","age","education")) %>%
  bruceR::EMMEANS(effect = "subtype_psycho", p.adjust = "fdr")
# ───────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df     t     p     Cohen’s d [95% CI of d]
# ───────────────────────────────────────────────────────────────────────────
# (NA-MDD) - HC    0.015 (0.010) 163 1.551  .123       0.276 [-0.075, 0.627]
# ───────────────────────────────────────────────────────────────────────────

p.adjust(c(.008, .793, .123), method = "fdr") # p-value adjust
# 0.0240 0.7930 0.1845

######################################################
#
# 3 correlation to agitation 
#
#####################################################

## time delay in somMot and agitation
model_lm <- dat_use %>% filter(subtype_psycho != "HC") %>%
  mutate_at(vars(td_SMN), ~ scale(., scale = F)) %>% # demean the time delay
  lm(HAMD_wave1_item9 ~ td_SMN + 
       age + gender + HAMA_wave1_total + education, data = .)
bruceR::GLM_summary(model_lm) # check the results
# ─────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p        [95% CI of β] r(partial) r(part)
# ─────────────────────────────────────────────────────────────────────────────────────
# td_SMN            -0.235 (0.088) -2.679  .009 **  [-0.408, -0.061]     -0.245  -0.228
# age               -0.164 (0.104) -1.577  .118     [-0.370,  0.042]     -0.147  -0.134
# gendermale         0.058 (0.085)  0.677  .500     [-0.111,  0.226]      0.064   0.058
# HAMA_wave1_total   0.421 (0.092)  4.566 <.001 *** [ 0.238,  0.604]      0.396   0.388
# education         -0.106 (0.096) -1.113  .268     [-0.296,  0.083]     -0.105  -0.095
# ─────────────────────────────────────────────────────────────────────────────────────

## trough duration and agitation
model_lm <- dat_use %>% filter(subtype_psycho != "HC") %>%
  mutate_at(vars(td_SMN), ~ scale(., scale = F)) %>% # demean the time delay
  lm(HAMD_wave1_item9 ~ duration_mean + 
       age + gender + HAMA_wave1_total + education, data = .)
bruceR::GLM_summary(model_lm) # check the results
# ────────────────────────────────────────────────────────────────────────────────────
# β    S.E.      t     p       [95% CI of β] r(partial) r(part)
# ────────────────────────────────────────────────────────────────────────────────────
# duration_mean     -0.118 (0.090) -1.309  .193     [-0.296, 0.060]     -0.123  -0.114
# age               -0.129 (0.106) -1.217  .226     [-0.339, 0.081]     -0.114  -0.106
# gendermale         0.056 (0.087)  0.637  .525     [-0.117, 0.229]      0.060   0.055
# HAMA_wave1_total   0.388 (0.094)  4.136 <.001 *** [ 0.202, 0.573]      0.364   0.360
# education         -0.095 (0.098) -0.974  .332     [-0.289, 0.099]     -0.092  -0.085
# ────────────────────────────────────────────────────────────────────────────────────