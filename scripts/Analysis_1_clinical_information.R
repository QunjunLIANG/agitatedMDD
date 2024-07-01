##############################################################################
#
# Subject selection
#
# In this script, we selected the data of the subjects who is suitable in 
# this study.
#
# Liang Qunjun 2023/11/13

library(tidyverse)
library(NbClust)
library(ggiraphExtra)
library(ggsci)
library(ggeasy)
library(tidyLPA)
library(psych)
library(mclust)
library(ggsignif)
library(patchwork)
library(scales)
source("scripts/function_PvalueForTable1.R")

# load the data
dat_sbj <- rio::import('inputs/Analysis1_subject_table.xlsx')
dat_sbj$psycho[which(is.na(dat_sbj$psycho))] <- "HC"
dat_hama <- dat_sbj %>% select(starts_with("HAMA"))
dat_sbj <- dat_sbj %>% mutate(HAMA_wave1_total = rowSums(dat_hama))

# demographic test ==================================================
bruceR::MANOVA(dat_sbj,
               subID = "participant_id", dv = "age",
               between = "psycho")%>%
  bruceR::EMMEANS(effect = "psycho", p.adjust = "fdr")
# ─────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────
# HC - (A-MDD)         -0.525 (2.212) 210 -0.237  .813      -0.043 [-0.484, 0.397]
# (NA-MDD) - (A-MDD)   -2.326 (2.299) 210 -1.012  .506      -0.192 [-0.650, 0.266]
# (NA-MDD) - HC        -1.802 (1.875) 210 -0.961  .506      -0.149 [-0.522, 0.225]
# ─────────────────────────────────────────────────────────────────────────────────
# Pooled SD for computing Cohen’s d: 12.108
# P-value adjustment: FDR method for 3 tests.

gender_table <- table(dat_sbj$psycho, dat_sbj$gender)
kruskal.test(gender ~ psycho, data = dat_sbj)
# Kruskal-Wallis chi-squared = 36.785, df = 2, p-value = 1.029e-08
PMCMRplus::kwAllPairsDunnTest(x = dat_sbj$gender, 
                              g = as.factor(dat_sbj$psycho), p.adjust.method = "fdr")
# Pairwise comparisons using Dunn's all-pairs test
#          A-MDD   HC     
# HC     4.5e-06     -      
# NA-MDD 0.89    2.7e-07

# test the difference of clinical signatures ========================
## subgroup difference in HAMD total score 
bruceR::MANOVA(dat_sbj,
               subID = "participant_id", dv = "HAMD_wave1_total",
               between = "psycho") %>%
  bruceR::EMMEANS(effect = "psycho")
# ─────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD)   -5.628 (1.044) 117 -5.392 <.001 *** -1.024 [-1.400, -0.648]
# ─────────────────────────────────────────────────────────────────────────────────

dat_sbj %>% filter(group == "MDD") %>%
  select(participant_id, gender, age, psycho, starts_with("HAMD_wave1_item")) %>%
  pivot_longer(cols = c(5:21), values_to = "score", names_to = "item") %>%
  bruceR::MANOVA(subID = "participant_id", dv = "score",
                 between = "psycho", within = "item",sph.correction="GG") %>%
  bruceR::EMMEANS(effect = "psycho", by = "item", p.adjust = "fdr")
# Pairwise Comparisons of "psycho":
# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Contrast            "item" Estimate    S.E.  df       t     p     Cohen’s d [95% CI of d]
# ────────────────────────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD) HAMD_wave1_item1    -0.055 (0.128) 117  -0.427  .670     -0.049 [-0.278,  0.179]
# (NA-MDD) - (A-MDD) HAMD_wave1_item10   -0.480 (0.169) 117  -2.839  .005 **  -0.433 [-0.736, -0.131]
# (NA-MDD) - (A-MDD) HAMD_wave1_item11   -0.371 (0.192) 117  -1.926  .057 .   -0.335 [-0.679,  0.009]
# (NA-MDD) - (A-MDD) HAMD_wave1_item12   -0.666 (0.135) 117  -4.919 <.001 *** -0.602 [-0.844, -0.359]
# (NA-MDD) - (A-MDD) HAMD_wave1_item13   -0.212 (0.133) 117  -1.593  .114     -0.192 [-0.430,  0.047]
# (NA-MDD) - (A-MDD) HAMD_wave1_item14   -0.345 (0.112) 117  -3.077  .003 **  -0.312 [-0.513, -0.111]
# (NA-MDD) - (A-MDD) HAMD_wave1_item15   -0.338 (0.192) 117  -1.764  .080 .   -0.305 [-0.648,  0.037]
# (NA-MDD) - (A-MDD) HAMD_wave1_item16    0.024 (0.140) 117   0.171  .864      0.022 [-0.229,  0.272]
# (NA-MDD) - (A-MDD) HAMD_wave1_item17    0.028 (0.123) 117   0.226  .821      0.025 [-0.195,  0.246]
# (NA-MDD) - (A-MDD) HAMD_wave1_item2    -0.306 (0.186) 117  -1.644  .103     -0.276 [-0.609,  0.056]
# (NA-MDD) - (A-MDD) HAMD_wave1_item3    -0.132 (0.243) 117  -0.542  .589     -0.119 [-0.553,  0.315]
# (NA-MDD) - (A-MDD) HAMD_wave1_item4     0.017 (0.123) 117   0.138  .891      0.015 [-0.205,  0.235]
# (NA-MDD) - (A-MDD) HAMD_wave1_item5    -0.028 (0.117) 117  -0.243  .808     -0.026 [-0.235,  0.184]
# (NA-MDD) - (A-MDD) HAMD_wave1_item6     0.063 (0.141) 117   0.447  .656      0.057 [-0.195,  0.309]
# (NA-MDD) - (A-MDD) HAMD_wave1_item7    -0.210 (0.139) 117  -1.505  .135     -0.190 [-0.439,  0.060]
# (NA-MDD) - (A-MDD) HAMD_wave1_item8    -0.607 (0.217) 117  -2.802  .006 **  -0.548 [-0.935, -0.161]
# (NA-MDD) - (A-MDD) HAMD_wave1_item9    -2.009 (0.100) 117 -20.032 <.001 *** -1åå.814 [-1.994, -1.635]
# ────────────────────────────────────────────────────────────────────────────────────────────────────

## subgroup difference in HAMA total score 
bruceR::MANOVA(dat_sbj,
               subID = "participant_id", dv = "HAMA_wave1_total",
               between = "psycho") %>%
  bruceR::EMMEANS(effect = "psycho")
# ────────────────────────────────────────────────────────────────────────────────
# Contrast Estimate    S.E.  df      t     p     Cohen’s d [95% CI of d]
# ─────────────────────────────────────────────────────────────────────────────────
# (NA-MDD) - (A-MDD)   -5.165 (1.250) 117 -4.133 <.001 *** -0.785 [-1.161, -0.409]
# ─────────────────────────────────────────────────────────────────────────────────

# table for dataset information ========================

table1::table1(~ age + gender + educations 
               | group, 
               data = dat_sbj)

table1::table1(~ age + gender + educations +
                 HAMD_wave1_total +
                 HAMA_wave1_total| subtype_psycho, 
               overall = F, extra.col=list("P-value"=pvalue),
               data = dat_sbj %>% filter(group == "MDD"))

table1::table1(~ age + gender + educations 
               | subtype_psycho, 
               overall = F, extra.col=list("P-value"=pvalue),
               data = dat_sbj %>% filter(subtype_psycho != "NA-MDD"))

# plot the item-level difference of HAMD ---------------------------------

label_item <- c("Depressive mood","Feelings of guilt","Suicide",
                "Early insomnia","Middle insomnia","Late insomnia",
                "Work and activities","Retardation","Agitation",
                "Psychic anxiety","Somatic anxiety","Gas. symptoms",
                "General somatic", "genital symptom", "Hypochondriasis",
                "Weight loss","Insight")

p_HAMD <- dat_sbj %>% filter(group == "MDD") %>%
  select(participant_id, gender, age, psycho, starts_with("HAMD_wave1_item")) %>%
  pivot_longer(cols = c(5:21), values_to = "score", names_to = "item") %>%
  group_by(psycho, item) %>%
  summarize(score_mean = mean(score)) %>%
  mutate(item = factor(item, level = paste0("HAMD_wave1_item",1:17), labels = label_item)) %>%
  ggplot(aes(x = item, y = score_mean, group =psycho)) +
  geom_line(aes(color = psycho)) +
  geom_point(aes(color = psycho)) +
  geom_signif(
    annotations = c("**","***","**","***","**"),
    textsize = 7, vjust = .7,
    y_position = c(3.5,3.5,3.5,3.5,3.5),
    xmin = c(8,9,10,12,14),
    xmax = c(8,9,10,12,14),
    tip_length = 1
  ) +
  ylab("Score") +
  theme_classic() +
  easy_rotate_x_labels(angle = -90) +
  easy_text_size(14) +
  easy_remove_x_axis(what = "title") +
  easy_move_legend(to = "bottom") +
  easy_add_legend_title("Subgroup")

ggsave(plot = p_HAMD, filename = "outputs/Analysis1_HAMD_subgroup.png", 
       height = 5, width = 8)
ggsave(plot = p_HAMD, filename = "outputs/Analysis1_HAMD_subgroup.tiff", 
       height = 5, width = 8)
