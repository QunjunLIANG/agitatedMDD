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
dat_sbj <- rio::import('inputs/subject_merge_MDD_HC.xlsx')

# fix the education 
dat_sbj <- dat_sbj %>% mutate(educations = ifelse(education == 1, 'Illiterate',
                                                  ifelse(education == 2, 'Primary education',
                                                         ifelse(education == 3, 'Junior high school',
                                                                ifelse(education == 4, 'Senior high school',
                                                                       ifelse(education == 5, 'Undergraduate', 'Graduate'))))))

# summary the clinical scales
dat_mdd <- dat_sbj %>% filter(group == "MDD") 
dat_hc <- dat_sbj %>% filter(group == "HC")

#############################################################
#
# calculate the factor and total score
#
#############################################################

# define the items
item_depression <- c(1,2,7,8,13)
item_anxiety <- c(15,9,10,11,16)
item_sleepness <- c(4,5,6)

item_A <- c(1,2,7,8,10,13)
item_B <- c(4,5,6,9,11,12,14,15,17)
item_C <- c(3,16)

# calculate the factor and total score
hamd_wave1 <- dat_mdd %>% select(participant_id, contains('HAMD_wave1_item')) 
hamd_wave1 <- hamd_wave1 %>%
  mutate(HAMD_wave1_total = rowSums(hamd_wave1[,2:ncol(hamd_wave1)])) %>%
  mutate(HAMD_wave1_depression = rowSums(hamd_wave1[, 1+item_depression])) %>%
  mutate(HAMD_wave1_anxiety = rowSums(hamd_wave1[, 1+item_anxiety])) %>%
  mutate(HAMD_wave1_sleepness = rowSums(hamd_wave1[, 1+item_sleepness])) %>%
  mutate(HAMD_wave1_A = rowSums(hamd_wave1[, 1+item_A])) %>%
  mutate(HAMD_wave1_B = rowSums(hamd_wave1[, 1+item_B])) %>%
  mutate(HAMD_wave1_C = rowSums(hamd_wave1[, 1+item_C]))

hamd_wave2 <- dat_mdd %>% select(participant_id, contains('HAMD_wave2_item')) 
hamd_wave2 <- hamd_wave2 %>%
  mutate(HAMD_wave2_total = rowSums(hamd_wave2[,2:ncol(hamd_wave2)])) %>%
  mutate(HAMD_wave2_depression = rowSums(hamd_wave2[, 1+item_depression])) %>%
  mutate(HAMD_wave2_anxiety = rowSums(hamd_wave2[, 1+item_anxiety])) %>%
  mutate(HAMD_wave2_sleepness = rowSums(hamd_wave2[, 1+item_sleepness])) %>%
  mutate(HAMD_wave2_A = rowSums(hamd_wave2[, 1+item_A])) %>%
  mutate(HAMD_wave2_B = rowSums(hamd_wave2[, 1+item_B])) %>%
  mutate(HAMD_wave2_C = rowSums(hamd_wave2[, 1+item_C]))

hamd_wave3 <- dat_mdd %>% select(participant_id, contains('HAMD_wave3_item')) 
hamd_wave3 <- hamd_wave3 %>%
  mutate(HAMD_wave3_total = rowSums(hamd_wave3[,2:ncol(hamd_wave3)])) %>%
  mutate(HAMD_wave3_depression = rowSums(hamd_wave3[, 1+item_depression])) %>%
  mutate(HAMD_wave3_anxiety = rowSums(hamd_wave3[, 1+item_anxiety])) %>%
  mutate(HAMD_wave3_sleepness = rowSums(hamd_wave3[, 1+item_sleepness])) %>%
  mutate(HAMD_wave3_A = rowSums(hamd_wave3[, 1+item_A])) %>%
  mutate(HAMD_wave3_B = rowSums(hamd_wave3[, 1+item_B])) %>%
  mutate(HAMD_wave3_C = rowSums(hamd_wave3[, 1+item_C]))

# calculate HAMA wave1 total score
hama_wave1 <- dat_mdd %>% select(participant_id, contains('HAMA_wave1_item')) 
hama_wave1 <- hama_wave1 %>% mutate(HAMA_wave1_total = rowSums(hama_wave1[,2:ncol(hama_wave1)]) )

# merger the factor to the main data table
dat_mdd_factor <- dat_mdd %>% left_join(hamd_wave1) %>% left_join(hamd_wave2) %>%
  left_join(hamd_wave3) %>% left_join(hama_wave1)

dat_mdd_factor <- dat_mdd_factor %>% filter(HAMD_wave1_total >= 17)

##########################################################
#
# divide patients into agited and non-agited MDD
#
##########################################################

## threshold as >1 
dat_mdd_factor_treat_psycho <- dat_mdd_factor_treat %>% 
  mutate(psycho = ifelse(HAMD_wave1_item9 > 1, "A-MDD","NA-MDD"))

## subgroup difference in agitated 
bruceR::MANOVA(dat_mdd_factor_treat_psycho,
               subID = "participant_id", dv = "HAMD_wave1_total",
               between = "psycho")
# ────────────────────────────────────────────────────────────────────────────────
# MS    MSE df1 df2      F     p     η²p [90% CI of η²p]  η²G
# ────────────────────────────────────────────────────────────────────────────────
# subtype_psycho  878.233 30.209   1 117 29.072 <.001 ***   .199 [.102, .303] .199
# ────────────────────────────────────────────────────────────────────────────────

## subgroup difference in agitated 
bruceR::MANOVA(dat_mdd_factor_treat_psycho,
               subID = "participant_id", dv = "HAMD_wave1_item8",
               between = "psycho")
# ─────────────────────────────────────────────────────────────────────────────
# MS   MSE df1 df2     F     p     η²p [90% CI of η²p]  η²G
# ─────────────────────────────────────────────────────────────────────────────
# subtype_psycho  10.206 1.300   1 117 7.849  .006 **    .063 [.011, .146] .063
# ─────────────────────────────────────────────────────────────────────────────

## subgroup difference in agitated 
bruceR::MANOVA(dat_mdd_factor_treat_psycho,
               subID = "participant_id", dv = "HAMA_wave1_total",
               between = "psycho")
# ────────────────────────────────────────────────────────────────────────────────
# MS    MSE df1 df2      F     p     η²p [90% CI of η²p]  η²G
# ────────────────────────────────────────────────────────────────────────────────
# subtype_psycho  739.745 43.312   1 117 17.079 <.001 ***   .127 [.048, .225] .127
# ────────────────────────────────────────────────────────────────────────────────


##############################################################
#
# Export the new subject information
#
##############################################################
# merge mdd and HC
dat_merge <- data.table::rbindlist(list(dat_mdd_factor_treat_psycho, dat_hc), fill = T)
dat_merge <- dat_merge %>% 
  mutate(subtype_psycho=ifelse(group == 'MDD', psycho, group))

# summarise the data
dat_merge$educations <- factor(dat_merge$educations, 
                               levels = c("Illiterate","Primary education", 
                                          "Junior high school", "Senior high school",
                                          "Undergraduate", "Graduate"))

table1::table1(~ age + gender + educations 
               | group, 
               data = dat_merge)


table1::table1(~ age + gender + educations +
                 HAMD_wave1_total +
                 HAMD_wave1_A +
                 HAMD_wave1_B +
                 HAMD_wave1_C + 
                 HAMD_wave1_item8 +
                 HAMD_wave1_item3 +
                 HAMD_wave1_sleepness +
                 HAMA_wave1_total| subtype_psycho, 
               overall = F, extra.col=list("P-value"=pvalue),
               data = dat_merge %>% filter(group == "MDD"))

table1::table1(~ age + gender + educations 
               | subtype_psycho, 
               overall = F, extra.col=list("P-value"=pvalue),
               data = dat_merge %>% filter(subtype_psycho != "NA-MDD"))


rio::export(dat_merge, file = "inputs/Analysis1_subject_table.xlsx")

