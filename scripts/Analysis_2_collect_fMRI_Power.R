############################################################################
#
# Collecting fMRI metrics
#
# This script is used to generated the metrics for the sequential analysis,
# including:
#   1. edge time series (ETS)
#   2. time delay
#
# Liang Qunjun 2023/11/11

library(tidyverse)
library(bruceR)
library(ggstatsplot)
library(ggridges)
library(psych)
library(brainconn)
library(RColorBrewer)
library(emmeans)
library(ggeasy)
library(ggsci)
library(cowplot)
library(scales)
library(ggsignif)
library(patchwork)
library(sjPlot)
source('scripts/function_ObtainNetID.R')
source('scripts/function_ObtainBrainData_individual.R')
source('scripts/function_AddNetworkScore.R')
source('scripts/function_EstPSMD.R')
source("scripts/functoin_ETSmetrices.R")

# indicate the path to the inputs and participant information
path_td <- "inputs/timeseries_Power/time_lag_estimation/"
path_ets <- "inputs/timeseries_Power/Edge_time_series_estimation/"
sbj_info <- rio::import('inputs/Analysis1_subject_table.xlsx')
net_anna <- rio::import("inputs/Power264_Yeo7.xlsx")
outfile_name <- "inputs/Analysis2_collect_fMRI.xlsx"

######################################################
#
# Collect all fMRI metrices
#
####################################################

#############   Collect for time delay   ##############
# load the data
file_names <- list.files(path_td, pattern = 'sub-[0-9]*_projection_map_weighted', full.names = T)
file_use <- file_names[grep(paste(sbj_info$participant_id, collapse = "|"), file_names)]
dat_td_raw <- map_dfr(data.frame(file_use), ObtainBrainData_individual, .progress = T)
colnames(dat_td_raw)[2:ncol(dat_td_raw)] <- paste0("R",formatC(1:214, width = 3, flag = "0"))
## add network-level mean TDp
dat_td_net <- AddNetworkScore_Power(dat_td_raw, net_anna)

##############   Collect for PSMD    ####################
# network identification
dat_psmd <- EstPSMD(dat_td_net, net_annotation = net_anna)

#################  Collect for ETS   ####################
# load the data
file_names <- list.files(path_ets, pattern = 'sub-[0-9]*_rms.csv', full.names = T)
file_use <- file_names[grep(paste(sbj_info$participant_id, collapse = "|"), file_names)]

dat_ets_raw <- data.frame()
for (i in file_use) {
  dat_tmp <- ETSmetrices(i)
  dat_ets_raw <- rbind(dat_ets_raw, dat_tmp)
}

############### merge all fMRI metrices ##################
dat_use <- sbj_info %>% left_join(dat_td_net) %>%
  left_join(dat_psmd) %>% left_join(dat_ets_raw)

######################################################
#
# Export to the outer
#
####################################################
rio::export(dat_use, file = outfile_name)
