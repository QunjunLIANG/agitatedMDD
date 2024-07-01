####################################
#
# This function used to extract intra-
# and inter- connectivity among SMN,
# VN and DMN

ExtractConnectivity <- function(file_name, net_annotation) {
  # define the position of each network
  id_smn <- which(net_annotation$network_new == "somMot")
  id_visual <- which(net_annotation$network_new == "visual")
  id_dmn <- which(net_annotation$network_new == "DMN")
  # load the correlation matrix and extract the connectivity
  mat_tmp <- readr::read_csv(file_name, col_names = F, show_col_types = F) %>% as.matrix()
  sbj_name <- stringr::str_extract(file_name, pattern = "sub-[0-9]*")
  
  intra_smn <- mat_tmp[id_smn,id_smn] %>% .[which(lower.tri(.))] %>% mean()
  intra_visual <- mat_tmp[id_visual,id_visual] %>% .[which(lower.tri(.))] %>% mean()
  smn_dmn <- mat_tmp[id_smn,id_dmn] %>% mean()
  smn_visual <- mat_tmp[id_smn,id_visual] %>% mean()
  global_fc <- mat_tmp[which(lower.tri(mat_tmp))] %>% mean
  
  result_dat <- data.frame(
    participant_id = sbj_name,
    globaFC = global_fc,
    intraSMN = intra_smn,
    intraVisual = intra_visual,
    SMN_DMN = smn_dmn,
    SMN_visual = smn_visual
  )
  
  return(result_dat)
}

ExtractSMNDMN <- function(file_name, net_annotation) {
  # define the position of each network
  id_smn <- which(net_annotation$network_new == "somMot")
  id_dmn <- which(net_annotation$network_new == "DMN")
  # load the correlation matrix and extract the connectivity
  mat_tmp <- readr::read_csv(file_name, col_names = F, show_col_types = F) %>% as.matrix()
  mat_tmp[id_smn,id_smn] <- 0
  mat_tmp[id_dmn,id_dmn] <- 0
  smn_dmn <- mat_tmp[c(id_smn,id_dmn),c(id_smn,id_dmn)] 
  
  return(smn_dmn)
}