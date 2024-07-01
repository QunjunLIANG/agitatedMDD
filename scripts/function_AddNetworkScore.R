##################################
#
# This function is used to add the 
# network score to the TDp
#
# Net_annotation is the dataframe obtained
# from ObtainNetID
#

AddNetworkScore_Schaefer <- function(dat_TDp, net_anna){
  dat_tmp <- dat_TDp
  net_annotation <- net_anna
  network_name <- unique(net_annotation$network) %>% as.character()
  network_pos <- net_annotation$network
  
  dat_tmp <- dat_tmp %>% 
    # Mean
    mutate(visual_mean = rowMeans(across(all_of(1 + which(network_pos=="visual"))))) %>%
    mutate(SMN_mean = rowMeans(across(all_of(1 + which(network_pos=="somMot"))))) %>%
    mutate(dorsalAttn_mean = rowMeans(across(all_of(1 + which(network_pos=="dorsalAttn"))))) %>%
    mutate(salience_mean = rowMeans(across(all_of(1 + which(network_pos=="salience"))))) %>%
    mutate(limbic_mean = rowMeans(across(all_of(1 + which(network_pos=="limbic"))))) %>%
    mutate(control_mean = rowMeans(across(all_of(1 + which(network_pos=="control"))))) %>%
    mutate(DMN_mean = rowMeans(across(all_of(1 + which(network_pos=="DMN"))))) %>%
    mutate(TD_mean =  rowMeans(across(all_of(2:401)))) %>%
    # SD
    rowwise() %>%
    mutate(visual_sd = sd(across(all_of(1 + which(network_pos=="visual"))))) %>%
    mutate(SMN_sd = sd(across(all_of(1 + which(network_pos=="somMot"))))) %>%
    mutate(dorsalAttn_sd = sd(across(all_of(1 + which(network_pos=="dorsalAttn"))))) %>%
    mutate(salience_sd = sd(across(all_of(1 + which(network_pos=="salience"))))) %>%
    mutate(limbic_sd = sd(across(all_of(1 + which(network_pos=="limbic"))))) %>%
    mutate(control_sd = sd(across(all_of(1 + which(network_pos=="control"))))) %>%
    mutate(DMN_sd = sd(across(all_of(1 + which(network_pos=="DMN"))))) %>%
    mutate(TD_sd =  sd(across(all_of(2:401)))) 
  
  return(dat_tmp)
}

AddNetworkScore_Power <- function(dat_TDp, net_annotation){
  dat_tmp <- dat_TDp
  net_id <- net_annotation
  
  dat_tmp <- dat_tmp %>%
    mutate(td_SMN = rowMeans(across(all_of(1+which(net_id$network_new == "somMot"))), na.rm = T)) %>%
    mutate(td_visual = rowMeans(across(all_of(1+which(net_id$network_new == "visual"))), na.rm = T)) %>%
    mutate(td_salience = rowMeans(across(all_of(1+which(net_id$network_new == "salience"))), na.rm = T)) %>%
    mutate(td_ATN = rowMeans(across(all_of(1+which(net_id$network_new == "ATN"))), na.rm = T)) %>%
    mutate(td_FPN = rowMeans(across(all_of(1+which(net_id$network_new == "FPN"))), na.rm = T)) %>%
    mutate(td_DMN = rowMeans(across(all_of(1+which(net_id$network_new == "DMN"))), na.rm = T))
    
    return(dat_tmp)
}