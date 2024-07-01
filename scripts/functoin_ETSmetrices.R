#############################################
#
# These functions are used to estimate the 
# duration between peak and the hight after 
# ETS rss estimation

ETSmetrices <- function(rss_rile){
  
  sbj <- stringr::str_extract(rss_rile, pattern = 'sub-[0-9]*')
  
  rss_raw <- readr::read_csv(rss_rile, col_names = "rss", show_col_types = FALSE)
  rss <- rss_raw$rss
    
  # initial vectors
  t_peak <- c()
  t_trough <- c()
  # detect the peak and trough in loop
  for (i in 2:(length(rss)-1)) {
    if (rss[i] > rss[i-1] & rss[i] > rss[i+1]) {
      t_peak <- c(t_peak, i) # obtain the position
    }
    if (rss[i] < rss[i-1] & rss[i] < rss[i+1]) {
      t_trough <- c(t_trough, i)
    }
  }
  # calculate the difference between trough
  peak_duration <- diff(t_trough)
  peak_high <- rss[t_peak]

  res_table <- data.frame(
  participant_id = sbj,
  
  duration_n = length(peak_duration),
  duration_mean = mean(peak_duration), 
  duration_sd = sd(peak_duration),
  
  peak_n = length(peak_high),
  peak_mean = mean(peak_high), 
  peak_sd = sd(peak_high)
  )
    
  return(res_table)
}
