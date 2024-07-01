#################################
#
# This function used to illustrate 
# the correlation matrix
#

ShowCorMat <- function(result_mat, lower = -1, upper = 1, title_text = "Correlation Matrix"){
  library(ggeasy)
  # define the color palette
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  viridis_palette <- viridis::viridis(50)
  # Create a mask to hide upper triangular elements
  mask <- upper.tri(result_mat)
  # Set upper triangular elements to NA
  result_mat[mask] <- NA
  diag(result_mat) <- NA
  cor_df <- melt(result_mat)
  # heatmap 
  p_heat <- ggplot(cor_df, aes(Var1, Var2, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = jet.colors(100), limits = c(lower, upper),
                         name = "Correlation", na.value = NA) +
    theme_blank() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title_text) +
    easy_remove_y_axis(what = c("title", "ticks", "text")) +
    easy_remove_x_axis(what = "title") + easy_center_title() +
    easy_text_size(13)
  
  return(p_heat)
}

ShowCorMat_whole <- function(result_mat, lower = -1, upper = 1, title_text = "Correlation Matrix"){
  library(ggeasy)
  # define the color palette
  jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
  viridis_palette <- viridis::viridis(50)

  diag(result_mat) <- NA
  cor_df <- melt(result_mat)
  
  # heatmap 
  p_heat <- ggplot(cor_df, aes(Var1, Var2, fill = value)) +
    geom_tile(color = NA) +
    scale_fill_gradientn(colors = jet.colors(100), limits = c(lower, upper),
                         name = "Correlation", na.value = NA) +
    theme_blank() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = title_text) +
    easy_remove_y_axis(what = c("title", "ticks", "text")) +
    easy_remove_x_axis(what = c("title","text")) + easy_center_title() +
    easy_text_size(13)
  
  return(p_heat)
}