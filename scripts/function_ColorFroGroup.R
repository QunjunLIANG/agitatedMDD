# color identification --------------------------------------------------------
#
ColorForGroup <- function(){
  # obtain the color representations of two group of subjects
  isfahan <- MetBrewer::met.brewer("Isfahan1")
  length(isfahan)
  isfahan[1]
  
  hc_color <- colorspace::lighten(isfahan[2], 0.35)
  mdd_color <-  colorspace::lighten(isfahan[6], 0.35)
  aSMN_color <- "#00468BFF"
  dSMN_color <- "#ED0000FF"
  
  return(data.frame(hc_col = hc_color, mdd_col = mdd_color,
                    aSMN_col = aSMN_color, dSMN_color = dSMN_color))
}
