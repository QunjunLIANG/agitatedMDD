#############################################################################
#
# Plot results 
#
#
# Liang Qunjun 2024/06/13

library(tidyverse)
library(ggeasy)
library(ggsignif)
library(ggsci)
library(patchwork)
# load the data
dat_use <- rio::import("inputs/Analysis2_collect_fMRI.xlsx")

# define the cut-off 
cut_hama <- 24
cut_hamd <- 24

#########################################################
#
# Main findings Result visualization
#
#########################################################
dat_use$subtype_psycho <- factor(dat_use$subtype_psycho, 
                                 levels = c("A-MDD","HC","NA-MDD"))
text_size <- 18
plot_height <- 4
plot_width <- 6

# group difference in duration -----------------------------------------
p_diff_duration <- dat_use %>%
  ggplot(aes(y = duration_mean, x = subtype_psycho, fill = subtype_psycho)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(size =3, alpha = .6, position = position_jitter(.03)) +
  geom_signif(annotations = rep("*", times = 2), 
              textsize = 7, vjust = .7,
              y_position = c(5.8, 5.9), 
              xmin = c(1, 1), 
              xmax = c(2, 3),
              tip_length = 0) +
  scale_fill_tron() +
  ylab("Trough duration (a.u.)") +
  theme_classic() +
  easy_text_size(text_size) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_diff_duration
ggsave(plot = p_diff_duration, 
       filename = "outputs/Fig2A.png", dpi = 300, 
       height = plot_height, width = plot_width)
ggsave(plot = p_diff_duration, 
       filename = "outputs/Fig2A.tiff", dpi = 300, 
       height = plot_height, width = plot_width)

# group difference in peak height -----------------------------------------
p_diff_peak <- dat_use %>%
  ggplot(aes(y = peak_mean, x = subtype_psycho, fill = subtype_psycho)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(size =3, alpha = .6, position = position_jitter(.03)) +
  scale_fill_tron() +
  ylab("Peak height (a.u.)") +
  ylim(c(175, 205)) +
  theme_classic() +
  easy_text_size(text_size) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_diff_peak
ggsave(filename = "outputs/Fig2B.png", dpi = 300, 
       height = plot_height, width = plot_width)
ggsave(filename = "outputs/Fig2B.tiff", dpi = 300, 
       height = plot_height, width = plot_width)

# group difference in Time delay in SMN -------------------------------------
p_diff_SMN <- dat_use %>%
  ggplot(aes(y = td_SMN, x = subtype_psycho, fill = subtype_psycho)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(size =3, alpha = .6, position = position_jitter(.03)) +
  geom_signif(annotations = rep("*", times = 1),
              textsize = 7, vjust = .7,
              y_position = c(0.12),
              xmin = c(1),
              xmax = c(3),
              tip_length = 0) +
  scale_fill_tron() +
  # coord_flip() +
  ylim(c(-.18, 0.13)) +
  ylab("Time delay in SMN (a.u.)") +
  theme_classic() +
  easy_text_size(text_size) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_diff_SMN
ggsave(filename = "outputs/Fig2D.png", dpi = 300, 
       height = plot_height, width = plot_width)
ggsave(filename = "outputs/Fig2D.tiff", dpi = 300, 
       height = plot_height, width = plot_width)

# correlation between TD in SMN and agitation ---------------------------------
p_cor_SMN <- dat_use %>%
  ggplot(aes(x = HAMD_wave1_item9, y = td_SMN)) +
  geom_point(size =3 , alpha = .6, position = position_jitter(.03)) +
  geom_smooth(method = "lm", color = "darkred") +
  xlab("Agitation score (a.u.)") + ylab("Time delay in SMN (a.u.)") +
  ylim(c(-0.14, 0.14)) +
  theme_classic() +
  easy_text_size(text_size)
p_cor_SMN
ggsave(filename = "outputs/Fig2E.png", dpi = 300, 
       height = plot_height, width = 6)
ggsave(filename = "outputs/Fig2E.tiff", dpi = 300, 
       height = plot_height, width = plot_width)

#########################################################
#
# control analysis result visualization - anxiety level
#
#########################################################
dat_hama <- dat_use %>% select(starts_with("HAMA"))
dat_use <- dat_use %>% mutate(HAMA_wave1_total = rowSums(dat_hama))
dat_use <- dat_use %>% 
  mutate(anxiety_group = ifelse(group == "HC", "HC",
                                ifelse(HAMA_wave1_total >= cut_hama, "severe", "moderated")))

p_control_duration <- dat_use %>%
  ggplot(aes(y = duration_mean, x = anxiety_group, fill = anxiety_group)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(size =3, alpha = .6, position = position_jitter(.03)) +
  ylab("Trough duraiton (a.u.)") + ggtitle("Anxiety subgroup") +
  scale_fill_tron() +
  theme_classic() +
  easy_text_size(text_size) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_control_duration

p_control_tdSMN <- dat_use %>%
  ggplot(aes(y = td_SMN, x = anxiety_group, fill = anxiety_group)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(size =3, alpha = .6, position = position_jitter(.03)) +
  ylab("Time delay in SMN (a.u.)") + ggtitle("Anxiety subgroup") +
  scale_fill_tron() +
  theme_classic() +
  easy_text_size(text_size) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_control_tdSMN

#########################################################
#
# control analysis - severity subgroups
#
#########################################################
dat_use <- dat_use %>%
  mutate(severity = ifelse(group=="HC","HC",
                           ifelse(HAMD_wave1_total >= cut_hamd,"severe","moderated")))

p_severe_duration <- dat_use %>%
  ggplot(aes(y = duration_mean, x = severity, fill = severity)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(size =3, alpha = .6, position = position_jitter(.03)) +
  scale_fill_tron() +
  ylab("Trough duration (a.u.)") + ggtitle("Depression subgroup") +
  theme_classic() +
  easy_text_size(text_size) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_severe_duration

p_severe_tdSMN <- dat_use %>%
  ggplot(aes(y = td_SMN, x = severity, fill = severity)) +
  geom_boxplot(size =1 , width = .4, alpha = .6, outlier.colour = NA) +
  geom_point(size =3, alpha = .6, position = position_jitter(.03)) +
  scale_fill_tron() +
  ylab("Time delay in SMN (a.u.)") + ggtitle("Depression subgroup") +
  theme_classic() +
  easy_text_size(text_size) + easy_remove_x_axis(what = "title") +
  easy_remove_legend()
p_severe_tdSMN

p_control_all <- p_control_duration + p_control_tdSMN + p_severe_duration + p_severe_tdSMN +
  patchwork::plot_annotation(tag_levels = "A")

p_control_all
ggsave(plot = p_control_all, filename = "outputs/FigS3.png",
       height = 8, width = 10, dpi = 300)
ggsave(plot = p_control_all, filename = "outputs/FigS3.tiff",
       height = 8, width = 10, dpi = 300)
