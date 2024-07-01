library(ggplot2)
library(ggeasy)

# 设置频率和时间范围
frequency1 <- 0.1  # 0.1 Hz
frequency2 <- 0.03  # 0.03 Hz
time <- seq(0, 4 / frequency1, length.out = 1000)  # 两个周期，时间从0到4/frequency1

# 计算两个余弦曲线的值
cosine_values1 <- cos(2 * pi * frequency1 * time)
cosine_values2 <- cos(2 * pi * frequency2 * time)

# 创建数据框
df <- data.frame(time, cosine_values1, cosine_values2)

# 使用 ggplot2 绘制两个余弦曲线和它们的组合
p_freq <- ggplot(df, aes(x = time)) +
  geom_line(aes(y = cosine_values1), color = "blue", linetype = "dotted", size = 1, alpha = 0.8) +
  geom_line(aes(y = cosine_values2), color = "red", linetype = "dashed", size = 1, alpha = 0.8) +
  geom_line(aes(y = cosine_values1 + cosine_values2), color = "green", linetype = "solid", size = 1, alpha = 0.8) +
  labs(title = "Combined Cosine Curves (0.1 & 0.03 Hz)",
       x = "Time (seconds)",
       y = "Cosine Value") +
  theme_minimal() +
  easy_text_size(13)
p_freq

ggsave(plot = p_freq, filename = "outputs/supplement_Combnine_curve.png", 
       width = 13, height = 6)
