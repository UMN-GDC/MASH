library(tidyverse)

df <- read.csv("SimProfilesDiffering.csv")

df %>%
  select(c(EstCombined, Est1, Est2)) %>%
  pivot_longer(everything()) %>%
  # filter(0.15 < value & value < 0.85) %>%
  ggplot(aes(x = value, group = name, fill = name)) + 
  geom_density(bw = 0.075, alpha = 0.25 ) + 
  scale_x_continuous(limits = c(0.2, 0.8))
  # 
  
ggsave("SimProfilesDiffering.png", width = 6, height = 4, units = "in", dpi = 300)
