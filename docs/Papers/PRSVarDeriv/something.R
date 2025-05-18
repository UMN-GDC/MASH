library(tidyverse)

data.frame(p_1 = c(1:100)/(200)) %>%
  mutate(y = (1/9) * (p_1- 1/4)^2 ) %>%
  ggplot(aes(x = p_1, y = y)) +
  geom_point()
