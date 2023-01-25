library(tidyverse)

size = c(500, 1000, 2500, 5000)
AdjHE = c(0.05, 0.08, 0.77, 5.6)
AdjHERE = c(0.08, 0.33, 4.66, 35)
GCTA = c(0.29, 1.42, 17, 109)

df <- data.frame(size, AdjHE, AdjHERE, GCTA)


df %>%
  pivot_longer(cols = c("AdjHE", "AdjHERE", "GCTA"), names_to = "Estimator", values_to = "time (s)") %>%
  ggplot(aes(x= size, y = `time (s)`, col = Estimator)) +
  geom_line() + 
  geom_point() +
  scale_y_log10() +
  # scale_x_log10() +
  xlab("Sample size")

