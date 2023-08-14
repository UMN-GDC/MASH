library(tidyverse)

df <- read.csv("5k_SWD_COMBAT_0EQUAL.csv")

df %>%
  filter(Sites == "S2") %>%
  select(-c(Sites, Distribution, Rep)) %>%
  pivot_longer(everything()) %>%
    


  
