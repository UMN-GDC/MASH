library(tidyverse)


apply(df, 2, function(x) length(unique(x)))

  

df %>% 
  filter(site_comp == "EQUAL", theta_alleles== 0.9) %>%
  ggplot(aes(x = Estimator, y = Estimate)) +
  geom_boxplot() +
  facet_grid(cols = vars(sg), rows = vars(ss), scales = "free_y") +
  theme(axis.text.x = element_text(angle = 90))
  


df %>% 
  filter(sg == 0.5, ss == 0.25, Estimator == "AdjHE_RE") %>%
  ggplot(aes(x = theta_alleles, y = Estimate, group = theta_alleles)) +
  geom_boxplot() + 
  facet_wrap(~ site_comp) +
  theme(axis.text.x = element_text(angle = 90))



# For fewer sites
df <- read_csv("Research/Stat_gen/tools/Basu_herit/Simulations/Sim_3_sites.csv")

df %>%
  ggplot(aes(x = Estimator, y = Estimate)) +
  geom_boxplot() +
  facet_grid(rows = vars(ss), cols = vars(sg)) + 
  theme(axis.text.x = element_text(angle = 90))



# Like this one for pulication
df <- read_csv("Research/Stat_gen/tools/Basu_herit/Simulations/Sim_explore3.csv")

df %>%
  ggplot(aes(x = Estimator, y = Estimate)) +
  geom_boxplot() +
  facet_grid(rows = vars(ss), cols = vars(sg)) + 
  theme(axis.text.x = element_text(angle = 90))

