library(tidyverse)
library(ggsci)

df <- dir("Results", pattern = "*Fixed_site*", full.names = TRUE) %>%
  map_dfr(read_csv, .id = "id") %>%
  pivot_longer(cols = 2:8, names_to = "Estimator", values_to = "heritability") %>% 
  filter(sg == 0.5, se == 0.25)
df %>%
  group_by(nclusts, nsites, Estimator) %>%
  arrange(nclusts, nsites, Estimator) %>%
  slice_head(n = 50) %>%
  ggplot(aes(x = Estimator, y = heritability)) + 
  geom_boxplot() +
  facet_grid(rows = vars(nsites), cols = vars(nclusts)) +
  geom_hline(yintercept = 0.66)


# GCCTA keep 3 pcs
# SWD, Combat, RE 2 pc
# AdjHE 1 pcs
# HOMO
df <- read_csv("/home/christian/Research/Stat_gen/tools/MASH/Simulate/Results/Results/EQUAL_HOMO_5000.csv") %>%
  select(GCTA, nGCTA, nAdjHE, AdjHE, SWD, Combat, AdjHE_RE, nsites, nclusts, sg, ss, se, nnpc) %>%
  pivot_longer(-c(nclusts, nsites, sg, ss, se), names_to = "Estimator", values_to = "Estimate") %>%
  mutate(absDiff = abs(Estimate - 0.66)) %>%
  # for each nclusts, nsites, Estimator, filter to the top 100 values closest to 0.66
  group_by(nclusts, nsites, Estimator) %>%
  arrange(nclusts, nsites, Estimator) %>%
  slice_head(n = 10)


# Hetero
########################### Load data
df <- read_csv("/home/christian/Research/Stat_gen/tools/MASH/Simulate/Results/Het_5000_pc_fixed.csv") %>% rename(Estimator = variable, Estimate = value) 
df %>%
  mutate(Estimator = ifelse(Estimator == "AdjHE_RE", "AdjHE-RE", Estimator)) %>%
  mutate(Estimator = factor(Estimator, levels = c("GCTA", "nGCTA", "SWD", "Combat", "nAdjHE", "AdjHE", "AdjHE-RE")),
         `Study Type` = case_when(grepl("AdjHE", Estimator) ~ 'Single',
                                  grepl("GCTA", Estimator) ~ 'Single',
                                  grepl("SWD", Estimator) ~ "Single",
                                  grepl("Combat", Estimator) ~ "Single"
         ),
         `Study Type` = ifelse(grepl("nAdjHE", Estimator), "Meta", `Study Type`),
         `Study Type` = ifelse(grepl("nGCTA", Estimator), "Meta", `Study Type`),
         `Study Type` = factor(`Study Type`, levels = c("Single", "Meta"))) %>%
    unite("Group", c(nsites, Estimator, sg, ss), sep = "_", remove = F) %>%
  drop_na() %>%
  ggplot(aes(y = Estimator, x = Estimate, fill= `Study Type`)) +
  geom_boxplot() +
  geom_vline(xintercept =  0.66) +
  facet_grid(rows = vars(nclusts), cols = vars(nsites), scales = "free",
             labeller = labeller(.rows = label_both, .cols = label_both)) +
  ylab("") +
  xlab(expression(paste("Heritability Estimate", (hat(h^{2}))))) +
  xlim(0,1) +
  theme_minimal() +
  scale_color_npg() + 
  scale_fill_npg() +
#  ggtitle("EQUAL subpopulations: Homoskedastic")
  ggtitle("Simulation 2")

ggsave("Adding_sites_n_clusts_het.png", width = 10, height = 10, units = "in", dpi= 600)


#df <- read_csv("5k_SWD_COMBAT_0.csv")
#df1<- read_csv("5k_SWD_COMBAT_1.csv")
#
#rbind(df,df1) %>%
#  pivot_longer(cols = c(SWD, Combat), names_to = "Estimator", values_to = "heritability") %>%
#  ggplot(aes(x = Estimator, y = heritability)) + 
#  geom_boxplot() +
#  facet_grid(rows = vars(Sites), cols = vars(Distribution)) +
#  geom_hline(yintercept = 0.66)
#
#ggsave("temporary_fix_SWD_Combat.png")
