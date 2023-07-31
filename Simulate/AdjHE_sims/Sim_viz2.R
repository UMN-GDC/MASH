library(tidyverse)

# GCCTA keep 3 pcs
# SWD, Combat, RE 2 pc
# AdjHE 1 pcs

########################### Load data
df <- read_csv("/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/Het_5000_pc_fixed.csv") %>%
  rename(Estimator = variable, Estimate = value) %>%
  mutate(Estimator = factor(Estimator, levels = c("GCTA", "nGCTA", "SWD", "Combat", "nAdjHE", "AdjHE", "AdjHE_RE")),
         `Study Type` = case_when(grepl("AdjHE", Estimator) ~ 'Single',
                                  grepl("GCTA", Estimator) ~ 'Single',
                                  grepl("SWD", Estimator) ~ "Single",
                                  grepl("Combat", Estimator) ~ "Single"
         ),
         `Study Type` = ifelse(grepl("nAdjHE", Estimator), "Meta", `Study Type`),
         `Study Type` = ifelse(grepl("nGCTA", Estimator), "Meta", `Study Type`),
         `Study Type` = factor(`Study Type`, levels = c("Single", "Meta"))) %>%
    unite("Group", c(nsites, Estimator, sg, ss), sep = "_", remove = F) %>%
  drop_na()




# apply(df, MARGIN = 2, unique)


df %>% 
  ggplot(aes(y = Estimator, x = Estimate, fill= `Study Type`)) +
  geom_boxplot() +
  geom_vline(xintercept =  0.66) +
  facet_grid(rows = vars(nclusts), cols = vars(nsites), scales = "free",
             labeller = labeller(.rows = label_both, .cols = label_both)) +
  ylab("") +
  xlab(expression(paste("Heritability Estimate", (hat(h^{2}))))) +
  xlim(0,1) +
  ggtitle("IID subpopulations: Heteroskedastic")

ggsave("Adding_sites_n_clusts_het.png", 
       dpi= 600)



df <- read_csv("5k_SWD_COMBAT_0.csv")
df1<- read_csv("5k_SWD_COMBAT_1.csv")

rbind(df,df1) %>%
  pivot_longer(cols = c(SWD, Combat), names_to = "Estimator", values_to = "heritability") %>%
  filter(case_when(
    (Distribution == "EQUAL") & (Sites == "S2") ~ between(heritability, 0.35,0.99),
    (Distribution == "EQUAL") & (Sites == "S25") ~ between(heritability, 0, 1),
    T ~ between(heritability, 0, 1))) %>% 
  ggplot(aes(x = Estimator, y = heritability)) + 
  geom_boxplot() +
  facet_grid(rows = vars(Sites), cols = vars(Distribution)) +
  geom_hline(yintercept = 0.66)

ggsave("SWD_Combat_0.png")
