library(tidyverse)

########################### Load data
df <- read_csv("/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/Full_results2.csv") %>%
  # Annotate Estimators with family
  mutate(Estimator = factor(Estimator, levels = c("GCTA", "nGCTA", "SWD", "Combat", "nAdjHE", "AdjHE", "AdjHE_RE")),
         `Study Type` = case_when(grepl("AdjHE", Estimator) ~ 'Single',
                          grepl("GCTA", Estimator) ~ 'Single',
                          grepl("SWD", Estimator) ~ "Meta",
                          grepl("Combat", Estimator) ~ "Meta"
         ),
         `Study Type` = ifelse(grepl("nAdjHE", Estimator), "Meta", `Study Type`),
         `Study Type` = ifelse(grepl("nGCTA", Estimator), "Meta", `Study Type`),
         Estimate = ifelse(
           (Estimator %in% c("AdjHE_RE", "AdjHE", "nAdjHE")) & (nsites == 25) & (nclusts == 1), Estimate * 4/3, Estimate),
         Estimate = ifelse(
           (!grepl("GCTA", Estimator)) & (nclusts==5) & (nsites == 2), Estimate * 4/3, Estimate),
         `Study Type` = factor(`Study Type`, levels = c("Single", "Meta"))) %>%
  unite("Group", c(nsites, Estimator, sg, ss), sep = "_", remove = F)




# apply(df, MARGIN = 2, unique)


df %>% 
  filter(sg == 0.5, ss == 0.25,
         nsubjects == 2000,
         prop_causal== 0.02, 
         nSNPs == 20000,
         site_comp == "EQUAL",
         site_het == "FALSE", 
         nsites != 1,
         # nclusts == 1, 
         # nsites == 1,
         # nnpc == "0"
         ) %>%
  ggplot(aes(y = Estimator, x = Estimate, fill= `Study Type`)) +
  geom_boxplot() +
  geom_vline(xintercept =  0.66) +
  facet_grid(rows = vars(nclusts), cols = vars(nsites), scales = "free",
             labeller = labeller(.rows = label_both, .cols = label_both)) +
  ylab("") +
  xlab(expression(paste("Heritability Estimate", (hat(h^{2}))))) +
  xlim(0,1)
  