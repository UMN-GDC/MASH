library(tidyverse)
library(ggpattern)


# Read data
files <- dir("/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/Het/Results/5000/", 
             recursive = T, pattern = "*.csv", full.names = T)

df <- files %>%
  map(read_csv) %>%
  reduce(rbind)


# Annotate Estimators with family
df <- df %>%
  mutate(Estimator = str_replace(Estimator, "_est", ""),
         Estimator = factor(Estimator, levels = c("GCTA", "nGCTA", "SWD", "Combat", "AdjHE_FE", "nAdjHE", "AdjHE_RE")),
         Type = case_when(grepl("AdjHE", Estimator) ~ 'AdjHE',
                          grepl("GCTA", Estimator) ~ 'GCTA',
                          grepl("SWD", Estimator) ~ "Combat",
                          grepl("Combat", Estimator) ~ "Combat"
         ))
# Dataframe for simulation sigmas
sgs <- data.frame(sg = rep(c(1e-06, 0.25, 0.5), each = 2),
                  ss = rep(c(1e-06, 0.25), 3))


# Visualization theme
Custom_vis <- function(df) {
  df %>%
    unite("Group", c(nsites, Estimator, sg, ss), sep = "_", remove = F) %>%
    ggplot(aes(x = Estimator, y = Estimate, fill = Type, group = Group, pattern= nsites)) +
    geom_boxplot_pattern() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_hline(aes(yintercept = 0.5)) +
    xlab("") +
    ylab(expression(paste("Heritability Estimate", (hat(h^{2}))))) +
    guides(fill=guide_legend(title="Estimation Type"),
           pattern = guide_legend(title = "# clusters"))
  
}

# check which varies
apply(df, MARGIN =2 , function(x){length(unique(x))} )


# Filter for simulation
FILTER <- function(df) {
  df %>%
    filter(nclusts == 5, site_het == "TRUE") 
}

# Single panel
df %>%
  FILTER() %>%
  filter(sg == 0.5, ss == 0.25) %>%
  Custom_vis()

# Multipanel
##############################################
# Like this one for pulication


df %>% 
  FILTER() %>%
  Custom_vis() +
  facet_grid(rows = vars(ss), cols = vars(sg),
             labeller = labeller(.rows = label_both, .cols = label_both))




  





