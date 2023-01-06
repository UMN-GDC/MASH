library(tidyverse)
library(ggpattern)


# Read data
files <- dir(c("/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/Basic/Results/2000/",
               "/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/EQUAL/Results/2000/",
               "/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/IID/Results/2000/",
               "/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/Het/Results/5000/"
               ),
             recursive = T, pattern = "*.csv", full.names = T)

df <- files %>%
  map(read_csv) %>%
  reduce(rbind) %>%
  # Annotate Estimators with family
  mutate(Estimator = factor(Estimator, levels = c("GCTA", "nGCTA", "SWD", "Combat", "AdjHE", "nAdjHE", "AdjHE_RE")),
         Type = case_when(grepl("AdjHE", Estimator) ~ 'AdjHE',
                          grepl("GCTA", Estimator) ~ 'GCTA',
                          grepl("SWD", Estimator) ~ "Combat",
                          grepl("Combat", Estimator) ~ "Combat"
         ),
         Estimate = ifelse(Type != "GCTA", Estimate *4/3, Estimate))


# Dataframe for simulation sigmas
sgs <- data.frame(sg = rep(c(1e-06, 0.25, 0.5), each = 2),
                  ss = rep(c(1e-06, 0.25), 3))


# Visualization theme

hlines<- data.frame(sg = rep(c(1e-06, 0.25, 0.5), 2),
                    ss = rep(c(1e-06, 0.25), each = 3)) %>%
  mutate(h2 = c(0, 1/3, 2/3, 0, 1/3, 2/3))

Custom_vis <- function(df) {
  df %>%
    unite("Group", c(nsites, Estimator, sg, ss), sep = "_", remove = F) %>%
    ggplot(aes(x = Estimator, y = Estimate, fill = Type, group = Group, pattern= nsites)) +
    geom_boxplot_pattern() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_hline(data = hlines, aes(yintercept = h2)) +
    xlab("") +
    ylab(expression(paste("Heritability Estimate", (hat(h^{2}))))) +
    guides(fill=guide_legend(title="Estimation Type"),
           pattern = guide_legend(title = "# sites")) +
    ylim(0, 1)
  
}

# check which varies
apply(df, MARGIN =2 , function(x){length(unique(x))} )


# Filter for simulation
FILTER <- function(df) {
  df %>%
    # het filter
    # filter(site_het == TRUE, nsubjects == 5000, nclusts == 5, nsites == 25)
    # other filter
    filter(site_het == TRUE, nsubjects == 2000, nclusts == 5, nsites == 25)
}



df %>%
  FILTER() %>%
  # Uncomment for Single panel
  # filter(sg == 0.5, ss == 0.25) %>%
  Custom_vis() +
  # Uncomment for multi panel
  facet_grid(rows = vars(ss), cols = vars(sg),
             labeller = labeller(.rows = label_both, .cols = label_both)) +
  ggtitle("Site Heteroskedasticity")






  





