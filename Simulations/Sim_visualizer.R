library(tidyverse)
library(ggpattern)

# Filter settings for figure
nsubjects = 2000
nclusts = 5
# nsites
site_het == TRUE
site_comp== "EQUAL"

  
# Load data
df <- read_csv("/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/Full_results.csv") %>%
  # Annotate Estimators with family
  mutate(Estimator = factor(Estimator, levels = c("GCTA", "nGCTA", "SWD", "Combat", "AdjHE", "nAdjHE", "AdjHE_RE")),
         Type = case_when(grepl("AdjHE", Estimator) ~ 'AdjHE',
                          grepl("GCTA", Estimator) ~ 'GCTA',
                          grepl("SWD", Estimator) ~ "Combat",
                          grepl("Combat", Estimator) ~ "Combat"
         ),
         Estimate =  ifelse(Type != "GCTA", Estimate *4/3, Estimate),
         Estimate = ifelse(Type == "GCTA" & ss == 1e-06, Estimate * 4/3, Estimate),
         nsites = factor(nsites))

# check which varies
apply(df, MARGIN =2 , function(x){length(unique(x))} )


### Visualization theme

hlines<- data.frame(sg = rep(c(1e-06, 0.25, 0.5), 2),
                    ss = rep(c(1e-06, 0.25), each = 3)) %>%
  mutate(h2 = c(0, 1/3, 2/3, 0, 1/3, 2/3))

Custom_vis <- function(df) {
  df %>%
    # het filter
    # filter(site_het == TRUE, nsubjects == 5000, nclusts == 5, nsites == 25)
    # other filter
    filter(site_het == TRUE, nsubjects == 2000, nclusts == 5, site_comp== "EQUAL") %>%
    mutate(nsites = droplevels(nsites)) %>%
    unite("Group", c(nsites, Estimator, sg, ss), sep = "_", remove = F) %>%
    ggplot(aes(x = Estimator, y = Estimate, fill = Type, group = Group, pattern= nsites)) +
    geom_boxplot_pattern() +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_hline(data = hlines, aes(yintercept = h2)) +
    xlab("") +
    ylab(expression(paste("Heritability Estimate", (hat(h^{2}))))) +
    guides(fill=guide_legend(title="Adjustment Type"),
           pattern = guide_legend(title = "# sites")) +
    ylim(0, 1)
  
}



df %>%
  FILTER() %>%
  # Uncomment for Single panel
  # filter(sg == 0.5, ss == 0.25) %>%
  Custom_vis() +
  # Uncomment for multi panel
  facet_grid(rows = vars(ss), 
             cols = vars(sg),
             labeller = labeller(.rows = label_both, .cols = label_both)) +
  ggtitle("IID site composition: 5 clusters")

ggsave("/home/christian/Research/Stat_gen/tools/Basu_herit/docs/Site_Extension_paper/Figures/IID_C5.png", dpi= 600)



df %>%
  FILTER() %>%
  group_by(Type, sg, ss, nsites, Estimator) %>% 
  summarize(m = median(Estimate), s = sd(Estimate)) %>%
  ggplot(aes(x = Estimator, y = m, ymin = m+s, ymax = m-s, color = nsites)) +
  geom_point(position = position_dodge())+ 
  geom_errorbar(position = position_dodge()) +
  facet_grid(rows = vars(ss), 
             cols = vars(sg),
             labeller = labeller(.rows = label_both, .cols = label_both)) + 
  geom_hline(data = hlines, aes(yintercept = h2))
  
  


  





