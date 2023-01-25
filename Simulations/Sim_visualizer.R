library(tidyverse)
library(ggpattern)


########################### Load data
df <- read_csv("/home/christian/Research/Stat_gen/tools/Basu_herit/Simulations/Full_results2.csv") %>%
  # Annotate Estimators with family
  mutate(Estimator = factor(Estimator, levels = c("GCTA", "nGCTA", "SWD", "Combat", "nAdjHE", "AdjHE", "AdjHE_RE")),
         Type = case_when(grepl("AdjHE", Estimator) ~ 'Single',
                          grepl("GCTA", Estimator) ~ 'Single',
                          grepl("SWD", Estimator) ~ "Meta",
                          grepl("Combat", Estimator) ~ "Meta"
         ),
         Type = ifelse(grepl("nAdjHE", Estimator), "Meta", Type),
         Type = ifelse(grepl("nGCTA", Estimator), "Meta", Type),
         Estimate =  ifelse(Type != "GCTA", Estimate *4/3, Estimate),
         Estimate = ifelse(Type == "GCTA" & ss == 1e-06, Estimate * 4/3, Estimate),
         Estimate = ifelse(nsites ==1, Estimate *3/4, Estimate),
         nsites = factor(nsites)) %>%
  unite("Group", c(nsites, Estimator, sg, ss), sep = "_", remove = F)
  

# check which varies
apply(df, MARGIN =2 , function(x){length(unique(x))} )


######################### DEFAULT plotter
# Hetitability markers
oneline<- data.frame(sg = 0.5,
                    ss = 0.25,
                    h2 = 2/3)

hlines<- data.frame(sg = rep(c(1e-06, 0.25, 0.5), 2),
                    ss = rep(c(1e-06, 0.25), each = 3)) %>%
  mutate(h2 = c(0, 1/3, 2/3, 0, 1/3, 2/3))
hlines1 <- data.frame(sg = c(1e-06, 0.25, 0.5), 
                      h2 = c(1e-06, 0.25, 0.5))


#' Boxplot_template
#' The general framerwork for making consistet boxplots for the simulation studies
#'
#' @param df dataframe
#' @param x string specifying the variable along the x axis
#' @param y sting specifying the variable along the y axis
#' @param type sting Family of site adjustments
#' @param Group string specifying a discrete varaible speciying the smallest unit for the boxplot counting 
#' @param pattern sting specifying a discrete vaariable to make patterns over
#' @param comp sting specifying a discrete vaariable to make patterns over
#' @param het sting specifying whether error was heterskedastic
#' @param facet_col string specifying which variable to facet
#' @param facet_row string specifying which variable to facet
#' @param hlines Dataframe variances and heritabilty for drawing true line 
#'
#' @return
#' @export
#'
#' @examples
Estimator_boxplot <- function(df, x="Estimator", y= "Estimate", type= "Type" , Group = "Group", pattern = "nsites", comp= "EQUAL",
                              het =F, facet_col = NULL, facet_row = NULL, hlines=NULL) {

  # for faceting
  if (is.null(facet_row)) {
    facet_row = NULL
  } else {
    facet_row = vars(!!as.name(facet_row))
  }
  if (is.null(facet_col)) {
    facet_col = NULL
  } else {
    facet_col= vars(!!as.name(facet_col))
  }
  
  g <- df %>%
    ggplot(aes(x = !!as.name(x), y = !!as.name(y), fill = !!as.name(type), group = !!as.name(Group))) +
    theme(axis.text.x = element_text(angle = 90)) +
    geom_hline(data = hlines, aes(yintercept = h2)) +
    xlab("") +
    ylab(expression(paste("Heritability Estimate", (hat(h^{2}))))) +
    guides(fill=guide_legend(title="Study Type"), pattern = guide_legend(title = "# sites")) +
    ylim(0, 1)+
    {if(is.null(pattern)) geom_boxplot()} +
    {if(!is.null(pattern)) geom_boxplot_pattern(aes(pattern = !!as.name(pattern)))} +
    facet_grid(rows = facet_row, cols = facet_col,
      labeller = labeller(.rows = label_both, .cols = label_both))
  #print(g)
  g
}

figure_maker <- function(gg, Het, N) {
  if (het) {
    Het = "Heteroskedastic"
    Het2 = "Het"
  } else{
    Het = "Homoskedastic"
    Het2 = "Homo"
  }
  
  
  
  gg <- gg + 
    ggtitle(paste(comp, "site composition,", Het, "error: \n", "N=", N))
  
  prefix <- "/home/christian/Research/Stat_gen/tools/Basu_herit/docs/Figures/Estimates/Simulations/"
  fileout <- paste0("N", N, "_", Het2) 
  # fileout <- "N2000_C1_S1_HOMO"
  
  ggsave(paste0(prefix, fileout, ".png"), 
         dpi= 600)
  

}


g <- df %>%
  filter(nsubjects==2000, site_het == F, 
         sg == 0.5, ss ==0.25, se == 0.25,
         nSNPs == 20000, theta_alleles == 0.1, prop_causal == 0.02) %>%
  Estimator_boxplot(x="Estimator", y= "Estimate", type= "Type" , Group = "Group", pattern = "nsites", comp= "EQUAL",
                het =F, facet_col = "nsites", facet_row = "nclusts", hlines=oneline)
g
figure_maker(g, Het = F, N = 2000)



