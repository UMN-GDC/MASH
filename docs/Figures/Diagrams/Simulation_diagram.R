library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)



# Code base (Black).
# Interfaces (Orange outline)
# Visualizers (Purple outline)
# Markdown-formatters (Blue outline)

# Basic diagram

grViz("
 digraph dot {
 
Idea -> Execute -> Profit
 }
", engine = "dot")

# Saving simple graph
graph_to_save <- grViz("
 digraph dot {
 
Idea -> Execute -> Profit
 }
", engine = "dot")

graph_to_save %>% export_svg %>% charToRaw %>% rsvg_png("Research/Stat_gen/tools/Basu_herit/docs/Simulation_diagram.png")




# More complex graph
grViz("
 digraph dot {
 
Idea[label = 'Have a \n good idea', shape = black]
Exceute[label= 'Execute well', shape = circle, color = blue]
profit[label=  'Make money', color = green]
something[label= 'something else']

something, Idea -> Exceute -> profit, something



 }
", engine = "dot")





grViz("
 digraph dot {
 Ancest[shape = square, color = black, 
       label = 'Anc SNP (s)  p@_{s} ~ Uni(0.1, 0.9)']
       
 Subpops[label = 'Subpop(k) specific SNPs freqs p@_{s,k} ~ Beta(p@_{s,k}(1-&theta;@_{k})/&theta;@_{k}, (1-p@_{s,k})(1-&theta;@_{k})/&theta;@_{k} ']
 Genes[label='SNPs for subject (i) in pop (k), x@_{s,i}~ Bin(2, p@_{s,k},)']
 
 loci[label='m causal loci selected from genome']
 gene_effects[label='Effects, u~N(0, h@^{2} /m(I@_{m}))']
 
 sites_assign[label='assigned equally to nsites']
 site_contrib[label='X@_{s} ~ N(0, f@_{s}EmpVar(Xu))']

 model[label='Y=Xu + X@_{s}u@_{s} + &epsilon;']
       
 Ancest -> Subpops -> Genes -> loci
 gene_effects -> loci
 
 sites_assign -> site_contrib
 
 loci, site_contrib -> model
 
 
 
 }
", engine = "dot")




simulation <- grViz("
 digraph dot {
 node[shape = rectangle, color = blue]
Ancest[label = 'Common ancestral \n prototype']
Subpops[label = 'Subpopulation \n prototypes']
Genes[label='Subjects \n genotypes']

node[color = red]
Site_pops[label = 'Site \n Subpopulations']
Subj_pop[label= 'Subject \n subpopulation']
 
  node[shape = rectangle, color = black]
loci[label='Causal SNPs \n and effects']
 site_contrib[label='Site effects']

 model[label='Phenotype']
     
       
 Ancest -> Subpops -> Genes 
 Site_pops -> Subj_pop 
 
 
 Subj_pop -> Genes
 
 loci, site_contrib, Genes-> model
 

 }
")





simulation %>% export_svg %>% charToRaw %>% rsvg_png("Research/Stat_gen/tools/Basu_herit/docs/Simulation_diagram.png")


####################
timeline <- grViz("
 digraph {
 graph[layout = dot, rankdir = LR]
Estim[label= 'Finish \n Estimator \n June']
git[label = 'Establish Github repo \n June']
stand[label= 'Standardize \n code and I/O \n August']
review[label = 'Team \n code review \n September']

Estim -> git -> stand -> review

review -> stand ->git

 }
")

timeline %>% export_svg %>% charToRaw %>% rsvg_png("Research/Stat_gen/tools/Basu_herit/docs/timeline_diagram.png")

