library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)



# Code base (Black).
# Interfaces (Orange outline)
# Visualizers (Purple outline)
# Markdown-formatters (Blue outline)


analysis <- grViz(diagram = "digraph flowchart {
  node[shape=parallelogram]
  Covs[label='Covariates \\n (.csv)']
  Phens[label='Phenotypes \\n (.phen)']
  Eigs[label='Eigenvectors \\n (.eigenvec)']
  Estimator[label='Select Estimator']
  Selector[label='PC + cov + mpheno \\n
  choice ']

  node[shape=rectangle]
  load[label='load Phenos and \\n PCs/covars if specified \\n 
  O: Joined dataframe',color = Black, style = filled, fontcolor=  White]


  Loop[label='Loop over Covs?']
  List[label='List of feature sets (X,y) \\n to estimate over']
  Features[label='X, y']
  Estimate[label='Estimate heritability \\n O: h2, SE, computer resrouces, X ,y']
  
  
  Covs, Phens, Eigs -> load 
  load, Estimator, Selector -> Loop -> List-> Features -> Estimate
  Estimate -> List
}
")

analysis
analysis %>% export_svg %>% charToRaw %>% rsvg_png("analysis.png")


####################################################################
## Data loading ####################################################
####################################################################
loading <- grViz(diagram = "digraph flowchart {
  node[shape=parallelogram]
  infiles[label='Input files']
  infile[label='Input file']

  node[shape=diamond, color = Black, style = filled, fontcolor=  White]
  filetype[label='Parse filetype \\n (.csv, .txt, .phen, .eigenvec, files,']

  node[shape=rectangle, color = Black, style = filled, fontcolor=  White]
  easytype[label='.csv, .txt, .phen, .eigenvec']
  clean[label='Add IID, FID']

  
  hardtype[label='.files']
  readnii[label='Read .nii files \\n extract upper triange \\n stack subjects']
  id_hard[label= 'Add IID, ???FID???, task']
  
  merge[label='pd.merge data']

  infiles -> infile -> filetype -> easytype, hardtype
  hardtype -> readnii -> id_hard
  easytype -> clean -> infile
  
  id_hard -> infile
  
  clean, id_hard -> merge
  
}
")

loading

loading %>% export_svg %>% charToRaw %>% rsvg_png("loading.png")
