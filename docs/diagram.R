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
  in1[label='File paths']
  in2[label='Pheno specification']
  out[label='Store results', shape = ellipse]
  
  node[shape=rectangle, color = Black, style = filled, fontcolor=  White]
  load[label='load Phenos and \\n PCs/covars if specified']
  subset[label='subset to 1 pheno \\n remove missing']
  predlmm[label='PredLMM']
  adjhe[label='AdjHE']
  

  select[label='select estimator', shape= diamond]
  
  
  
  in1 -> load
  in2 -> load
  
  load -> subset -> select -> predlmm -> out
  select->adjhe -> out -> subset
  
  
  
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

  

  node[shape=parallelogram]
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
