library(DiagrammeR)
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)



# Code base (Black).
# Interfaces (Orange outline)
# Visualizers (Purple outline)
# Markdown-formatters (Blue outline)


graph <- grViz(diagram = "digraph flowchart {
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

graph
# graph %>% export_svg %>% charToRaw %>% rsvg_png("diagram.png")


