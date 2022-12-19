library(tidyverse)
library(optparse)
 
option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="base directory of simulation results", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="Simulations/out.txt", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# Read data
files <- dir(opt$d, full.names = T, pattern = "*.csv", recursive = T)

df <- files %>% 
  map(read_csv) %>%
  reduce(rbind) %>%
  write_csv(opt$out)


# Remove the folder that used to hold the separate results
unlink(opt$d, recursive = T)

 
