
## load BayesNMF functions
if (!require(mutsigNMF)) 
  devtools::install_github('jburos/mutsigNMF')
library(mutsigNMF)
library(tidyr)
library(ggplot2)
library(dplyr)
library(rstan)
library(httr)

## read data 
## (extracted from ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/mutational_catalogs/exomes/Bladder/)
d <- read.table('sample_data/Bladder_exomes_mutational_catalog_96_subs.txt', sep = "\t", header = T, as.is = T)
m <- as.matrix(d[,-1])
rownames(m) <- gsub(d$Mutation.Type, pattern = "[[:punct:]]", replacement = '', fixed = F) 

## test bayesNMF method
## as described here: https://www.broadinstitute.org/cancer/cga/msp
bayesNMF <- mutsigNMF::run_bayesNMF(mat = m, n_runs = 2, iter_per_run = 100)

## 
