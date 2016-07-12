
## load BayesNMF functions
source('functions/bayesnmf.signature_discovery.r')

library(tidyr)
library(ggplot2)
library(dplyr)

# read data (extracted from ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/mutational_catalogs/exomes/Bladder/)
d <- read.table('sample_data/Bladder_exomes_mutational_catalog_96_subs.txt', sep = "\t", header = T, as.is = T)
m <- as.matrix(d[,-1])
## replace punctuation in mutation descriptions, per 
rownames(m) <- gsub(d$Mutation.Type, pattern = "[[:punct:]]", replacement = '', fixed = F) 

#####################################################
############## BayesNMF parameters
############## n.iter = number of independent simulations
############## Kcol = number of initial signatures
############## tol = tolerance for convergence
############## a0 = hyper-parameter
n.iter <- 100000
Kcol <- 10
tol <- 1.e-07
a0 <- 10
res <- BayesNMF.L1KL(V0 = m, n.iter = n.iter, a0 = a0, tol = tol, K = Kcol, K0 = Kcol, phi = 1.0)
names(res) <- c('W','H','n.like','n.evid','n.lambda','n.error')
str(res['W'])

## inspect raw 'W' object returned
results <- 
  as.data.frame(res['W']) %>%
  dplyr::mutate(mutation_type = d$Mutation.Type,
                rowname = rownames(m)
                ) %>%
  tidyr::gather(signature, weight, starts_with('W.')) %>%
  dplyr::mutate(mutation_context = gsub(mutation_type, pattern = '(\\D)\\[\\D>\\D\\](\\D)', replacement = '\\1x\\2')
                , mutation_replacement = gsub(mutation_type, pattern = '(\\D)\\[(\\D>\\D)\\](\\D)', replacement = '\\2')
                )

## plot/summarize 'W'
ggplot(data = results, aes(x = mutation_context, y = weight, fill = mutation_replacement)) + 
  stat_summary(fun.y = identity, geom = 'bar', position = 'stack') + facet_grid(signature~mutation_replacement) +
  theme_minimal() + theme(axis.text.x=element_text(angle=90, hjust=1))

## identify K
summarize_k <- function(res) {
  lambda <- res[['n.lambda']]
  lambda <- lambda[[length(lambda)]]
  lambda.min <- min(lambda)
  lambda <- lambda-lambda.min
  lambda.norm <- lambda/sum(lambda)
  index <- lambda.norm > 0.01
  K <- sum(index)
  return(K)
}
print(summarize_k(res))

## re-create summary produced by example_usage.R code
summarize_results <- function(res = NULL,
                              W = res[['W']],
                              H = res[['H']],
                              lambda = res[['n.lambda']],
                              tumor.type = 'dataset'
) { 
  W1 <- W
  H1 <- H
  W.norm <- apply(W,2,function(x) x/sum(x))
  for (i in 1:ncol(W)) {
    H1[i,] <- H1[i,]*colSums(W)[i]
    W1[,i] <- W1[,i]*rowSums(H)[i]
  }
  df <- get.df.solution(W1, H, lambda, tumor.type)
  K <- length(table(df$class))
  return(plot.lego.observed.barplot(df,paste("Mutation Signatures in ",tumor.type,sep="")))
}
summarize_results(res, tumor.type = 'BLCA') +
  theme(axis.text.x=element_text(size = 8))



