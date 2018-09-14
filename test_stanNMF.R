library(httr)
library(dplyr)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = min(3,parallel::detectCores()))

## test Stan model implementation

## ---- load TCGA-BLCA data ---- 

# read data (extracted from ftp://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/mutational_catalogs/exomes/Bladder/)
d <- read.table('sample_data/Bladder_exomes_mutational_catalog_96_subs.txt', sep = "\t", header = T, as.is = T)
m <- as.matrix(d[,-1])
rownames(m) <- d$Mutation.Type


## ---- download COSMIC signature probabilities ----

cosmic_content <- httr::content(httr::GET('http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt'), as = 'text', encoding = 'UTF-8')
cosmic_sigs <- read.table(textConnection(cosmic_content), sep = '\t', header = T)
cosmic_sigs <- cosmic_sigs %>%
  dplyr::select(-starts_with('X'))


blca_sig_ids <- c(1,2,5,10,13)
blca_sig_names <- paste('Signature',blca_sig_ids, sep = '.')

init_sigs <- cosmic_sigs %>% 
  dplyr::select(`Somatic.Mutation.Type`, one_of(blca_sig_names)) %>%
  dplyr::arrange(`Somatic.Mutation.Type`)

Winit <- as.matrix(init_sigs %>% dplyr::select(-Somatic.Mutation.Type))
rownames(Winit) <- init_sigs$Somatic.Mutation.Type

## confirm that column ordering of Winit & X are identical
assertthat::are_equal(rownames(Winit), rownames(m))

## ---- prep data for input to stan -----

###
# int<lower=1> N;  // num patients
# int<lower=1> M;  // num trinucleotides
# int<lower=1> K;  // num signatures
# real<lower=0> Wconc;        // concentration parameter on template priors - a large number means relative certainty, stick close to the prior
# real<lower=0> Hconc;        // concentration parameter on activations - a large number encourages sparsity
# matrix<lower=0>[M,K] Winit; // prior spectral templates -- these should each be different (or else it'd be unidentifiable). avoid zeroes too.
# matrix<lower=0>[M,N] X;     // spectrogram
###
stan_data <- list(
  'N' = dim(m)[2],
  'M' = dim(m)[1],
  'K' = dim(Winit)[2],
  'Wconc' = 1,
  'Hconc' = 10,
  'Winit' = Winit + 10e-5, ## add small value since stan model doesn't like exactly 0 values
  'X' = m
)

## ---- test stan model ---- 

test_is <-  rstan::stan('nmf_is_transcribe.stan', data = stan_data, iter = 10, chains = 1, init_r = 0.2)
stanNMFis <- rstan::stan('nmf_is_transcribe.stan', data = stan_data, iter = 1000, chains = 4)

test_plca <-  rstan::stan('nmf_plca_transcribe.stan', data = stan_data, iter = 10, chains = 1, init_r = 0.2)
stanNMFplca <- rstan::stan('nmf_plca_transcribe.stan', data = stan_data, iter = 1000, chains = 4)


