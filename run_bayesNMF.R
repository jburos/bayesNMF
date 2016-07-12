
source('functions/bayesnmf.signature_discovery.r')

library(tidyr)
library(ggplot2)
library(dplyr)
library(purrr)

#' execute bayesNMF for signature deconvolution
#' 
#' @param mat (matrix) N*M matrix of ints with N mutation frequencies for each of M samples
#' @param tumor_name (chr) name of tumor, used for plotting headers, etc. default: 'tumor.type'
#' @param max_k (int) max hypothesized number of signatures to look for. default: 10
#' @param n_runs (int) number of times to run optimization code. default: 10
#' @param iter_per_run (int) number of iter for each run. default: 100000
#' @param tol (float) tolerance used to assess convergence
#' @param a0 (int) value of hyperparameter `a`. default: 10
#' @param phi (int) value of hyperparameter `phi`. default: 1 (appropriate for l1-loss with exponential prior)
#' @param prior (chr) either L1KL (expoential priors) or L2KL (half-normal priors). default: L1KL
#' @param hyper (logical) if TRUE, reduce the effect of hyper-mutant samples in the signature discovery. default: FALSE
#' 
#' @return list of objects
#' 
run_bayesNMF <- function(mat, tumor_name = 'tumor.type', max_k = 10, n_runs = 10,
                         iter_per_run = 100000, tol = 1.e-07, a0 = 10, phi = 1.0,
                         prior = 'L1KL', hyper = FALSE
                         ) {
  
  ## TODO check input mat
  lego96 <- mat
  
  ## check params
  if (!prior %in% c("L1KL","L2KL")) {
    stop(paste("prior must be one of L1KL or L2KL. Provided value was ",prior))
  }

  method <- prior

  # --- process hypermutations --- #
  if (hyper) {
    lego96 <- get.lego96.hyper(lego96)
    method <- paste(method,"hyper",sep=".")
  }


  # --- run the algorithm n.iter times --- #  
  results <- list()
  for (i in 1:n.iter) {
    if (prior=="L1KL") {
      results[[i]] <- BayesNMF.L1KL(as.matrix(lego96),iter_per_run,a0,tol,max_k,max_k,phi)
    } else {
      results[[i]] <- BayesNMF.L2KL(as.matrix(lego96),iter_per_run,a0,tol,max_k,max_k,phi)
    }
  }

  # --- summarize results --- #  
  res.WES <- summarize_results(results)
  
  ############## frequency figure 
  #pdf(file=paste(OUTPUT,paste(method,a0,"signature.freq.pdf",sep="."),sep=""),width=4,height=5)
  s1 <- 1.5
  s2 <- 2.0
  par(mfrow=c(1,1))
  par(mar=c(5,5,2,1))
  barplot(table(unlist(res.WES[['K']])),cex=s1,cex.axis=s1,cex.main=s1,cex.names=s1,cex.lab=s1,xlab="# of signatures",ylab="Freq.",main=paste(tumor.type,sep="."))

  ##########################################################
  ############## select the best solution (maximum posteria solution) for given K
  ##########################################################
  tmpK <- unlist(res.WES[['K']])
  unique.K <- sort(unique(tmpK))
  n.K <- length(unique.K)
  MAP <- list()
  for (i in 1:n.K) {
    tmpX <- res.WES[['evid']]
    tmpX[tmpK != unique.K[i]] <- 0
    MAP[[i]] <- which.min(tmpX)
  }
  
  MAP <- unlist(MAP)
  names(MAP) <- unique.K
  MAP.nontrivial <- MAP[names(MAP)!=1]
  ##########################################################
  ##########################################################
  
  n.K <- length(MAP.nontrivial)
  signatures <- list()
  if (n.K > 0) { 
    for (j in 1:n.K) {
      this <- list()
      res <- results[[j]]
      W <- res[[1]]
      H <- res[[2]]
      W1 <- W
      H1 <- H
      W.norm <- apply(W,2,function(x) x/sum(x))
      for (i in 1:ncol(W)) {
        H1[i,] <- H1[i,]*colSums(W)[i]
        W1[,i] <- W1[,i]*rowSums(H)[i]
      }
      lambda <- res[[5]]
      df <- get.df.solution(W1,H,lambda,tumor.type)
      K <- length(table(df$class))
      
      ############# Signature plot
      p <- plot.lego.observed.barplot(df,paste("Mutation Signatures in ",tumor.type,sep=""))
      this['plot'] <- p
      plot(p)

      index.nonzero <- colSums(W) != 0
      lambda <- unlist(res[[5]][length(res[[5]])])
      lambda <- lambda/min(lambda)
      lambda <- lambda[index.nonzero]
      names(lambda) <- paste("W",seq(1:length(lambda)),sep="")
      K <- sum(index.nonzero)
      if (K != 1) {
        W0.tumor <- W[,index.nonzero]
        H0.tumor <- H[index.nonzero,]
        x <- W0.tumor
        y <- H0.tumor
        x.norm <- apply(x,2,function(x) x/sum(x))
        W1.tumor <- x.norm
        for (j in 1:K) y[j,] <- y[j,]*colSums(x)[j]
        H1.tumor <- y
        y.norm <- apply(y,2,function(x) x/sum(x))
        H2.tumor <- y.norm
      } else {
        stop("No non-trivial solutations; All simulations converged to a trivial solution with one signature")
      }
      
      ############# Reconstructing the activity of signatures
      W.mid <- W1.tumor
      H.mid <- H1.tumor
      H.norm <- H2.tumor
      if (length(grep("__",colnames(H.mid))) != 0) {
        hyper <- colnames(H.mid)[grep("__",colnames(H.mid))]
        H.hyper <- H.mid[,colnames(H.mid) %in% hyper]
        H.nonhyper <- H.mid[,!(colnames(H.mid) %in% hyper)]
        sample.hyper <- colnames(H.hyper)
        sample.hyper <- sapply(sample.hyper,function(x) strsplit(x,"__")[[1]][[1]])
        unique.hyper <- unique(sample.hyper)
        n.hyper <- length(unique.hyper)
        x.hyper <- array(0,dim=c(nrow(H.hyper),n.hyper))
        for (i in 1:n.hyper) {
          x.hyper[,i] <- rowSums(H.hyper[,sample.hyper %in% unique.hyper[i]])
        }
        colnames(x.hyper) <- unique.hyper
        rownames(x.hyper) <- rownames(H.mid)
        H.mid <- (cbind(H.nonhyper,x.hyper))
        H.norm <- apply(H.mid,2,function(x) x/sum(x))
      }
      W.norm <- apply(W.mid,2,function(x) x/sum(x))
      
      ##########################################################
      ############# W.norm = extracted signatures normalized to one
      ############# H.mid = activity of signatures across samples (expected mutations associated with signatures)
      ############# H.norm = normalized signature activity 
      ##########################################################
      WH <- list('W.norm' = W.norm, 'H.mid' = H.mid, 'H.norm' = H.norm)
      this['signatures'] = W.norm
      this['signature_mutations'] = H.mid
      this['normalized_mutations'] = H.norm

      ############# Activity plot
      p1 <- plot.activity.barplot(H.mid,H.norm,1.0,tumor.type)
      plot(p1)
      this['plot_activity'] <- p1
      this['K'] <- K
      signatures[i] <- this
    }
  }
  
  return(signatures)
}

summarize_bayesNMF_result <- function(result) {
  W <- result[[1]]
  H <- result[[2]]
  evid <- result[[4]]
  evid <- evid[[length(evid)]]
  lambda <- result[[5]]
  lambda <- lambda[[length(lambda)]]
  lambda.min <- min(lambda)
  lambda <- lambda-lambda.min
  lambda.norm <- lambda/sum(lambda)
  index <- lambda.norm > 0.01
  K <- sum(index)
  return(list('W' = W, 'H' = H, 'lambda' = lambda, 'K' = K, 'evid' = evid))
}

summarize_results <- function(result_list) {
  W.ALL <- list()
  H.ALL <- list()
  K.ALL <- list()
  lambda.ALL <- list()
  evid.ALL <- list()
  n.iter <- length(result_list)
  for (i in 1:n.iter) {
    this_res <- result_list[[i]]
    W <- this_res[[1]]
    H <- this_res[[2]]
    evid <- this_res[[4]]
    evid <- evid[[length(evid)]]
    lambda <- this_res[[5]]
    lambda <- lambda[[length(lambda)]]
    lambda.min <- min(lambda)
    lambda <- lambda-lambda.min
    lambda.norm <- lambda/sum(lambda)
    index <- lambda.norm > 0.01
    K <- sum(index)
    W.ALL[[i]] <- W
    H.ALL[[i]] <- H
    lambda.ALL[[i]] <- lambda
    K.ALL[[i]] <- K
    evid.ALL[[i]] <- evid
  }
  return(list('W' = W.ALL,'H' = H.ALL,'lambda' = lambda.ALL,'K' = K.ALL, 'evid' = evid.ALL))
}
