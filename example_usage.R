
library(gplots)
library(ggplot2)
library(reshape2)

#CURRENT <- paste(getwd(),'/',sep="")
#SOURCE <- paste(CURRENT,"SOURCE/",sep="")
#PARAMETER <- SOURCE
#source(paste(SOURCE,"get.context96_annotated.from.maf.R",sep=""))
#source(paste(SOURCE,"get.util.signature.R",sep=""))
#source(paste(SOURCE,"plot.signature.tools.Alexsandrov.R",sep=""))

CURRENT <- paste(getwd(),"/",sep="")
OUTPUT <- paste(CURRENT,"OUTPUT_lego96/",sep="")
system(paste("mkdir",OUTPUT,sep=" "))

############## INPUT ######################################
## lego96 - mutation counts matrix (96 by # of samples) 
## This matrix should contain mutation counts along 96 tri-nucleotide mutation contexts (rows) across samples (columns). 
## Rownames of the lego matrix should be 4-letters ex) CGTA (C to G mutation at 5'-TCA-3'contexts) (see the acoompanied example lego matrix).
###########################################################

#####################################################
############## BayesNMF parameters
############## n.iter = number of independent simulations
############## Kcol = number of initial signatures
############## tol = tolerance for convergence
############## a0 = hyper-parameter
n.iter <- 10 
Kcol <- 96   
tol <- 1.e-07
a0 <- 10
tumor.type <- "tumor.type" ### please specify your cohort name here
##################################

##################################
############### Choose pirors for W and H
############### Default = L1KL (expoential priors); L2KL (half-normal priors)
##################################
prior <- "L1KL" 
if (prior=="L1KL") {
  method <- paste("L1KL.lego96",tumor.type,sep=".")
} else {
  method <- paste("L2KL.lego96",tumor.type,sep=".")
}
##################################

##################################
############## Default = FALSE ; TRUE - to reduce the effect of hyper-mutant samples in the signature discovery
##################################
hyper <- FALSE
if (hyper) {
  lego96 <- get.lego96.hyper(lego96)
  method <- paste(method,"hyper",sep=".")
}
##################################

##########################################################
###################### Running the algorithms ############
##########################################################
for (i in 1:n.iter) {
  if (prior=="L1KL") {
    res <- BayesNMF.L1KL(as.matrix(lego96),100000,a0,tol,Kcol,Kcol,1.0)
  } else {
    res <- BayesNMF.L2KL(as.matrix(lego96),100000,a0,tol,Kcol,Kcol,1.0)
  }
  save(res,file=paste(OUTPUT,paste(method,a0,i,"RData",sep="."),sep=""))
}

##########################################################
###################### Analysis ##########################
##########################################################
res.WES <- get.stats.simulation(tumor.type,n.iter,OUTPUT)

############## frequency figure 
pdf(file=paste(OUTPUT,paste(method,a0,"signature.freq.pdf",sep="."),sep=""),width=4,height=5)
s1 <- 1.5
s2 <- 2.0
par(mfrow=c(1,1))
par(mar=c(5,5,2,1))
barplot(table(unlist(res.WES[[4]])),cex=s1,cex.axis=s1,cex.main=s1,cex.names=s1,cex.lab=s1,xlab="# of signatures",ylab="Freq.",main=paste(tumor.type,sep="."))
dev.off()

##########################################################
############## select the best solution (maximum posteria solution) for given K
##########################################################
tmpK <- unlist(res.WES[[4]])
unique.K <- sort(unique(tmpK))
n.K <- length(unique.K)
MAP <- list()
for (i in 1:n.K) {
  tmpX <- res.WES[[5]]
  tmpX[tmpK != unique.K[i]] <- 0
  MAP[[i]] <- which.min(tmpX)
}

tmpK <- unlist(res.WES[[4]])
unique.K <- sort(unique(tmpK))
n.K <- length(unique.K)
MAP <- list()
for (i in 1:n.K) {
  tmpX <- res.WES[[5]]
  tmpX[tmpK != unique.K[i]] <- 0
  MAP[[i]] <- which.min(tmpX)
}
MAP <- unlist(MAP)
names(MAP) <- unique.K
MAP.nontrivial <- MAP[names(MAP)!=1]
##########################################################
##########################################################

n.K <- length(MAP.nontrivial)
if (n.K > 0) { 
  for (j in 1:n.K) {
    load(paste(OUTPUT,paste(method,a0,MAP.nontrivial[[j]],"RData",sep="."),sep=""))
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
    width <- 16
    height <- ifelse(K==1,3,K*2)
    pdf(file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"signature.pdf",sep="."),sep=""),width=width,height=height)
    p <- plot.lego.observed.barplot(df,paste("Mutation Signatures in ",tumor.type,sep=""))
    plot(p)
    dev.off()
    
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
    WH <- list(W.norm,H.mid,H.norm)
    save(WH,file=paste(OUTPUT,paste(method,a0,paste("MAP",K,sep=""),"WH.RData",sep="."),sep=""))
    
    ############# Activity plot
    p1 <- plot.activity.barplot(H.mid,H.norm,1.0,tumor.type)
    pdf(file = paste(OUTPUT,paste(method,a0,"activity.barplot1",K,"pdf",sep="."),sep=""),width=15,height=12)
    plot(p1)
    dev.off()
  }
}
