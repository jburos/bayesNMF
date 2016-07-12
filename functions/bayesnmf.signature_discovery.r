############################################################################################
############################################################################################
# Copyright (c) 2016, Broad Institute
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#   
#   Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 
# Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in
# the documentation and/or other materials provided with the
# distribution.
# 
# Neither the name of the Broad Institute nor the names of its
# contributors may be used to endorse or promote products derived
# from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#                                               LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#                                               DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################################
############################################################################################

######################################################################################################
####### Mutation Signature Profiling using Bayesian NMF algorithms 
######################################################################################################
####### For details on the implementation 
####### see J Kim, Mouw K, P Polak et al, Somatic ERCC2 mutations are associated with a distinct genomic signature in urothelial tumors 
####### Nat. Genet. DOI: 10.1038/ng.3557
####### For details on the original algorithms 
####### see Tan, V.Y. & Févotte, C. Automatic relevance determination in nonnegative matrix factorization with the beta-divergence.
####### IEEE Trans. Pattern Anal. Mach. Intell. 35, 1592–1605 (2013).
######################################################################################################

##########################################################
###################### R functions  ######################
##########################################################

## Bayesian NMF algorithm with exponential priors for W and H
BayesNMF.L1KL <- function(V0,n.iter,a0,tol,K,K0,phi) {
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[,active_nodes]
  V <- V0-min(V0) + eps
  Vmin <- min(V)
  Vmax <- max(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  W <- matrix(runif(N * K)*sqrt(Vmax),ncol=K)
  H <- matrix(runif(M * K)*sqrt(Vmax),ncol=M)
  V.ap <- W %*% H + eps
  I <- array(1,dim=c(N,M))
  
  C <- N + M + a0 + 1
  b0 <- sqrt((a0-1)*(a0-2)*mean(V,na.rm=T)/K0)
  lambda.bound <- b0/C
  lambda <- (colSums(W)+rowSums(H)+b0)/C
  lambda.cut <- 1.5*lambda.bound
  
  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  
  iter <- 2
  count <- 1
  while ((del >= tol) & (iter < n.iter)) {
    H <- H * (t(W) %*% (V/V.ap))/(matrix(rep(colSums(W)+phi/lambda,M),ncol=M) + eps)
    V.ap <- W %*% H + eps
    W <- W * ((V/V.ap) %*% t(H))/t(matrix(rep(rowSums(H)+phi/lambda,N),ncol=N) + eps)
    V.ap <- W %*% H + eps
    lambda <- (colSums(W) + rowSums(H) + b0) / C
    del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
    like <- sum(V*log(V/V.ap)+V.ap-V)
    n.like[[iter]] <- like
    n.evid[[iter]] <- like+phi*sum((colSums(W)+rowSums(H)+b0)/lambda+C*log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V-V.ap)^2)
    if (iter %% 100 == 0) {
      cat(iter,n.evid[[iter]],n.like[[iter]],n.error[[iter]],del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
    }
    iter <- iter+1
  }
  return(list(W,H,n.like,n.evid,n.lambda,n.error))
}

## Bayesian NMF algorithm with hlaf-normal priors for W and H
BayesNMF.L2KL <- function(V0,n.iter,a0,tol,K,K0,phi) {
  eps <- 1.e-50
  del <- 1.0
  active_nodes <- colSums(V0) != 0
  V0 <- V0[,active_nodes]
  V <- V0-min(V0)
  Vmax <- mean(V)
  N <- dim(V)[1]
  M <- dim(V)[2]
  
  W <- matrix(runif(N * K)*Vmax,ncol=K)
  H <- matrix(runif(M * K)*Vmax,ncol=M)
  V.ap <- W%*%H+eps
  
  C <- (N+M)/2+a0+1
  b0 <- 3.14*(a0-1)*mean(V)/(2*K0)
  lambda <- (0.5*colSums(W*W)+0.5*rowSums(H*H)+b0)/C
  lambda.bound <- b0/C
  lambda.cut <- b0/C*1.25
  
  I <- array(1,dim=c(N,M))
  n.like <- list()
  n.evid <- list()
  n.error <- list()
  n.lambda <- list()
  n.lambda[[1]] <- lambda
  iter <- 2
  count <- 1
  while ((del >= tol) & (iter < n.iter)) {
    H <- H*(t(W)%*%(V/V.ap)/(t(W)%*%I+phi*H/matrix(rep(lambda,M),ncol=M)+eps))^0.5
    V.ap <- W%*%H+eps
    W <- W*((V/V.ap)%*%t(H)/(I%*%t(H)+phi*W/t(matrix(rep(lambda,N),ncol=N))+eps))^0.5
    lambda <- (0.5*colSums(W*W)+0.5*rowSums(H*H)+b0)/C
    V.ap <- W%*%H+eps
    del <- max(abs(lambda-n.lambda[[iter-1]])/n.lambda[[iter-1]])
    like <- sum(V * log((V+eps)/(V.ap+eps)) + V.ap - V)
    n.like[[iter]] <- like
    n.evid[[iter]] <- like + phi*sum((0.5*colSums(W^2)+0.5*rowSums(H^2)+b0)/lambda+C*log(lambda))
    n.lambda[[iter]] <- lambda
    n.error[[iter]] <- sum((V-V.ap)^2)
    if (iter %% 100 == 0) {
      cat(iter,n.evid[[iter]],n.like[[iter]],n.error[[iter]],del,sum(colSums(W)!=0),sum(lambda>=lambda.cut),'\n')
    }
    iter <- iter+1
  }
  return(list(W,H,n.like,n.evid,n.lambda,n.error))
}

######### Handling hypermutant samples; See J Kim et al Nat. Genet. DOI: 10.1038/ng.3557 for details.
get.lego96.hyper <- function(lego96) {
  x <- lego96
  for (i in 1:100) {
    SNV <- colSums(x)
    q1 <- quantile(SNV,prob=1/4)
    q3 <- quantile(SNV,prob=3/4)
    sample.hyper <- colnames(x)[SNV > (median(SNV)+1.5*(q3-q1))]
    if (length(sample.hyper)==0) break
    lego96.hyper <- as.matrix(x[,(colnames(x) %in% sample.hyper)])
    colnames(lego96.hyper) <- sample.hyper
    lego96.nonhyper <- x[,!(colnames(x) %in% sample.hyper)]
    lego96.hyper1 <- apply(lego96.hyper,2,function(x) x/2)
    lego96.hyper2 <- lego96.hyper1
    colnames(lego96.hyper1) <- paste(colnames(lego96.hyper1),1,sep="__")
    colnames(lego96.hyper2) <- paste(colnames(lego96.hyper2),2,sep="__")
    x <- cbind(lego96.nonhyper,lego96.hyper1,lego96.hyper2)
  }
  return(x)
}

######### Visualizing the activity of signatures
plot.activity.barplot <- function(H.mid,H.norm,scale,tumor.type) {
  .theme_ss <- theme_bw(base_size=14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=8*scale, family="mono"),
          axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
          axis.text = element_text(size = 12*scale, family = "mono"))
  ordering <- order(colSums(H.mid),decreasing=T)
  H.mid <- H.mid[,ordering]
  rownames(H.mid) <- paste("W",seq(1:nrow(H.mid)),sep="")
  H.norm <- H.norm[,ordering]
  rownames(H.norm) <- paste("W",seq(1:nrow(H.norm)),sep="")
  sample.ordering <- colnames(H.mid)
  x1 <- melt(H.mid)
  x2 <- melt(H.norm)
  colnames(x1) <- c("Signature","Sample","Activity")
  colnames(x2) <- c("Signature","Sample","Activity")
  x1[,"class0"] <- c("Counts")
  x2[,"class0"] <- c("Fractions")
  df2 <- rbind(x1,x2)
  df2$class0 <- factor(df2$class0,c("Counts","Fractions"))
  df2$Sample <- factor(df2$Sample,sample.ordering)
  scale <- 1
  p = ggplot(df2,aes(x=factor(Sample),y=Activity,fill=factor(Signature)))
  p = p+geom_bar(stat="identity",position='stack',color='black',alpha=0.9)
  p = p + scale_fill_manual(values=c("red","cyan","yellow","blue","magenta","gray50","orange","darkgreen","brown","black",rainbow(10)[4:10]))
  p = p + facet_grid(class0 ~ ., scale = "free_y")
  p = p + ggtitle(paste("Siganture Activities in",tumor.type,sep=" "))
  p = p + theme(plot.title=element_text(lineheight=1.0,face="bold",size=14*scale))
  p = p + xlab("Samples") + ylab("Signature Activities")
  p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
  p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
  p = p + theme(axis.text.x = element_text(angle=90,vjust=0.5,size=8*scale,face="bold",colour="black"))
  p = p + theme(axis.text.y = element_text(size=10*scale,face="bold",colour="black"))
  p = p + theme(legend.title=element_blank())
  p = p + .theme_ss
  p = p + theme(legend.position="top")
  return(p)
}

######### Collecting data from several independent runs
get.stats.simulation <- function(tumor,n.iter,OUTPUT) {
  W.ALL <- list()
  H.ALL <- list()
  K.ALL <- list()
  lambda.ALL <- list()
  evid.ALL <- list()
  for (i in 1:n.iter) {
    x <- load(paste(OUTPUT,paste(method,a0,i,"RData",sep="."),sep=""))
    W <- res[[1]]
    H <- res[[2]]
    evid <- res[[4]]
    evid <- evid[[length(evid)]]
    lambda <- res[[5]]
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
  return(list(W.ALL,H.ALL,lambda.ALL,K.ALL,evid.ALL))
}

######### Visualizing the profile of signatures
plot.lego.observed.barplot <- function(mat,title) {
  scale <- 1.5
  .theme_ss <- theme_bw(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size=10*scale, family="mono"),
          axis.text.y = element_text(hjust = 0.5,size=12*scale, family="mono"),
          axis.text = element_text(size = 16*scale, family = "mono"))
  p = ggplot(mat)
  p = p + geom_bar(aes_string(x="context3",y="signature",fill="base.sub"),stat="identity",position="identity",colour="gray50")
  p = p + facet_grid(class ~ base.sub, scale = "free_y")
  p = p + .theme_ss
  p = p + scale_fill_manual(values=c("cyan","red","yellow","purple","green","blue","black","gray")) #p = p + scale_fill_brewer(palette = "Set1")
  p = p + guides(fill=FALSE) #p = p + theme(legend.position = "none")
  p = p + ggtitle(title)
  p = p + xlab("Motifs") + ylab("Contributions")
  p = p + theme(axis.title.x = element_text(face="bold",colour="black",size=14*scale))
  p = p + theme(axis.title.y = element_text(face="bold",colour="black",size=14*scale))
  return(p)
}

get.df.solution <- function(W,H,lambda,tumor) {
  lambda <- lambda[[length(lambda)]]
  lambda.min <- min(lambda)
  lambda <- (lambda-lambda.min)/max(lambda)
  norm.lambda <- lambda/sum(lambda)
  index <- lambda > 0.01
  x <- colSums(W)/max(colSums(W))
  index <- x > 0.01
  K <- sum(index)
  if (K != 1) {
    W <- W[,index]
    H <- H[index,]
    norm.W <- t(apply(W,1,function(x) x/sum(x)))
    norm.H <- apply(H,2,function(x) x/sum(x))
  }
  colnames(W) <- paste("W",seq(1:ncol(W)),sep="")
  rownames(H) <- paste("H",seq(1:ncol(W)),sep="")
  
  context96 <- rownames(W)[1:96]
  context96 <- sub("->","",context96)
  context96 <- gsub("[.]","",context96)
  index96 <- c(order(substring(context96,1,2),decreasing=F))
  W.SNP <- as.vector(W[index96,])
  seq <- c(context96[index96])
  context4 <- rep(seq,K)
  context3 <- paste(substring(context4,3,3),"-",substring(context4,4,4),sep="")
  base.sub <- paste(substring(context4,1,1),"->",substring(context4,2,2),sep="")
  signature.class <- as.vector(t(matrix(t(rep(1:K,96)),nrow=K)))
  signature.class <- paste("W",signature.class,sep="")
  df <- data.frame(W.SNP,context4,context3,base.sub,signature.class,rep(tumor,nrow(W)*ncol(W)))
  colnames(df) <- c("signature","context4","context3","base.sub","class","tumor")
  return(df)
}
##########################################################
##########################################################
##########################################################
