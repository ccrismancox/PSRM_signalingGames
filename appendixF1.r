rm(list=ls())

suppressMessages(library(rootSolve))
suppressMessages(library(pbivnorm))
suppressMessages(library(Formula))
suppressMessages(library(mc2d))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
suppressMessages(library(foreach))
suppressMessages(library(ggplot2))
suppressMessages(library(matrixStats))

source("signalingFunctions_main.r")
source("gradientFunctions.r")


#### Edit selectEq for M=1 #### 
selectEq <- function(X){
  M <- length(X)
  # select equilibria
  Pstar <- c(0.2964518, 0.4715766, 0.8740314) # these are the only eq
  Pstar <- rowSums(matrix(rep(Pstar, each=M),nrow=M) * cbind(X<1/3,X<2/3 & X>1/3, X>2/3))
  
  # compute equilibra 
  f <- function(p){const.jo(p,U)}
  grf <- function(p){diag(1-eval_gr_fh(Pstar, U), nrow=length(1-eval_gr_fh(Pstar, U)))}
  out <- multiroot(f,Pstar, jacfunc=grf, jactype="fullusr", ctol=1e-10,rtol=1e-10,atol=1e-10)
  Pstar <- out$root	
  
  return(Pstar)
}

set.seed(1)
nPerGame <- 200
M <- 1
X <- data.frame(X=runif(M))

f1 <- # create regression matrices
  as.Formula(~0|#SA
               1|#VA
               0|#CB
               1 | #barWA
               X| #barWB
               1| #bara
               1) #VB

truth <- c(1, -1.9, -2.9,.1, -1.2,1)
regr <- list()                   
for(j in 1:7){
  regr[[j]] <- model.matrix(f1, data=X, rhs=j)
}
names(regr) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")  

U <- vec2U.regr(truth, regr)
Pstar <- selectEq(X$X)
Y <- genData.jo(nPerGame, Pstar, U)






LL.nfxp2 <- function(x, Y,regr, control=list(ftol=1e-6,maxit=500), ngrid=10000){
  M <- dim(Y)[2]
  U <-vec2U.regr(x,regr)
  U$sig <- rep(1, M)
  
  root <- matrix(0,M,3)
  eps <- sqrt(.Machine$double.eps)
  grid = seq(from=eps,to=1-eps,length.out=ngrid)
  for (i in 1:M){
    Ui <- as.list(sapply(U, "[", i))
    fgrid <- const.jo(grid,Ui)
    sols <- which(tail(fgrid,-1)*head(fgrid,-1) <= 0)
    sols <- matrix(grid[c(sols, sols+1)], nrow=length(sols))
    sols <- rowMeans(sols)
    
    if(length(sols)==0){
      root[i,] <- rep(sols,3)
    }else{
      root[i,] <- sort(sols)
    }
  }
  
  out <- root
  LL <- matrix(0,M,3)
  for(i in 1:3){
    root <- as.numeric(out[,i])
    EQ <- eqProbs(root,U, RemoveZeros = T)
    OUT <- cbind(1-EQ[,2], EQ[,2]*(1-EQ[,1]), 	
                 EQ[,2]*EQ[,1]*EQ[,3], EQ[,2]*EQ[,1]*(1-EQ[,3]))
    OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
    LL[,i] <- colSums(log(t(OUT))*Y)
  }
  LL <- sum(rowMaxs(LL))
  return(unique(LL))
}





B <- 501
values <- matrix(truth, nrow=B, ncol=6, byrow=T)
values[,3] <- seq(-2.75, -3.05, length=B)


nworkers <- detectCores()-1
workers <- makeCluster(nworkers)
registerDoParallel(workers)


naiveML <- function(x){LL.nfxp2(x,Y=Y,regr=regr)}
LLout1 <-foreach(b=1:B,
                 .packages=c("pbivnorm","rootSolve", "matrixStats"),
                 .inorder=T
) %dorng%{
  xout <- naiveML(values[b,])
}


set.seed(1)
nPerGame <- 200
M <- 10
X <- data.frame(X=runif(M))

f1 <- # create regression matrices
  as.Formula(~0|#SA
               1|#VA
               0|#CB
               1 | #barWA
               X| #barWB
               1| #bara
               1) #VB

truth <- c(1, -1.9, -2.9,.1, -1.2,1)
regr <- list()                   
for(j in 1:7){
  regr[[j]] <- model.matrix(f1, data=X, rhs=j)
}
names(regr) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")  

U <- vec2U.regr(truth, regr)
Pstar <- selectEq(X$X)
Y <- genData.jo(nPerGame, Pstar, U)


LLout2 <-foreach(b=1:B,
                 .packages=c("pbivnorm","rootSolve", "matrixStats"),
                 # .combine = rbind,
                 .inorder=T
) %dorng%{
  xout <- naiveML(values[b,])
}
stopCluster(workers)


barWB <- values[,3]
LLdata <- data.frame(log.likelihood = c(unlist(LLout1), unlist(LLout2)),
                     barWB = barWB,
                     D = rep(c("1 Dyad", paste(M," Dyads",sep="")), each=B))






g2 <-ggplot(LLdata)+
  geom_vline(xintercept=-2.9, size=2, color="orangered1", linetype="dashed")+
  geom_point(aes(y=log.likelihood, x=barWB), color="navyblue", size=2.5)+
  facet_wrap(~D, scales = "free_y", nrow=1)+
  scale_colour_manual(values=c("navyblue", "orangered1"))+
  theme_bw(16)+
  ylab('Log-likelihood')+
  xlab(expression(italic(hat(beta)[bar(W)[B]]^1)))+
  theme(legend.position="bottom",
        legend.title = element_text(size=18, 
                                    face="bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(.65,"in"))


ggsave(g2, file="figure20.pdf", height=6, width=10)



