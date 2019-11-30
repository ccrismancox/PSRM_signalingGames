rm(list=ls())

suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
suppressMessages(library(maxLik))
suppressMessages(library(foreach))
suppressMessages(library(rootSolve))
suppressMessages(library(pbivnorm))
suppressMessages(library(Formula))
suppressMessages(library(mc2d))
suppressMessages(library(randomForest))


source("signalingFunctions_main.r")
source("gradientFunctions.r")

####data####
load("SanctionsDataSet.rdata")


# choose 1 year, 5 year, or 10 year (1,2,3, respectively)
year <- 3

# remove irrelevant and missing dyads values
missing <- apply(is.na(Xall), 1, max)

Y <- Yall[,missing==0]
X <- Xall[missing==0,]


# create regression matrices
f1 <- as.Formula(~ sqrt(senderecondep) + senderdemocracy + contig + ally -1|#SA
                   anticipatedsendercosts|#VA
                   sqrt(targetecondep) + anticipatedtargetcosts + contig + ally|#CB
                   sqrt(senderecondep) + senderdemocracy + lncaprat | #barWA
                   targetdemocracy + lncaprat| #barWB
                   senderdemocracy| #bara
                   -1) #VB
regr <- list()
for(i in 1:7){
  regr[[i]] <- model.matrix(f1, data=X, rhs=i)
}
names(regr) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")



####CMLE####
load('CMLE_estimation_output.rdata')
results <- c(output)
LL.cmle <- -value
truth <- results[1:20]


## Setup the parameters of the experiment
nPerGame <- 120
nGames <- 418
B <- 1000
nworkers <- detectCores()-1



## Initialize cluster
workers <- makeCluster(nworkers)


####Actual Cluster####
registerDoParallel(workers)



set.seed(1)
Results <- foreach(b=1:B,
                        .packages=c("pbivnorm",
                                    "maxLik",
                                    "rootSolve",
                                    "Matrix",
                                    "randomForest",
                                    "Formula",
                                    "mc2d"),
                        .combine=cbind,
                        .multicombine=T,
                        .inorder=F
) %dorng% {

  U <- vec2U.regr(truth, regr)
  Pstar <- plogis(results[-c(1:20)])
  Y <- genData.jo(nPerGame, Pstar, U)
  
  X <- do.call(cbind, regr)
  X1 <- X2 <-  unique(X, MARGIN=2)
 
  X1 <- X1[,which(colnames(X1) != "(Intercept)")]
  X2 <- X2[,which(colnames(X2) != "(Intercept)")]
  
  index1 <- colSums(Y[2:4,]) >= 1
  index2 <- colSums(Y[3:4,]) >= 1
  
  phat.X <- list(X1 = cbind((colSums(Y[3:4,])/colSums(Y[2:4,]))[index1], X1[index1,]),
                 X2 = cbind(((Y[3,])/colSums(Y[3:4,]))[index2], X2[index2,]))
  phat.X <- lapply(phat.X, as.data.frame)
 
  m1 <- randomForest(y=phat.X$X1[,1], x=phat.X$X1[,-1])
  m2 <- randomForest(y=phat.X$X2[,1], x=phat.X$X2[,-1])
  Phat <- list(PRhat = predict(m1, newdata=X1, type="response"),
               PFhat = predict(m2, newdata=X2, type="response"))
  
  ## estimate two step
  fqll <- function(x){
    -QLL.jo(x,Phat$PRhat,Phat$PFhat,Y,regr)
  }
  gr.qll <- function(x){
    -eval_gr_qll(x,Phat$PRhat,Phat$PFhat,Y,regr)
  }
  x0 <- c(rnorm(20)*.05,qlogis(Phat$PRhat))
  
 
  ptm <- proc.time()[3]
  out.2step <- try(maxLik(start=x0[1:20], logLik=fqll, grad=gr.qll, method="NR"))
  out.2step$time <- proc.time()[3] - ptm
  
  if(class(out.2step[[1]])=="character"){
    out.2step$par <- rep(NA, 20)
    out.2step$convergence <- -99
    out.2step$iter <- -99
  }else{
    out.2step$par <- out.2step$est
    out.2step$convergence <- out.2step$code
    out.2step$iter <- out.2step$counts[2]
  } 
  
 
  
  ## estimate nested-fixpoint
  fnfxp <- function(x){LL.nfxp(x,Y,regr)}
  
  ptm <- proc.time()[3]
  out.nfxp  <- try(optim(fn=fnfxp, par=x0[1:20], method='Nelder-Mead',
                         control=list(maxit=15000)))
  out.nfxp$time <- proc.time()[3] - ptm
  out.nfxp$iter <- out.nfxp$counts[1]
  if(class(out.nfxp[[1]])=="character"){
    out.nfxp$par <- rep(NA, 20)
    out.nfxp$convergence <- -99
    out.nfxp$iter <- -99
  }

  ptm <- proc.time()[3]
  out.NPL <- out.2step 
  if(class(out.NPL[[1]])!="character"){
    eval = 1000; tol = 1e-6; maxit=100
    iter <- 0
    Phat0 <- Phat
    while(eval > tol & iter < maxit){
      
      Uk <- vec2U.regr(out.NPL$estimate, regr)
      Pk.F <- eqProbs(Phat$PRhat, Uk, RemoveZeros = T)[,3]
      Pk.R <- pnorm((Phat$PFhat*Uk$barWB + (1-Phat$PFhat)*Uk$VB - Uk$CB)/Phat$PFhat)
      Phat.k_1 <- Phat
      Phat <- list(PRhat = Pk.R, PFhat = Pk.F)
      Phat$PRhat <-  pmin(pmax(Phat$PRhat,
                               0.0001), .9999)
      Phat$PFhat <-  pmin(pmax(Phat$PFhat,
                               0.0001), .9999)
      
      out.NPL.k <- try(maxLik(start=out.NPL$par, logLik=fqll, grad=gr.qll, method="NR"))
      if(class(out.NPL.k[[1]])=="character" || out.NPL.k$code==100){
        out.NPL <- out.NPL.k
        break
      }
      out.NPL.k$par <- out.NPL.k$est
      out.NPL.k$convergence <- out.NPL.k$code
      eval <- max( abs(out.NPL.k$par -out.NPL$par))
      print(eval)
      out.NPL <- out.NPL.k
      iter <- iter + 1 
    }
    out.NPL$time <- proc.time()[3] - ptm
    out.NPL$iter <- iter
    
    if(class(out.NPL[[1]])=="character"|| out.NPL.k$code==100){
      out.NPL$par <- rep(NA, 20)
      out.NPL$convergence <- -99
      out.NPL$iter <- -99
    }else{
      out.NPL$convergence <- ifelse(iter==maxit,
                                    -69,
                                    out.NPL$convergence)
      out.NPL$convergence <- ifelse(eval==0,
                                    -99,
                                    out.NPL$convergence)
    }
  }else{
    out.NPL$par <- rep(NA, 20)
    out.NPL$convergence <- -99
    out.NPL$iter <- -99
  }  
  
  
  c(out.2step$par, out.2step$convergence, out.2step$time, out.2step$iter,
    out.NPL$par, out.NPL$convergence,  out.NPL$time, out.NPL$iter,
    out.nfxp$par, out.nfxp$convergence, out.nfxp$time, out.nfxp$iter
  )

}
save(Results, 
     file = "appendixC4_results.RData")  

stopCluster(workers)

