############################################
############################################
######## Estimating Signaling Games ########
########   Monte Carlo Experiment   ########
############################################
############################################


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


B <- 1000
pUnstable <- seq(0, 1, length=11)
Results <- list()
nworkers <- detectCores()-1

workers <- makeCluster(nworkers)
registerDoParallel(workers)


set.seed(1)
for (i in 1:11){
  Results[[i]] <- foreach(b=1:B,
                          .packages=c("pbivnorm","rootSolve","Matrix", "Formula", "randomForest", "mc2d", "maxLik", "numDeriv"),
                          .combine=cbind,
                          .multicombine=T,
                          .inorder=F
  ) %dorng% {
    nPerGame <- 1000
    M <- 200
    unstable <- M*pUnstable[i]
    unstable <- round(unstable)
    
    
    X <- data.frame(X = c(runif(unstable, 1/3,2/3),
                          ifelse(runif(M-unstable)>1/2,
                                 runif(M-unstable,0,1/3),
                                 runif(M-unstable,2/3,1))))
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
    
    
    # compute PRhat for starting values
    index1 <- colSums(Y[2:4,]) >= 1 #observations where B chooses resist or not
    index2 <- colSums(Y[3:4,]) >= 1 #observations where A chooses SF or BD
    data1 <- cbind.data.frame((colSums(Y[3:4,])/colSums(Y[2:4,]))[index1], X[index1,])  #data where B chooses resist or not
    data2 <- cbind.data.frame(((Y[3,])/colSums(Y[3:4,]))[index2], X[index2,])  #data where A chooses SF or BD
    colnames(data1) <- c("Yr", "X")
    colnames(data2) <- c("Yf", "X")
    

    m1 <- randomForest(Yr~X, data = data1, na.action=na.omit, ntree=1000)
    m2 <- randomForest(Yf~X, data = data2, na.action=na.omit, ntree=1000)
    Phat <- list(PRhat = pmin(pmax(predict(m1, newdata=X), 0.0001), .9999),
                 PFhat = pmin(pmax(predict(m2, newdata=X), 0.0001), .9999))
    
    
    
    ## estimate two step
    fqll <- function(x){
      -QLL.jo(x,Phat$PRhat,Phat$PFhat,Y,regr)
    }
    gr.qll <- function(x){
      -eval_gr_qll(x,Phat$PRhat,Phat$PFhat,Y,regr)
    }
    x0 <- c(runif(6), qlogis(Phat$PRhat))
    
    ## Test to make sure that this data will actually work
    ## Minor issue in smaller samples
    while(anyNA(gr.qll(x0))|| any(is.infinite(gr.qll(x0)))){
      X <- data.frame(X = c(runif(unstable, 1/3,2/3),
                            ifelse(runif(M-unstable)>1/2,
                                   runif(M-unstable,0,1/3),
                                   runif(M-unstable,2/3,1))))
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
      
      
      
      # compute PRhat for starting values
      index1 <- colSums(Y[2:4,]) >= 1 #observations where B chooses resist or not
      index2 <- colSums(Y[3:4,]) >= 1 #observations where A chooses SF or BD
      data1 <- cbind.data.frame((colSums(Y[3:4,])/colSums(Y[2:4,]))[index1], X[index1,])  #data where B chooses resist or not
      data2 <- cbind.data.frame(((Y[3,])/colSums(Y[3:4,]))[index2], X[index2,])  #data where A chooses SF or BD
      colnames(data1) <- c("Yr", "X")
      colnames(data2) <- c("Yf", "X")
      
      m1 <- randomForest(Yr~X, data = data1, na.action=na.omit, ntree=1000)
      m2 <- randomForest(Yf~X, data = data2, na.action=na.omit, ntree=1000)
      Phat <- list(PRhat = pmin(pmax(predict(m1, newdata=X), 0.0001), .9999),
                   PFhat = pmin(pmax(predict(m2, newdata=X), 0.0001), .9999))
      
      
      x0 <- c(runif(6), qlogis(Phat$PRhat))
    }

    #PL 
    ptm <- proc.time()[3]
    out.2step <- try(maxLik(start=x0[1:6], logLik=fqll, grad=gr.qll, method="NR"))
    out.2step$time <- proc.time()[3] - ptm
    
    if(class(out.2step[[1]])=="character"){
      out.2step$par <- rep(NA, 6)
      out.2step$convergence <- -99
      out.2step$iter <- -99
    }else{
      out.2step$par <- out.2step$est
      out.2step$convergence <- out.2step$code
      out.2step$iter <- out.2step$counts[2]
    } 
    
    # estimate nested-fixpoint
    fnfxp <- function(x){LL.nfxp(x,Y,regr)}
    
    ptm <- proc.time()[3]
    out.nfxp  <- try(optim(fn=fnfxp, par=x0[1:6], method='Nelder-Mead',
                           control=list(maxit=5000)))
    out.nfxp$time <- proc.time()[3] - ptm
    out.nfxp$iter <- out.nfxp$counts[1]
    
    if(class(out.nfxp[[1]])=="character"){
      out.nfxp$par <- rep(NA, 6)
      out.nfxp$convergence <- -99
      out.nfxp$iter <- -99
    }

    ptm <- proc.time()[3]
    out.NPL <- out.2step 
    if(class(out.NPL[[1]])!="character"){
      eval = 1000; tol = 1e-5; maxit=500
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
        eval <- mean((c(out.NPL.k$par, unlist(Phat)) -c(out.NPL$par,unlist(Phat.k_1)))^2)
        out.NPL <- out.NPL.k
        iter <- iter + 1
      }
      out.NPL$time <- proc.time()[3] - ptm
      out.NPL$iter <- iter
      
      if(class(out.NPL[[1]])=="character"|| out.NPL.k$code==100){
        out.NPL$par <- rep(NA, 6)
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
      out.NPL$par <- rep(NA, 6)
      out.NPL$convergence <- -99
      out.NPL$iter <- -99
    }
    
 
    c(out.2step$par, out.2step$convergence, out.2step$time, out.2step$iter,
      out.NPL$par, out.NPL$convergence,  out.NPL$time, out.NPL$iter,
      out.nfxp$par, out.nfxp$convergence, out.nfxp$time, out.nfxp$iter
    )

  }
  save(list="Results", file="MonteCarloUnstable.rdata")  
}
stopCluster(workers)




