############################################
############################################
######## Estimating Signaling Games ########
########   Monte Carlo Experiment   ########
########     Unique Equilibrium     ########
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



## Setup the parameters of the experiment
nPerGame <- c(5, 25, 50, 100, 200)
nGames <- c(25, 50, 100, 200)
B <- 1000
nworkers <- detectCores()-1

Bparams <- expand.grid(nGames,nPerGame)
Results <- list()

workers <- makeCluster(nworkers) #in parallel
registerDoParallel(workers)


# Begin the simulation
set.seed(2) #seed 1 hanged at one point 
for (i in 1:nrow(Bparams)){
  Results[[i]] <- foreach(b=1:B,
                          .packages=c("pbivnorm","rootSolve", "Formula", "randomForest", "mc2d", "maxLik"),
                          .combine=cbind,
                          .multicombine=T,
                          .inorder=F
  ) %dorng% {
    nPerGame <- Bparams[i,2]
    M <- Bparams[i,1]
    
    X <- data.frame(X=runif(M))
    f1 <- # create regression matrices
      as.Formula(~0|#SA
                   1|#VA
                   0|#CB
                   1 | #barWA
                   X| #barWB
                   1| #bara
                   1) #VB
    
    truth1 <- c(1, -1.9, -2.9,.1, -1.2,1)
    truth2 <- c(1, -1.7, -2,.1, -1.2,1)
    truth <- .5*truth1 + .5*truth2
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


	# Yr <- colSums(Y[3:4,])/colSums(Y[2:4,])
	# Yf <- Y[3,]/colSums(Y[3:4,])
	# dat <- cbind(X, Yr, Yf)
	m1 <- randomForest(Yr~X, data = data1, na.action=na.omit, ntree=1000)
	m2 <- randomForest(Yf~X, data = data2, na.action=na.omit, ntree=1000)
	Phat <- list(PRhat = pmin(pmax(predict(m1, newdata=X), 0.0001), .9999),
			     PFhat = pmin(pmax(predict(m2, newdata=X), 0.0001), .9999))

    
    
    fqll <- function(x){
      -QLL.jo(x,Phat$PRhat,Phat$PFhat,Y,regr)
    }
    gr.qll <- function(x){
      -eval_gr_qll(x,Phat$PRhat,Phat$PFhat,Y,regr)
    }
    x0 <- c(runif(6), qlogis(Phat$PRhat))
    
    
    ## Test to make sure that these data will actually work
    while(anyNA(gr.qll(x0))|| any(is.infinite(gr.qll(x0)))){
      X <- data.frame(X=runif(M))
      f1 <- # create regression matrices
        as.Formula(~0|#SA
                     1|#VA
                     0|#CB
                     1 | #barWA
                     X| #barWB
                     1| #bara
                     1) #VB
      
      truth1 <- c(1, -1.9, -2.9,.1, -1.2,1)
      truth2 <- c(1, -1.7, -2,.1, -1.2,1)
      truth <- .5*truth1 + .5*truth2
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


	  # Yr <- colSums(Y[3:4,])/colSums(Y[2:4,])
	  # Yf <- Y[3,]/colSums(Y[3:4,])
	  # dat <- cbind(X, Yr, Yf)
	  m1 <- randomForest(Yr~X, data = data1, na.action=na.omit, ntree=1000)
	  m2 <- randomForest(Yf~X, data = data2, na.action=na.omit, ntree=1000)
	  Phat <- list(PRhat = pmin(pmax(predict(m1, newdata=X), 0.0001), .9999),
		           PFhat = pmin(pmax(predict(m2, newdata=X), 0.0001), .9999))

          
      
      x0 <- c(runif(6),
              qlogis(pmin(pmax(Phat$PRhat,
                               0.0001), .9999)))
    }
    #Once we have a usuable dataset, move to estimation
    
    
    
    #### PL  ####
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
      out.2step$iter <- out.2step$iter
    } 
    
    #### tML ####
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
    
    #### NPL ####
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
  save(Results, Bparams, file="MonteCarloResults_Unique.rdata")  
}
stopCluster(workers)


