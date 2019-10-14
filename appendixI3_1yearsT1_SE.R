##########################################################
## Sanctions: Combine Models and Create Standard Errors ##
##########################################################

rm(list=ls())

suppressMessages(library(pbivnorm))
suppressMessages(library(Formula))
suppressMessages(library(foreign))
suppressMessages(library(randomForest))
suppressMessages(library(data.table))
suppressMessages(library(maxLik))
suppressMessages(library(rootSolve))
suppressMessages(library(doParallel))
suppressMessages(library(doRNG))

source("signalingFunctions_main.r")
source("gradientFunctions.r")
load("SanctionsDataSet_1yearsT1.rdata")

missing <- apply(is.na(Xall), 1, max)
Y <- Yall[,missing==0]
X <- Xall[missing==0,]

y1 <- factor(apply(Y, 2, function(y){ifelse(y[4]>0,4,
                                            ifelse(y[3]>0, 3,
                                                   ifelse(y[2]>0, 2,
                                                          1)))}))
Y <- model.matrix(~ y1 +0)
Y <- t(Y)
colnames(Y) <- colnames(Yall[,missing==0])
M <- ncol(Y)
# create regression matrices
f1 <- as.Formula(~ sqrt(senderecondep) + senderdemocracy + contig + ally -1|#SA
                   1|#VA
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

# compute PRhat for starting values
index1 <- colSums(Y[2:4,]) >= 1
polyX1 <- cbind((colSums(Y[3:4,])/colSums(Y[2:4,]))[index1], X[index1,])
index2 <- colSums(Y[3:4,]) >= 1
polyX2 <- cbind(((Y[3,])/colSums(Y[3:4,]))[index2], X[index2,])
colnames(polyX1)[1] <- "V1"
colnames(polyX2)[1] <- "V1"
polyX1[,`:=`(senderecondep=sqrt(senderecondep),
             targetecondep=sqrt(targetecondep))]

polyX2[,`:=`(senderecondep=sqrt(senderecondep),
             targetecondep=sqrt(targetecondep))]

####bootstrap the standard errors####
set.seed(12345)
B <- 50
PL.boot <- matrix(0, ncol=19, nrow=B)
NPL.boot <- matrix(0, ncol=19, nrow=B)
for(i in 1:B){
  
  use <- sample(1:M, replace=T)
  Y.boot <- Y[,use]
  X.boot <- X[use,]
  index1.boot <- colSums(Y.boot[2:4,]) >= 1
  polyX1.boot <- cbind((colSums(Y.boot[3:4,])/colSums(Y.boot[2:4,]))[index1.boot],
                       X.boot[index1.boot,])
  index2.boot <- colSums(Y.boot[3:4,]) >= 1
  polyX2.boot <- cbind(((Y.boot[3,])/colSums(Y.boot[3:4,]))[index2.boot], X.boot[index2.boot,])
  colnames(polyX1.boot)[1] <- "V1"
  colnames(polyX2.boot)[1] <- "V1"
  polyX1.boot[,`:=`(senderecondep=sqrt(senderecondep),
                    targetecondep=sqrt(targetecondep))]
  polyX2.boot[,`:=`(senderecondep=sqrt(senderecondep),
                    targetecondep=sqrt(targetecondep))]
  m1 <- randomForest(factor(V1)~ senderecondep + senderdemocracy + contig + ally +
                       targetecondep + anticipatedtargetcosts +
                       lncaprat + targetdemocracy,
                     data = polyX1.boot,
                     ntree=50, mtry=6
  )
  
  m2 <- randomForest(factor(V1)~ (senderecondep) + senderdemocracy + contig + ally +
                       (targetecondep) + anticipatedtargetcosts + lncaprat + targetdemocracy,
                     data = polyX2.boot,
                     ntree=1000, mtry=8)
  PRhat.boot <- predict(m1, newdata=X.boot, type = "prob")[,2]
  PFhat.boot <- predict(m2, newdata=X.boot, type = "prob")[,2]
  
  
  regr.boot <- list()                   
  for(j in 1:7){
    regr.boot[[j]] <- model.matrix(f1, data=X.boot, rhs=j)
  }
  names(regr.boot) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")  
  
  # estimate 2step
  fqll <- function(x){-QLL.jo(x,PRhat.boot,PFhat.boot,Y.boot,regr.boot)}
  grqll <- function(x){-eval_gr_qll(x, PRhat.boot, PFhat.boot,Y.boot,regr.boot)}
  
  #gradient methods
  x0 <- out$est
  out.pl.boot <- try(maxLik::maxLik(start=x0, fqll,
                                    gr=grqll,
                                    method='NR',
                                    control=list(tol=1e-10,
                                                 reltol=1e-10,
                                                 gradtol=1e-10,
                                                 iterlim=5000)))
  if(class(out.pl.boot[[1]])=="character" || (!out.pl.boot$code %in% 1:2)){
    PL.boot[i,] <- rep(NA, 19)
    NPL.boot[i,] <- rep(NA, 19)
  }else{
    PL.boot[i,] <- out.pl.boot$est
    
    
    
    
    out.NPL.boot <- out.pl.boot
    Phat <- list(PRhat=PRhat.boot, PFhat=PFhat.boot)
    
    eval = 1000; tol = 1e-5; maxit=25
    iter <- 0
    Phat0 <- Phat
    while(eval > tol & iter < maxit){
      
      Uk <- vec2U.regr(out.NPL.boot$estimate, regr.boot)
      Pk.F <- eqProbs(Phat$PRhat, Uk, RemoveZeros = T)[,3]
      Pk.R <- pnorm((Phat$PFhat*Uk$barWB + (1-Phat$PFhat)*Uk$VB - Uk$CB)/Phat$PFhat)
      Phat.k_1 <- Phat
      Phat <- list(PRhat = Pk.R, PFhat = Pk.F)
      Phat$PRhat <-  pmin(pmax(Phat$PRhat,
                               0.0001), .9999)
      Phat$PFhat <-  pmin(pmax(Phat$PFhat,
                               0.0001), .9999)
      fqll <- function(x){- QLL.jo(x,Phat$PRhat,Phat$PFhat,Y.boot,regr.boot)}
      grqll <- function(x){- eval_gr_qll(x, Phat$PRhat, Phat$PFhat,Y.boot,regr.boot)}
      out.NPL.k <- try(maxLik(start=out.NPL.boot$est, logLik=fqll, grad=grqll, method="NR"))
      if(class(out.NPL.k[[1]])=="character" || out.NPL.k$code==100){
        out.NPL.boot <- out.NPL.k
        break
      }
      out.NPL.k$convergence <- out.NPL.k$code
      eval <- max( abs(c(out.NPL.k$est, unlist(Phat.k_1)) -c(out.NPL.boot$est, unlist(Phat))))
      out.NPL.boot <- out.NPL.k
      iter <- iter + 1 
    }
    if(class(out.NPL.boot[[1]])=="character" || (!out.NPL.boot$code %in% 1:2)){
      NPL.boot[i,] <- rep(NA, 19)
    }else{
      NPL.boot[i,] <- out.NPL.boot$est
    } 
  } 
}

save(list=c("NPL.boot", "PL.boot"), file="appendixI3_bootstraps_1yearsT1.rdata")
load("appendixI3_bootstraps_1yearsT1.rdata")
######################### PICK IT UP HERE ###############################################
varnames <- paste(names(unlist(lapply(regr, colnames))), ": ", unlist(lapply(regr, colnames)), sep="")

load("appendixI3_output_1yearsT1.Rdata")
####PL####
out <- PL.out$model
vcov.PL <- var(PL.boot, na.rm=TRUE)
se.PL <- sqrt(diag(vcov.PL))
PL.hat <- out$estimate
z <- PL.hat/se.PL
p <- 2*pnorm(abs(z),lower.tail = F)<0.1
names(PL.hat) <- varnames
out.PL <- round(cbind(sign(PL.hat), PL.hat, se.PL, z, p),2)
colnames(out.PL) <- c("Sign", "Est.", "S.E.", "z-stat", "p<0.05")
cat("PL:\n") #The PL column ofTable 8
print(out.PL[19:20,])
cat(paste("PL log Likelihood:", round(out$max,2), "\n"))

####NPL####
out.NPL <- NPL.out$model

vcov.NPL <- var(NPL.boot, na.rm=TRUE)
se.NPL <- sqrt(diag(vcov.NPL))
NPL.hat <- out.NPL$estimate
z <- NPL.hat/se.NPL
p <- 2*pnorm(abs(z),lower.tail = F)<0.1
names(NPL.hat) <- varnames
out.npl <- round(cbind(sign(NPL.hat), NPL.hat, se.NPL, z, p),2)
colnames(out.npl) <- c("Sign", "Est.", "S.E.", "z-stat", "p<0.05")
cat("NPL:\n") # NPL column ofTable 8
print(out.npl[19:20,])
cat(paste("NPL log Likelihood:", round(out.NPL$max,2), "\n"))


