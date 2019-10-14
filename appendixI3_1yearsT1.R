############################################
############################################
######## Estimating Signaling Games ########
########     Economic Sanctions     ########
########       2step code           ########
############################################
############################################


rm(list=ls())

suppressMessages(library(pbivnorm))
suppressMessages(library(Formula))
suppressMessages(library(foreign))
suppressMessages(library(randomForest))
suppressMessages(library(data.table))
suppressMessages(library(maxLik))
suppressMessages(library(rootSolve))

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

####with the recommended settings####
set.seed(12345)
m1 <- randomForest(factor(V1)~ senderecondep + senderdemocracy + contig + ally + 
                      targetecondep + anticipatedtargetcosts + 
                     lncaprat + targetdemocracy,
                   data = polyX1,
                   ntree=50, mtry=6
)

set.seed(12345)
m2 <- randomForest(factor(V1)~ (senderecondep) + senderdemocracy + contig + ally + 
                     (targetecondep) + anticipatedtargetcosts + lncaprat + targetdemocracy,
                   data = polyX2,
                   ntree=1000, mtry=8
)

PRhat <- predict(m1, newdata=X, type = "prob")[,2]
PFhat <- predict(m2, newdata=X, type = "prob")[,2]

# estimate 2step
fqll <- function(x){-QLL.jo(x,PRhat,PFhat,Y,regr)}
grqll <- function(x){-eval_gr_qll(x, PRhat, PFhat,Y,regr)}

#gradient methods
set.seed(12345)
x0 <- rnorm(19)*.005
out <- maxLik::maxLik(start=x0, fqll,
                      gr=grqll,
                      method='NR',
                      control=list(tol=1e-10,
                                   reltol=1e-10,
                                   gradtol=1e-10,
                                   iterlim=5000))
Phat.PL=list(PRhat=PRhat, PFhat=PFhat)
PL.out <- list(model=out, start=x0, Phat=Phat.PL)



ptm <- proc.time()[3]
out.NPL <- out
Phat <- list(PRhat=PRhat, PFhat=PFhat)

eval = 1000; tol = 1e-7; maxit=100
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
  fqll <- function(x){- QLL.jo(x,Phat$PRhat,Phat$PFhat,Y,regr)}
  grqll <- function(x){- eval_gr_qll(x, Phat$PRhat, Phat$PFhat,Y,regr)}
  out.NPL.k <- try(maxLik(start=out.NPL$est, logLik=fqll, grad=grqll, method="NR"))
  if(class(out.NPL.k[[1]])=="character" || out.NPL.k$code==100){
    out.NPL <- out.NPL.k
    break
  }
  out.NPL.k$convergence <- out.NPL.k$code
  eval <- max( abs(c(out.NPL.k$est, unlist(Phat.k_1)) -c(out.NPL$est, unlist(Phat))))
  cat("NPL updating: ", eval, "Function Value", out.NPL.k$max, "\n")
  out.NPL <- out.NPL.k
  iter <- iter + 1 
}
out.NPL$time <- proc.time()[3] - ptm
out.NPL$iter <- iter

Uk <- vec2U.regr(out.NPL$estimate, regr)
Pk.F <- eqProbs(Phat$PRhat, Uk, RemoveZeros = T)[,3]
Pk.R <- pnorm((Phat$PFhat*Uk$barWB + (1-Phat$PFhat)*Uk$VB - Uk$CB)/Phat$PFhat)
Phat.k_1 <- Phat
Phat <- list(PRhat = Pk.R, PFhat = Pk.F)
NPL.out <- list(model=out.NPL, start=x0, Phat=Phat)

save(list=c("PL.out", "NPL.out"), file="appendixI3_output_1yearsT1.Rdata")

