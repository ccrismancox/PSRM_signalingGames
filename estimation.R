############################################
############################################
######## Estimating Signaling Games ########
########     Economic Sanctions     ########
########       PL, NPL, tML         ########
############################################
############################################


rm(list=ls())

library(pbivnorm)
library(Formula)
library(foreign)
library(randomForest)
library(data.table)
library(maxLik)
library(rootSolve)

source("signalingFunctions_main.r")
source("gradientFunctions.r")
load("SanctionsDataSet.rdata")


# Find and remove missing values
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



##### Estimating Phat ####
index1 <- colSums(Y[2:4,]) >= 1 #observations where B chooses resist or not
data1 <- cbind((colSums(Y[3:4,])/colSums(Y[2:4,]))[index1], X[index1,])  #data where B chooses resist or not
index2 <- colSums(Y[3:4,]) >= 1 #observations where A chooses SF or BD
data2 <- cbind(((Y[3,])/colSums(Y[3:4,]))[index2], X[index2,])  #data where A chooses SF or BD
colnames(data1)[1] <- "V1"
colnames(data2)[1] <- "V1"

#data transformations like in f1
data1[,`:=`(senderecondep=sqrt(senderecondep), 
             targetecondep=sqrt(targetecondep))]

data2[,`:=`(senderecondep=sqrt(senderecondep),
             targetecondep=sqrt(targetecondep))]

####based on tuning####
set.seed(12345)
m1 <- randomForest(V1~ senderecondep + senderdemocracy + contig + ally +
                          anticipatedsendercosts + targetecondep + anticipatedtargetcosts +
                          lncaprat + targetdemocracy,
                        data = data1,
                        ntree= 1000,
                        mtry=3)
set.seed(12345)
m2 <- randomForest(V1~ senderecondep + senderdemocracy + contig + ally +
                         anticipatedsendercosts + targetecondep + anticipatedtargetcosts +
                          lncaprat + targetdemocracy,
                        data = data2,
                        ntree=1000,
                        mtry=3)
PRhat <- predict(m1, newdata=X)
PFhat <- predict(m2, newdata=X)

##### estimate with the  PL #####
fqll <- function(x){-QLL.jo(x,PRhat,PFhat,Y,regr)}
grqll <- function(x){-eval_gr_qll(x, PRhat, PFhat,Y,regr)}

set.seed(12345)
x0 <- rnorm(20)*.05
out <- maxLik::maxLik(start=x0, fqll,
                      gr=grqll,
                      method='NR',
                      control=list(tol=1e-10,
                                   reltol=1e-10,
                                   gradtol=1e-10,
                                   iterlim=5000))
Phat.PL=list(PRhat=PRhat, PFhat=PFhat)
PL.out <- list(model=out, start=x0, Phat=Phat.PL)
print(out)


##### Estimate with the NPL ####
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
  fqll <- function(x){-QLL.jo(x,Phat$PRhat,Phat$PFhat,Y,regr)}
  grqll <- function(x){-eval_gr_qll(x, Phat$PRhat, Phat$PFhat,Y,regr)}

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
out.NPL$iter <- iter


Pk.F <- eqProbs(Phat$PRhat, Uk, RemoveZeros = T)[,3]
Pk.R <- pnorm((Phat$PFhat*Uk$barWB + (1-Phat$PFhat)*Uk$VB - Uk$CB)/Phat$PFhat)
Phat <- list(PRhat = Pk.R, PFhat = Pk.F)
NPL.out <- list(model=out.NPL, start=x0, Phat=Phat)
print(out.NPL)

##### estimate with traditional ML ####
tML <- function(x){-LL.nfxp(x,Y=Y,regr=regr)}
set.seed(15) #first seed value that returned a successful convergence code
x1 <- rnorm(20)*0.05
# load("tmlstart.rdata")
# x1 <- tml.start
out.ml <- maxLik::maxLik(start=x1,
                         tML,
                         method="NM",
                         control=list(tol=1e-15,
                                      reltol=1e-15,
                                      gradtol=1e-15,
                                      iterlim=300000)
)

tML.out <- list(model=out.ml, start=x1)
save(list=c("PL.out","tML.out", "NPL.out"), file="estimation_output.Rdata")

print(out.ml)
warning("End of file. Press enter if the system hangs here.")