############################################
############################################
######## Estimating Signaling Games ########
########     Economic Sanctions     ########
########       Standard errors      ########
############################################
############################################

rm(list=ls())

library(pbivnorm)
library(Formula)
library(foreign)
library(randomForest)
library(data.table)
library(numDeriv)
library(Matrix)	
library(maxLik)
library(rootSolve)
library(e1071)


source("signalingFunctions_main.r")
source("gradientFunctions.r")
load("SanctionsDataSet.rdata")
load("estimation_output.Rdata")
load("CMLE_estimation_output.rdata")


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

varnames <- paste(names(unlist(lapply(regr, colnames))), ": ", unlist(lapply(regr, colnames)), sep="")

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

####bootstrap the first-stage covariance####
set.seed(12345)
Phat.boot <- matrix(0, ncol=418*2, nrow=500)
for(i in 1:500){
   use <- sample(1:418, replace=TRUE)
   Y.boot <- Y[,use]
   X.boot <- X[use,]
   index1.boot <- colSums(Y.boot[2:4,]) >= 1
   data1.boot <- cbind((colSums(Y.boot[3:4,])/colSums(Y.boot[2:4,]))[index1.boot],
                        X.boot[index1.boot,])
   index2.boot <- colSums(Y.boot[3:4,]) >= 1
   data2.boot <- cbind(((Y.boot[3,])/colSums(Y.boot[3:4,]))[index2.boot], X.boot[index2.boot,])
   colnames(data1.boot)[1] <- "V1"
   colnames(data2.boot)[1] <- "V1"
   data1.boot[,`:=`(senderecondep=sqrt(senderecondep),
                     targetecondep=sqrt(targetecondep))]
   data2.boot[,`:=`(senderecondep=sqrt(senderecondep),
                     targetecondep=sqrt(targetecondep))]
 
   m1 <- tune.randomForest(V1~ senderecondep + senderdemocracy + contig + ally +
                      anticipatedsendercosts + targetecondep + anticipatedtargetcosts +
                      lncaprat + targetdemocracy,
                    data = data1.boot,
                    ntree=c(50, 500, 1000), mtry=c(3, 6,9)
   )$best.model
   m2 <- tune.randomForest(V1~ (senderecondep) + senderdemocracy + contig + ally +
                               anticipatedsendercosts + (targetecondep) +
                               anticipatedtargetcosts +lncaprat + targetdemocracy,
                    data = data2.boot,
                    ntree=c(50, 500, 1000), mtry=c(3, 6,9)
   )$best.model
 
   Phat.boot[i,] <- c(predict(m1, newdata=X), predict(m2, newdata=X))
}
  
SIGMA <- var(Phat.boot)
save(list="SIGMA", file="SIGMA.rdata")
load("SIGMA.rdata")


####nMLE####
out.ml <- tML.out$model
##functions
nMLE.i <- function(x){LL.nfxp.i(x, Y=Y, regr=regr)}
Dtheta <- jacobian(nMLE.i, out.ml$est)
Dtheta.theta <- crossprod(Dtheta)
vcov.ml <- solve(Dtheta.theta)

se.ml <- sqrt(diag(vcov.ml))

ML.hat <- out.ml$estimate
z <- ML.hat/se.ml
p <- 2*pnorm(abs(z),lower.tail = F)<0.05
names(ML.hat) <- varnames
out.ML <- round(cbind(sign(ML.hat), ML.hat, se.ml, z, p),2)
colnames(out.ML) <- c("Sign", "Est.", "S.E.", "z-stat", "p<0.05")
cat("traditional ML:\n") #tML column of Table 3
print(out.ML)
cat(paste("tML log Likelihood:", round(out.ml$max,2),"\n"))


####PL####
out <- PL.out$model
Phat.PL <- PL.out$Phat
#create matrices
Dtheta <- eval_gr_qll.i(out$estimate, Phat.PL$PRhat, Phat.PL$PFhat, Y, regr)
Dtheta.theta <- crossprod(Dtheta)
Dtheta.theta.inv <- solve(Dtheta.theta)

Dp <- eval_gr_qll.ip(Phat.PL$PRhat, Phat.PL$PFhat, out$est, Y, regr)
Dtheta.p <- crossprod(Dtheta,Dp)


vcov.PL <- Dtheta.theta.inv + Dtheta.theta.inv %*% Dtheta.p %*% SIGMA %*% t(Dtheta.p) %*% Dtheta.theta.inv
se.PL <- sqrt(diag(vcov.PL))
PL.hat <- out$estimate
z <- PL.hat/se.PL
p <- 2*pnorm(abs(z),lower.tail = F)<0.05
names(PL.hat) <- varnames
out.PL <- round(cbind(sign(PL.hat), PL.hat, se.PL, z, p),2)
colnames(out.PL) <- c("Sign", "Est.", "S.E.", "z-stat", "p<0.05")
cat("PL:\n") #The PL column of Table 3
print(out.PL)
cat(paste("PL log Likelihood:", round(out$max,2), "\n"))

####NPL####
out.NPL <- NPL.out$model
Phat <- NPL.out$Phat
Dtheta <- eval_gr_qll.i(out.NPL$estimate,Phat$PRhat, Phat$PFhat, Y, regr )
Dtheta.theta <- crossprod(Dtheta)

Dp <- eval_gr_qll.ip(Phat$PRhat, Phat$PFhat, out.NPL$est, Y, regr)
Dtheta.p <- crossprod(Dtheta,Dp)

JpPsi <- dPsiDp(Phat$PRhat, Phat$PFhat, out.NPL$est, Y, regr)
JtPsi <- dPsi.dTheta(out.NPL$estimate,Phat$PRhat, Phat$PFhat, Y, regr)

topSlice <- solve(Dtheta.theta  + Dtheta.p %*% solve(diag(nrow(JpPsi)) - t(JpPsi)) %*% JtPsi)
bottomSlice <- solve(Dtheta.theta  + t(JtPsi) %*% solve(diag(nrow(JpPsi)) - t(JpPsi)) %*% t(Dtheta.p))
vcov.NPL <- topSlice %*% Dtheta.theta %*% bottomSlice
se.NPL <- sqrt(diag(vcov.NPL))

NPL.hat <- out.NPL$estimate
z <- NPL.hat/se.NPL
p <- 2*pnorm(abs(z),lower.tail = F)<0.05
names(NPL.hat) <- varnames
out.NPL <- round(cbind(sign(NPL.hat), NPL.hat, se.NPL, z, p),2)
colnames(out.NPL) <- c("Sign", "Est.", "S.E.", "z-stat", "p<0.05")
cat("NPL:\n") # NPL column of Table 3
print(out.NPL)
cat(paste("NPL log Likelihood:", round(out.NPL$max,2), "\n"))




####CMLE####
results <- c(output)
LL.cmle <- -value

gr <- function(x){eval_gr_obj(x,Y,regr)}
heq.jac <- function(x){Jconst(x,Y,regr)}

# SE stuff
Btheta <- jacobian(gr,results) # this is hessian of likelihood.
Htheta <- t(heq.jac(results))
W <- Btheta + Htheta %*% t(Htheta)
top <- cbind(W, -Htheta)
bottom <- cbind(-t(Htheta), matrix(0, dim(Htheta)[2], dim(Htheta)[2]))
BorderHessian <- rbind(top,bottom)
BorderHessian <- Matrix(BorderHessian, sparse=T)
invBorderHessian <- solve(BorderHessian)
se.cmle <- sqrt(diag(invBorderHessian)[1:20])
CMLE.hat <- results[1:20]

z <- CMLE.hat/se.cmle
p <- 2*pnorm(abs(z),lower.tail = FALSE)<0.05
names(CMLE.hat) <- varnames
out.CMLE <- round(cbind(sign(CMLE.hat), CMLE.hat, se.cmle, z, p),2)
colnames(out.CMLE) <- c("Sign", "Est.", "S.E.", "z-stat", "p<0.05")
cat("CMLE:\n") #These are the values for the CMLE estimates in Table 3
print(out.CMLE)
cat(paste("CMLE log Likelihood:", round(LL.cmle,2),"\n"))

warning("End of file. Press enter if the system hangs here.")


