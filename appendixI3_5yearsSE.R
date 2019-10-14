##########################################################
## Appendix I3: Combine Models and Create Standard Errors
##########################################################

rm(list=ls())

suppressMessages(library(pbivnorm))
suppressMessages(library(Formula))
suppressMessages(library(foreign))
suppressMessages(library(randomForest))
suppressMessages(library(data.table))
suppressMessages(library(numDeriv))
suppressMessages(library(Matrix))	
suppressMessages(library(maxLik))
suppressMessages(library(rootSolve))
suppressMessages(library(e1071))

ignore.case <- function(string) {
  structure(string, ignore.case = TRUE)
}

source("signalingFunctions_main.r")
source("gradientFunctions.r")
load("SanctionsDataSet_5years.rdata")
load("appendixI3_output_5years.Rdata")
load("appendixI3_CMLEoutput_5years.Rdata")


missing <- apply(is.na(Xall), 1, max)

Y <- Yall[,missing==0]
X <- Xall[missing==0,]

M <- ncol(Y)
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


####bootstrap the first-stage covariance####
set.seed(12345)
Phat.boot <- matrix(0, ncol=M*2, nrow=50)
for(i in 1:50){
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

  m1 <- randomForest(V1~ senderecondep + senderdemocracy + contig + ally +
                     anticipatedsendercosts + targetecondep + anticipatedtargetcosts +
                     lncaprat + targetdemocracy,
                   data = polyX1.boot,
                   ntree=c(1000), mtry=c(3))
  m2 <- randomForest(V1~ (senderecondep) + senderdemocracy + contig + ally +
                              anticipatedsendercosts + (targetecondep) +
                              anticipatedtargetcosts +lncaprat + targetdemocracy,
                   data = polyX2.boot,
                   ntree=c(1000), mtry=c(3))

  Phat.boot[i,] <- c(predict(m1, newdata=X), predict(m2, newdata=X))
}

SIGMA <- cov(Phat.boot)
save(list="SIGMA", file="SIGMA_5year.rdata")
load("SIGMA_5year.rdata")


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
cat("PL:\n") #The PL column ofTable 8
print(out.PL[19:20,])
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
out.npl <- round(cbind(sign(NPL.hat), NPL.hat, se.NPL, z, p),2)
colnames(out.npl) <- c("Sign", "Est.", "S.E.", "z-stat", "p<0.05")
cat("NPL:\n") # NPL column ofTable 8
print(out.npl[19:20,])
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
cat("CMLE:\n") #These are the values for the CMLE estimates inTable 8
print(out.CMLE[19:20,])
cat(paste("CMLE log Likelihood:", round(LL.cmle,2),"\n"))


