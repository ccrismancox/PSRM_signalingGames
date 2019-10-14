############################################
############################################
######## Estimating Signaling Games ########
########     Economic Sanctions     ########
########       PL, NPL, tML         ########
############################################
############################################


rm(list=ls())

suppressMessages(library(pbivnorm))
suppressMessages(library(Formula))
suppressMessages(library(foreign))
suppressMessages(library(randomForest))
suppressMessages(library(data.table))
suppressMessages(library(maxLik))
suppressMessages(library(numDeriv))
suppressMessages(library(rootSolve))

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

varnames <- paste(names(unlist(lapply(regr, colnames))), ": ", unlist(lapply(regr, colnames)), sep="")



load("estimation_output.Rdata")


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
cat("traditional ML--Newton Solver:\n") #tML column of Table 3
print(out.ML[19:20,])
cat(paste("tML log Likelihood:", round(out.ml$max,2),"\n"))




LL.nfxp.max <- function(x, Y,regr, control=list(ftol=1e-6,maxit=500)){
  M <- dim(Y)[2]
  U <-vec2U.regr(x,regr)
  U$sig <- rep(1, M)
  
  root <- rep(0,M)
  eps <- 0 #sqrt(.Machine$double.eps)
  for(i in 1:M){
    fi <- function(p){const.jo(p, lapply(U,function(x){x[i]}))}
    out <- uniroot.all(fi, c(eps,1-eps))
    out <- ifelse(length(out)>1, max(out), out[1])
    root[i] <- out
  }

  EQ <- eqProbs(root,U, RemoveZeros = T)
  OUT <- cbind(1-EQ[,2], EQ[,2]*(1-EQ[,1]), 	
               EQ[,2]*EQ[,1]*EQ[,3], EQ[,2]*EQ[,1]*(1-EQ[,3]))
  OUT[OUT <= sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
  LL <- sum(log(t(OUT))*Y)
  return(-LL)
}


naiveML2 <- function(x){-LL.nfxp.max(x,Y=Y,regr=regr)}
set.seed(15) #successful convergence
x1 <- rnorm(20)*0.05
# x1 <- c(-0.380896556, 0.031169976, 0.276266653, 0.381151024, -0.143802653, -1.706572959, 1.077894058,
# 0.230771078, -0.120054410, -0.146945740, -0.372871231, -0.413424361, 0.302676740, 0.056092323,
# 0.066489059, 1.191130805, 0.004102281, -0.053892534, -0.756271095, 0.063399313)
out.ml.max <- maxLik::maxLik(start=x1,
                             naiveML2,
                             method="NM",
                             control=list(tol=1e-15,
                                          reltol=1e-15,
                                          gradtol=1e-15,
                                          iterlim=300000))




naiveML <- function(x){-LL.nfxp(x,Y=Y,regr=regr)}
set.seed(15) #successful convergence
x1 <- PL.out$model$estimate
out.ml.pl <- maxLik::maxLik(start=x1,
                             naiveML,
                             method="NM",
                             control=list(tol=1e-15,
                                          reltol=1e-15,
                                          gradtol=1e-15,
                                          iterlim=300000))




####nMLE--MAX####
out.ml <- out.ml.max
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
cat("traditional ML-- Select largest EQ:\n") #Model 5 of Table 5
print(out.ML[19:20,])
cat(paste("tML log Likelihood:", round(out.ml$max,2),"\n"))


####nMLE####
out.ml <- out.ml.pl
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
cat("traditional ML--PL start values:\n") #Model 6 of Table 5
print(out.ML[19:20,])
cat(paste("tML log Likelihood:", round(out.ml$max,2),"\n"))

