suppressMessages(library(doParallel))
suppressMessages(library(doRNG))
suppressMessages(library(foreach))
suppressMessages(library(rootSolve))
suppressMessages(library(pbivnorm))
suppressMessages(library(Formula))
suppressMessages(library(mc2d))
suppressMessages(library(maxLik))
suppressMessages(library(randomForest))

source("signalingFunctions_main.r")
source("gradientFunctions.r")

## Setup the parameters of the experiment
nPerGame <- c(5, 25, 50, 100, 200)
nGames <- c(25, 50, 100, 200)



Bparams <- expand.grid(nGames,nPerGame)
Results <- list()

set.seed(b)

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


# Yr <- colSums(Y[3:4,])/colSums(Y[2:4,])
# Yf <- Y[3,]/colSums(Y[3:4,])
# dat <- cbind(X, Yr, Yf)
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

  
  x0 <- c(runif(6), qlogis(Phat$PRhat))
}
   

xL <- runif(length(Phat$PRhat),-1, 1)

ptm <- proc.time()[3]
out.2step <- try(maxLik(start=x0[1:6], logLik=fqll, grad=gr.qll, method="NR"))
out.2step$time <- proc.time()[3] - ptm

if(class(out.2step[[1]])=="character"){
  out.2step$par <- rep(NA, 6)
  out.2step$convergence <- -99
}else{
  out.2step$par <- out.2step$est
  out.2step$convergence <- out.2step$code
} 
if(!any(is.na(out.2step$par))){
  x0 <- c(out.2step$par,qlogis(Phat$PRhat))
}
