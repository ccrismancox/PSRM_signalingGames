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



set.seed(b)


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
xL <- runif(length(Phat$PRhat),-1, 1)
ptm <- proc.time()[3]

out.2step <- try(maxLik(start=x0[1:20], logLik=fqll, grad=gr.qll, method="NR"))
out.2step$time <- proc.time()[3] - ptm

if(class(out.2step[[1]])=="character"){
  out.2step$par <- rep(NA, 20)
  out.2step$convergence <- -99
}else{
  out.2step$par <- out.2step$est
  out.2step$convergence <- out.2step$code
} 
if(!any(is.na(out.2step$par))){
  x0 <- c(out.2step$par,qlogis(Phat$PRhat))
}
