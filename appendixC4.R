############################################
############################################
######## Estimating Signaling Games ########
########   Monte Carlo Experiment   ########
########    Multiple variables      ########
############################################
############################################


suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
rm(list=ls())


load('CMLE_estimation_output.rdata')
results <- c(output)
LL.cmle <- -value
truth <- results[1:20]

load("CMLE_sanctionsMC.rdata")
IPOPT <- Results

IPOPT <- cbind(IPOPT, IPOPT[,21]); IPOPT <- IPOPT[,-21] 
IPOPT <- t(IPOPT)
rownames(IPOPT)[nrow(IPOPT)] <- "elapsed"

rmse <- function(bias, var){
  return(sqrt(sum(var)+ crossprod(bias)))
}

load("AppendixC4_results.RData")

len <- 1
K <- ncol(Results)
twoStepOutBias <- matrix(0, nrow=len, ncol=length(truth))
twoStepOutVar <- matrix(0, nrow=len, ncol=length(truth))
twoStepOutInfo <- matrix(0, nrow=len, ncol=3)

nplOutBias <- matrix(0, nrow=len, ncol=length(truth))
nplOutVar <- matrix(0, nrow=len, ncol=length(truth))
nplOutInfo <- matrix(0, nrow=len, ncol=3)


cmleOutBias <- matrix(0, nrow=len, ncol=length(truth))
cmleOutVar <- matrix(0, nrow=len, ncol=length(truth))
cmleOutInfo <- matrix(0, nrow=len, ncol=3)

nfxpOutBias <- matrix(0, nrow=len, ncol=length(truth))
nfxpOutVar <- matrix(0, nrow=len, ncol=length(truth))
nfxpOutInfo <- matrix(0, nrow=len, ncol=3)

# cap <- 10
twoStep <- Results[1:23,]
drop <- which(apply(twoStep[1:20,], 2, function(x){any(abs(x)>50)}))
twoStep[21,drop] <- -99 
twoStep <- twoStep[-21,twoStep[21,]==1 | twoStep[21,]==2] #0 for bfgs, 1-2 for NR
twoStepOutBias  <- (rowMeans(twoStep[1:length(truth),]-truth))
twoStepOutVar  <-  rowVars(twoStep[1:length(truth),])
twoStepRMSE <-  rmse(twoStepOutBias,twoStepOutVar)
twoStepOutInfo<-  c(ncol(twoStep)/K, median(twoStep['elapsed',]), twoStepRMSE)


npl <- Results[24:46,]
npl['elapsed',] <- npl['elapsed',] + Results[22,]
drop <- which(apply(npl[1:20,], 2, function(x){any(abs(x)>50)}))
npl[21,drop] <- -99
npl <- npl[-21,npl[21,]==1| npl[21,]==2 | npl[21,]==-69]
nplOutBias  <- (rowMeans(npl[1:length(truth),]-truth))
nplOutVar  <-  rowVars(npl[1:length(truth),])
nplRMSE <-  rmse(nplOutBias, nplOutVar)
nplOutInfo <-  c(ncol(npl)/K, median(npl['elapsed',]), nplRMSE)


cmle <- IPOPT
drop <- which(apply(cmle[1:20,], 2, function(x){any(abs(x)>50)}))
cmle[21,drop] <- -99
cmle <- cmle[-21,cmle[21,]==0]
cmleOutBias  <- (rowMeans(cmle[1:length(truth),]-truth))
cmleOutVar  <- rowVars(cmle[1:length(truth),])
cmleRMSE <- rmse(cmleOutBias,cmleOutVar)
cmleOutInfo <-  c(ncol(cmle)/K, median(cmle['elapsed',]), cmleRMSE)

nfxp <- Results[47:69,]
drop <- which(apply(nfxp[1:20,], 2, function(x){any(abs(x)>50)}))
nfxp[21,drop] <- -99
nfxp <- nfxp[-21,nfxp[21,]==0 | nfxp[21,]==1]
nfxpOutBias <- (rowMeans(nfxp[1:length(truth),]-truth))
nfxpOutVar <- rowVars(nfxp[1:length(truth),])
nfxpRMSE <- rmse(nfxpOutBias, nfxpOutVar)
nfxpOutInfo <-  c(ncol(nfxp)/K, median(nfxp['elapsed',]), nfxpRMSE)




plotData <- data.table(Bias = c((twoStepOutBias),
                                (nplOutBias),
                                (cmleOutBias),
                                (nfxpOutBias)
                       ),
                       Var =  c(((twoStepOutVar)),
                                ((nplOutVar)),
                                ((cmleOutVar)),
                                ((nfxpOutVar))
                       ),
                       Estimator =  rep(c("PL",
                                          "NPL",
                                          "CMLE", 
                                          "tMLE"
                       ), 
                       each=20))
plotData[,RMSE := sqrt(Bias^2 + Var)]
plotData[,relRmse := RMSE/plotData[Estimator=="tMLE", RMSE], by=Estimator]

relative <- data.frame(PL = plotData[Estimator=="PL", relRmse],
                       NPL = plotData[Estimator=="NPL", relRmse],
                       CMLE = plotData[Estimator=="CMLE", relRmse])
relative <- rbind(relative,  c(twoStepRMSE, nplRMSE, cmleRMSE)/c(nfxpRMSE))
row.names(relative) <- c("$S_A$: Econ. Dep$_A$", "$S_A$:  Dem$_A$","$S_A$: Contiguity","$S_A$: Alliance",
                         "$V_A$: Const.", "$V_A$: Costs$_A$",
                         "$C_B$: Const.", "$C_B$: Econ. Dep$_B$", "$C_B$: Costs$_B$", "$C_B$: Contiguity","$C_B$: Alliance",
                         "$\\bar{W}_A$: Const.", "$\\bar{W}_A$: Econ. Dep$_A$", "$\\bar{W}_A$: Dem$_A$", "$\\bar{W}_A$: Cap. Ratio",
                         "$\\bar{W}_B$: Const.", "$\\bar{W}_B$: Dem$_B$", "$\\bar{W}_B$: Cap. Ratio",
                         "$\\bar{a}$: Const.", "$\\bar{a}$: Dem$_A$", "Multivariate RMSE")

cat("Table 4:\n")
print(relative)




