############################################
############################################
######## Estimating Signaling Games ########
########   Monte Carlo Experiment   ########
########     Multiple Equilibria    ########
############################################
############################################

suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
suppressMessages(library(scales))

rm(list=ls())


## Read in the CMLE results
load("CMLE_meq.rdata")
IPOPT <- Results

## Reading PL, NPL, and tML results
load("MonteCarloResults_MEQ.rdata")

## Reshape the CMLE results to match the others
IPOPT <- lapply(IPOPT, function(x){x <- cbind(x, x[,7]); x<- x[,-7]; return(x)})
IPOPT <- lapply(IPOPT, function(x){return(t(x[,1:8]))})
IPOPT <- lapply(IPOPT, function(x){rownames(x)[nrow(x)] <- "elapsed"; return(x)})

truth <- c(1, -1.9, -2.9,.1, -1.2,1)


rmse <- function(bias, var){
  return(sqrt(sum(var)+ crossprod(bias)))
}



len <- length(Results)
K <- ncol(Results[[1]])
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

for(i in 1:length(Results)){
  twoStep <- Results[[i]][1:9,]
  drop <- which(apply(twoStep[1:6,], 2, function(x){any(abs(x)>50)})) #Check for bizarre results
  twoStep[7,drop] <- -99  
  twoStep <- twoStep[-7,twoStep[7,]==1 | twoStep[7,]==2] #Check for convergence
  twoStepOutBias[i, ]  <- (rowMeans(twoStep[1:length(truth),]-truth))
  twoStepOutVar[i, ]  <-  rowVars(twoStep[1:length(truth),])
  twoStepRMSE <-  rmse(twoStepOutBias[i, ],twoStepOutVar[i, ])
  twoStepOutInfo[i,] <-  c(ncol(twoStep)/K, median(twoStep['elapsed',]), twoStepRMSE)
  
  npl <- Results[[i]][10:18,]
  npl['elapsed',] <- npl['elapsed',] + Results[[i]][8,]
  drop <- which(apply(npl[1:6,], 2, function(x){any(abs(x)>50)})) #Check for bizarre results
  npl[7,drop] <- -99
  npl <- npl[-7,npl[7,]==1| npl[7,]==2 | npl[7,]==-69]  #Check for convergence
  nplOutBias[i, ]  <- (rowMeans(npl[1:length(truth),]-truth))
  nplOutVar[i, ]  <-  rowVars(npl[1:length(truth),])
  nplRMSE <-  rmse(nplOutBias[i, ], nplOutVar[i, ])
  nplOutInfo[i,] <-  c(ncol(npl)/K, median(npl['elapsed',]), nplRMSE)

  ipopt <- IPOPT[[i]]
  drop <- which(apply(ipopt[1:6,], 2, function(x){any(abs(x)>50)})) #Check for bizarre results
  ipopt[7,drop] <- -99
  ipopt <- ipopt[-7,ipopt[7,]==0] # Check for convergence
  cmleOutBias[i, ]  <- (rowMeans(ipopt[1:length(truth),]-truth))
  cmleOutVar[i, ]  <- rowVars(ipopt[1:length(truth),])
  cmleRMSE <- rmse(cmleOutBias[i, ], cmleOutVar[i, ])
  cmleOutInfo[i,] <-  c(ncol(ipopt)/1000, median(ipopt['elapsed',]),cmleRMSE)
  
  
  nfxp <- Results[[i]][19:27,]
  drop <- which(apply(nfxp[1:6,], 2, function(x){any(abs(x)>50)})) #Check for bizarre results
  nfxp[7,drop] <- -99
  nfxp <- nfxp[-7,nfxp[7,]==0]  # Check for convergence
  nfxpOutBias[i, ]  <- (rowMeans(nfxp[1:length(truth),]-truth))
  nfxpOutVar[i, ]  <- rowVars(nfxp[1:length(truth),])
  nfxpRMSE <- rmse(nfxpOutBias[i, ],nfxpOutVar[i, ])
  nfxpOutInfo[i,] <-  c(ncol(nfxp)/K, median(nfxp['elapsed',]), nfxpRMSE)
}




plotData <- data.table(Mnames = factor(paste(Bparams[,1], "Games"), 
                                       levels=c("25 Games", "50 Games",  "100 Games", "200 Games")),
                       M = Bparams[,1],
                       Observations = Bparams[,2],
                       Bias = c(abs(rowMeans(twoStepOutBias)),
                                abs(rowMeans(nplOutBias)),
                                abs(rowMeans(cmleOutBias)),
                                abs(rowMeans(nfxpOutBias))
                       ),
                       Var =  c(rowMeans((twoStepOutVar)),
                                rowMeans((nplOutVar)),
                                rowMeans((cmleOutVar)),
                                rowMeans((nfxpOutVar))
                       ),
                       Time = c(twoStepOutInfo[,2],
                                nplOutInfo[,2],
                                cmleOutInfo[,2],
                                nfxpOutInfo[,2]
                       ),
                       Converged = c(twoStepOutInfo[,1],
                                     nplOutInfo[,1],
                                     cmleOutInfo[,1],
                                     nfxpOutInfo[,1]
                       ),
                       Estimator =  rep(c("PL",
                                          "NPL",
                                          "CMLE", 
                                          "tML"
                       ), 
                       each=len))



plotData[,Estimator := factor(Estimator, 
                              levels=c("CMLE",
                                       "PL",
                                       "tML",
                                       "NPL"))]

#bias
pBias <- ggplot(plotData, aes(x=Observations, y=Bias))+
  geom_line(aes(color=Estimator, linetype=Estimator), size=1.25)+
  facet_wrap(~Mnames, ncol=2)+
  theme_bw(16)+
  ylab('Absolute Mean Bias')+
  xlab('Within-game observations')+
  scale_linetype_manual(values=c("dotted", "solid", "dashed",  "dotdash"))+
  scale_color_manual(values=hue_pal()(4))+
  scale_x_continuous(breaks=unique(Bparams[,2]))+
  theme(legend.position="bottom",
        legend.title = element_text(size=18, 
                                    face="bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(.65,"in"))

#Var
pVar <- ggplot(plotData, aes(x=Observations, y=log(Var)))+
  geom_line(aes(color=Estimator, linetype=Estimator), size=1.25)+
  facet_wrap(~Mnames, ncol=2)+
  theme_bw(16)+
  ylab('Logged Mean Variance')+
  xlab('Within-game observations')+
  scale_linetype_manual(values=c("dotted", "solid", "dashed",  "dotdash"))+
  scale_color_manual(values=hue_pal()(4))+
  scale_x_continuous(breaks=unique(Bparams[,2]))+
  theme(legend.position="bottom",
        legend.title = element_text(size=18, 
                                    face="bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(.65,"in"))


#Converged
pConv <- ggplot(plotData, aes(x=Observations, y=Converged))+
  geom_line(aes(color=Estimator, linetype=Estimator), size=1.25)+
  facet_wrap(~Mnames, ncol=2)+
  theme_bw(16)+
  ylab("Proportion Converged")+
  xlab('Within-game observations')+
  scale_linetype_manual(values=c("dotted", "solid", "dashed",  "dotdash"))+
  scale_color_manual(values=hue_pal()(4))+
  scale_x_continuous(breaks=unique(Bparams[,2]))+
  theme(legend.position="bottom",
        legend.title = element_text(size=18, 
                                    face="bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(.65,"in"))

#Time
pTime <- ggplot(plotData, aes(x=Observations, y=log(Time)))+
  geom_line(aes(color=Estimator, linetype=Estimator), size=1.25)+
  facet_wrap(~Mnames, ncol=2)+
  theme_bw(16)+
  ylab("logged CPU Time (second)")+
  xlab('Within-game observations')+
  scale_linetype_manual(values=c("dotted", "solid", "dashed",  "dotdash"))+
  scale_color_manual(values=hue_pal()(4))+
  scale_x_continuous(breaks=unique(Bparams[,2]))+
  theme(legend.position="bottom",
        legend.title = element_text(size=18, 
                                    face="bold"),
        legend.text = element_text(size = 18),
        legend.key.size = unit(.65,"in"))



ggsave(pBias, file="figure8.pdf", height=6, width=10)
ggsave(pVar,  file="figure9.pdf",  height=6, width=10)
ggsave(pConv, file="figure10.pdf", height=6, width=10)
ggsave(pTime, file="figure11.pdf", height=6, width=10)
