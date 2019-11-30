suppressMessages(library(ggplot2))
suppressMessages(library(data.table))
suppressMessages(library(matrixStats))
suppressMessages(library(scales))

rm(list=ls())

load("appendixE_Unique_out.rdata")
load("baseLinePlotDataUnique.rdata")

truth1 <- c(1, -1.9, -2.9,.1, -1.2,1)
truth2 <- c(1, -1.7, -2,.1, -1.2,1)
truth <- .5*truth1 + .5*truth2
rmse <- function(bias, var){
  return(sqrt(sum(var)+ crossprod(bias)))
}

len <- length(Results)
K <- ncol(Results[[1]])
twoStepOutBias <- matrix(0, nrow=len, ncol=length(truth))
twoStepOutVar <- matrix(0, nrow=len, ncol=length(truth))
twoStepOutInfo <- matrix(0, nrow=len, ncol=3)

nfxpOutBias <- matrix(0, nrow=len, ncol=length(truth))
nfxpOutVar <- matrix(0, nrow=len, ncol=length(truth))
nfxpOutInfo <- matrix(0, nrow=len, ncol=3)


nfxpTruthBias <- matrix(0, nrow=len, ncol=length(truth))
nfxpTruthVar <- matrix(0, nrow=len, ncol=length(truth))
nfxpTruthInfo <- matrix(0, nrow=len, ncol=3)


for(i in 1:length(Results)){
  twoStep <- Results[[i]][1:9,]
  drop <- which(apply(twoStep[1:6,], 2, function(x){any(abs(x)>100)}))
  twoStep[7,drop] <- -99 
  twoStep[7,is.na(twoStep[7,])] <- -99
  twoStep <- twoStep[-7,twoStep[7,]==1 | twoStep[7,]==2] #0 for bfgs, 1-2 for NR
  twoStepOutBias[i, ]  <- (rowMeans(twoStep[1:length(truth),]-truth))
  twoStepOutVar[i, ]  <-  rowVars(twoStep[1:length(truth),])
  twoStepRMSE <-  rmse(twoStepOutBias[i, ],twoStepOutVar[i, ])
  twoStepOutInfo[i,] <-  c(ncol(twoStep)/K, median(twoStep['elapsed',]), twoStepRMSE)
     
  nfxp <- Results[[i]][10:18,]
  drop <- which(apply(nfxp[1:6,], 2, function(x){any(abs(x)>100)}))
  nfxp[7,drop] <- -99 
  nfxp[7,is.na(nfxp[7,])] <- -99
  nfxp <- nfxp[-7,nfxp[7,]==0]
  nfxpOutBias[i, ]  <- (rowMeans(nfxp[1:length(truth),]-truth))
  nfxpOutVar[i, ]  <- rowVars(nfxp[1:length(truth),])
  nfxpRMSE <- rmse(nfxpOutBias[i, ],nfxpOutVar[i, ])
  nfxpOutInfo[i,] <-  c(ncol(nfxp)/K, median(nfxp['elapsed',]), nfxpRMSE)

  nfxp <- Results[[i]][19:27,]
  drop <- which(apply(nfxp[1:6,], 2, function(x){any(abs(x)>100)}))
  nfxp[7,drop] <- -99 
  nfxp[7,is.na(nfxp[7,])] <- -99
  nfxp <- nfxp[-7,nfxp[7,]==0]
  nfxpTruthBias[i, ]  <- (rowMeans(nfxp[1:length(truth),]-truth))
  nfxpTruthVar[i, ]  <- rowVars(nfxp[1:length(truth),])
  nfxpRMSE <- rmse(nfxpTruthBias[i, ],nfxpTruthVar[i, ])
  nfxpTruthInfo[i,] <-  c(ncol(nfxp)/K, median(nfxp['elapsed',]), nfxpRMSE)
  
}




plotData <- data.table(Mnames = factor(paste(Bparams[,1], "Games"), 
                                       levels=c("25 Games", "50 Games",  "100 Games", "200 Games")),
                       M = Bparams[,1],
                       Observations = Bparams[,2],
                       RMSE = c(twoStepOutInfo[,3],
                                nfxpTruthInfo[,3],
                                nfxpOutInfo[,3]
                       ),
                       Estimator =  rep(c("PL",
                                          "tML-Truth",
                                          "tML-PL"
                       ), 
                       each=len))

plotData <- rbind(plotData, baseLinePlotData.unique[Estimator %in% c("NPL", "CMLE")])
plotData[,Estimator := factor(Estimator, levels=c("CMLE", "PL", "tML-Truth", "tML-PL", "NPL"))]

#RMSE
pMSE <- ggplot(plotData, aes(x=Observations, y=log(RMSE)))+
  geom_line(aes(color=Estimator, linetype=Estimator), size=1.25)+
  facet_wrap(~Mnames, ncol=2)+
  theme_bw(16)+
  ylab('Logged RMSE')+
  xlab('Within-game Observations')+
  scale_linetype_manual(values=c("dotted", "solid", "dashed", "longdash",  "dotdash"))+
  scale_color_manual(values=hue_pal()(5))+
  scale_x_continuous(breaks=unique(Bparams[,2]))+
  theme(legend.position="bottom",
        legend.title = element_text(size=16, 
                                    face="bold"),
        legend.text = element_text(size = 16),
        legend.key.size = unit(.55,"in"))


ggsave(pMSE,  file="figure18.pdf", height=6, width=10)

