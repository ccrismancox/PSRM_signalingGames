rm(list=ls())


suppressMessages(library(data.table))
suppressMessages(library(Formula))

load("SanctionsDataSet_1yearsT12.rdata")
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

load("appendixI3_output_5years.Rdata")
PRhat <- PL.out$Phat$PRhat
out <- PL.out$model
set.seed(15) 
xL <- runif(length(PRhat),-1, 1)
x0 <- c(out$est, qlogis(PRhat))


