############################################
############################################
######## Estimating Signaling Games ########
########     Economic Sanctions     ########
########       2step code           ########
############################################
############################################


rm(list=ls())
library(sigInt)
library(data.table)
load("SanctionsDataSet_prelevant.rdata")
missing <- apply(is.na(Xall), 1, max)

Y <- Yall[,missing==0]
X <- Xall[missing==0,]


# create regression matrices
f1 <- sq+cd+sf+bd~ sqrt(senderecondep) + senderdemocracy + contig + ally -1|#SA
                   anticipatedsendercosts|#VA
                   sqrt(targetecondep) + anticipatedtargetcosts + contig + ally|#CB
                   sqrt(senderecondep) + senderdemocracy + lncaprat | #barWA
                   targetdemocracy + lncaprat| #barWB
                   senderdemocracy| #bara
                   -1 #VB

f2 <- ~sqrt(senderecondep) + senderdemocracy + contig + ally + 
  anticipatedsendercosts + sqrt(targetecondep) + anticipatedtargetcosts + 
  lncaprat + targetdemocracy 




X[,`:=`(sq=Y[1,],
        cd=Y[2,],
        sf=Y[3,],
        bd=Y[4,])]

cat("PL:\n") # PL column of Table 7
pl.out <- sigint(f1, data=X, phat.formulas = f2, method="pl", pl.vcov=25)
print(round(summary(pl.out)$coef[19:20,],2))
cat(paste("PL log Likelihood:", round(logLik(pl.out),2), "\n"))

cat("NPL:\n") # NPL column of Table 7
npl.out <- sigint(f1, data=X, phat.formulas = f2, method="npl")
print(round(summary(npl.out)$coef[19:20,],2))
cat(paste("NPL log Likelihood:", round(logLik(npl.out),2), "\n"))



