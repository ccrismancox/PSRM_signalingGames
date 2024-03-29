############################################
############################################
######## Estimating Signaling Games ########
########     Economic Sanctions     ########
########        Comp Stat           ########
############################################
############################################

rm(list=ls())

suppressMessages(library(pbivnorm))
suppressMessages(library(Formula))
suppressMessages(library(foreign))
suppressMessages(library(ggplot2))
suppressMessages(library(doRNG))
suppressMessages(library(rootSolve))
suppressMessages(library(doParallel))

############################################
# Functions
source("signalingFunctions_main.r")
source("gradientFunctions.r")

############################################
# Data
load("SanctionsDataSet.rdata")


# Find and remove missing values
missing <- apply(is.na(Xall), 1, max)
Y <- Yall[,missing==0]
X <- Xall[missing==0,]
M <- dim(Y)[2]

rm(Yall, Xall, missing)


############################################
# CMLE
load("CMLE_estimation_output.rdata") 
out.cmle <- c(output)
eqCMLE <- plogis(out.cmle[(length(out.cmle)-M+1):length(out.cmle)])
thetaCMLE <- out.cmle[1:(length(out.cmle)-M)]

rm(value, out.cmle, output)


############################################
# prepare data for counterfactuals

# create regression matrices
f1 <- as.Formula(~ sqrt(senderecondep) + senderdemocracy + contig + ally -1|#SA
                   anticipatedsendercosts|#VA
                   sqrt(targetecondep) + anticipatedtargetcosts + contig + ally|#CB
                   sqrt(senderecondep) + senderdemocracy + lncaprat | #barWA
                   targetdemocracy + lncaprat| #barWB
                   senderdemocracy| #bara
                   -1) #VB

Xhi <- Xlo <- X

regrhi <- list()   
regrlo <- list() 
regr <- list()   
for(i in 1:7){
  regrhi[[i]] <- model.matrix(f1, data=Xhi, rhs=i)
  regrlo[[i]] <- model.matrix(f1, data=Xlo, rhs=i)
  regr[[i]] <- model.matrix(f1, data=X, rhs=i)

}
names(regrhi) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")
names(regrlo) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")
names(regr) <- c("SA", "VA", "CB", "barWA", "barWB", "bara", "VB")





############################################
# prepare data for counterfactuals

# cmle
xThi <- xTlo <- thetaCMLE
xThi[19] <- 0
xTlo[19] <- -6
xThi[20] <- 0
xTlo[20] <- 0

Uhi.cmle <- vec2U.regr(xThi, regr) 
Ulo.cmle <- vec2U.regr(xTlo, regr)
U.cmle <- vec2U.regr(thetaCMLE, regr) 

Uhi.cmle$sig <- Ulo.cmle$sig <- U.cmle$sig <- rep(1,M) 

info <- data.frame(gameID=X$gameID[93],
                   sender = "USA",
                   receiver = "CHN",
                   year = X$tenyear[93])

############################################
# counterfactuals
use <- 93
es <- function(x,k){return(x[k])}
nstep <- 100
out.cmle <- signalEQC(lapply(Ulo.cmle,es,k=use), lapply(Uhi.cmle,es,k=use), length.out=nstep,gridsize=1e4,comp=T)


ggdata <- data.frame(x = c(out.cmle[,"i"]),
                     eq = c(out.cmle[,"sols"]),
                     pc = c(out.cmle[,"sols.pc"]),
                     pf = c(out.cmle[,"sols.pf"]),
                     pr = c(out.cmle[,"sols"]),
                     pars = c(rep("CMLE", dim(out.cmle)[1])))
ggdata$psanction <- ggdata$pc*ggdata$pr*ggdata$pf

vldata <- data.frame(vl = (nstep/6)*c(thetaCMLE[19] + thetaCMLE[20]*X$senderdemocracy[use]) + nstep,
                     pars = c("CMLE"))
                            
pR <- eqCMLE
pC <- g.jo(cStar.jo(pR,U.cmle),U.cmle)
pF <- h.jo(cStar.jo(pR,U.cmle),U.cmle)/pC
psanction <- pC*pR*pF

                     
pdf("figure21.pdf", height=6, width=10)
ggplot(ggdata, aes(x=x,y=pf)) + geom_point(color="navyblue",size=3.5) + 
  scale_x_continuous(breaks=seq(from=0,to=nstep,length.out=7), 
                     labels = round(seq(from=-6,to=0,by=1), digits=1)) +
  theme_bw(18) + 
  geom_point(aes(x=vldata$vl-1, y=pF[93]),color= "orangered1", shape=18, size=5) + 
  xlab( bquote(.(paste(info$sender, "Audience Costs,")) ~ bar(italic(a)))) +   
  ylab(expression(paste("Prob. of Fighting, ", italic(p[F]))))
dev.off()

pdf("figure22.pdf", height=6, width=10)
ggplot(ggdata, aes(x=x,y=pr)) + geom_point(color="navyblue",size=3.5) + 
  scale_x_continuous(breaks=seq(from=0,to=nstep,length.out=7), 
                     labels = round(seq(from=-6,to=0,by=1), digits=1)) +
  theme_bw(18) + 
  geom_point(aes(x=vldata$vl-1, y=pR[93]),color= "orangered1", shape=18, size=5) + 
  xlab( bquote(.(paste(info$sender, "Audience Costs,")) ~ bar(italic(a)))) +   
  ylab(expression(paste("Equilibria, ", italic(p[R]))))
dev.off()


pdf("figure23.pdf", height=6, width=10)
ggplot(ggdata, aes(x=x,y=psanction)) + geom_point(color="navyblue",size=3.5) + 
  scale_x_continuous(breaks=seq(from=0,to=nstep,length.out=7), 
                     labels = round(seq(from=-6,to=0,by=1), digits=1)) +
  theme_bw(18) + 
  geom_point(aes(x=vldata$vl-1, y=psanction[93]+.004),color= "orangered1", shape=18, size=5) + #small manual jitter
  xlab( bquote(.(paste(info$sender, "Audience Costs,")) ~ bar(italic(a)))) +   
  ylab("Probability of\nObserving Sanctions")
dev.off()

