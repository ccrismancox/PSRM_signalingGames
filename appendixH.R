##############################
# Analyze whether variables change
#
###############################

rm(list=ls())


suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(gridExtra))

load("SanctionsDataSet_1yearsT12.rdata")


# remove irrelevant and missing dyads values
Y <- Yall
X <- Xall
X <- X[which(colSums(Y) == 12)]
Y <- Y[,which(colSums(Y) == 12)]

X[, `:=` (senderecondep=sqrt(senderecondep),
          targetecondep=sqrt(targetecondep))]





#################################
# RESULTS AT THE COUNTRY LEVEL ##
#################################
# create a country-year matrix
InOut <- data.table(ccode = c(X$ccode1, X$ccode2),
                    year = c(X$year, X$year),
                    nChallenge=c(colSums(Y[-1,]), rep(0, ncol(Y))),
                    nResist=c(rep(0, ncol(Y)),colSums(Y[-c(1,2),])),
                    nSF=c(Y[4,], rep(0, ncol(Y))),
                    polity = c(X$senderdemocracy, X$targetdemocracy)
)
InOut[,decade:=as.numeric(str_c(year %/% 10, 0))]
InOut[,`:=`(nChallenge = sum(nChallenge),
            nResist = sum(nResist),
            nSF = sum(nSF)),
      by=list(ccode, year)]
InOut <- unique(InOut)

InOut[,polityDev := polity-mean(polity,na.rm=T), by=list(ccode, decade)]






########################
# DYADIC VARIABLES     #
########################


Mcowt <- X[,list(dyadID, year, ccode1, ccode2, anticipatedtargetcosts,anticipatedsendercosts,
                 senderecondep,targetecondep, lncaprat,ally)]
Mcowt[,decade:=as.numeric(str_c(year %/% 10, 0))]
Mcowt[,`:=`(targetcostDev = anticipatedtargetcosts-mean(anticipatedtargetcosts, na.rm=T),
            sendercostDev = anticipatedsendercosts-mean(anticipatedsendercosts, na.rm=T),
            senderecondepDev = senderecondep-mean(senderecondep, na.rm=T),
            allyDev = ally-mean(ally, na.rm=T),
            targetecondepDev = targetecondep-mean(targetecondep, na.rm=T),
            lncapratDev = lncaprat-mean(lncaprat, na.rm=T)
), by=list(dyadID, decade)]

Mcowt[,`:=`(nChallenge = colSums(Y[-1,]),
            nResist = colSums(Y[-c(1,2),]),
            nSF=c(Y[4,]))]

pdf("figure24.pdf", width=15)

grid.arrange(
  ggplot(InOut, aes(x=InOut$polityDev, y=..density..)) + 
    geom_histogram(color="black",fill="lightsteelblue",binwidth=0.5) + 
    theme_bw() + 
    xlab("Deviation from Mean Polity2") + 
    ylab("Density") + 
    theme(axis.title.x=element_text(size=14,vjust=-.2),
          axis.title.y=element_text(size=14),
          axis.text=element_text(size=12)
    )
  ,
  
  ggplot(Mcowt, aes(x=Mcowt$targetcostDev, y=..density..)) + 
    geom_histogram(color="black",fill="lightsteelblue",binwidth=.005) + 
    theme_bw() + 
    xlab("Deviation from Mean\n Anticipated Target Costs") + 
    ylab("Density") + 
    theme(axis.title.x=element_text(size=14,vjust=-.2),
          axis.title.y=element_text(size=14),
          axis.text=element_text(size=12)
    )  ,
  
  ggplot(Mcowt, aes(x=Mcowt$sendercostDev, y=..density..)) + 
    geom_histogram(color="black",fill="lightsteelblue",binwidth=.1) + 
    theme_bw() + 
    xlab("Deviation from Mean\n Anticipated Sender Costs") + 
    ylab("Density") + 
    theme(axis.title.x=element_text(size=14,vjust=-.2),
          axis.title.y=element_text(size=14),
          axis.text=element_text(size=12)
    )   ,
  
  
  ggplot(Mcowt, aes(x=Mcowt$senderecondepDev, y=..density..)) + 
    geom_histogram(color="black",fill="lightsteelblue",binwidth=.005) + 
    theme_bw() + 
    xlab("Deviation from Mean\n Trade Dependency (Sender)") + 
    ylab("Density") + 
    theme(axis.title.x=element_text(size=14,vjust=-.2),
          axis.title.y=element_text(size=14),
          axis.text=element_text(size=12)
    )  
  ,
  ggplot(Mcowt, aes(x=Mcowt$targetecondepDev, y=..density..)) + 
    geom_histogram(color="black",fill="lightsteelblue",binwidth=.005) + 
    theme_bw() + 
    xlab("Deviation from Mean\n Trade Dependency (Target)") + 
    ylab("Density") + 
    theme(axis.title.x=element_text(size=14,vjust=-.2),
          axis.title.y=element_text(size=14),
          axis.text=element_text(size=12)
    )  , 
  
  
  ggplot(Mcowt, aes(x=Mcowt$targetcostDev, y=..density..)) + 
    geom_histogram(color="black",fill="lightsteelblue",binwidth=.005) + 
    theme_bw() + 
    xlab("Deviation from Mean\n Anticipated Target Costs") + 
    ylab("Density") + 
    theme(axis.title.x=element_text(size=14,vjust=-.2),
          axis.title.y=element_text(size=14),
          axis.text=element_text(size=12)
    )  ,
  
  ggplot(Mcowt, aes(x=Mcowt$lncapratDev, y=..density..)) + 
    geom_histogram(color="black",fill="lightsteelblue",binwidth=.1) + 
    theme_bw() + 
    xlab("Deviation from Mean\n Logged Capability Ratios") + 
    ylab("Density") + 
    theme(axis.title.x=element_text(size=14,vjust=-.2),
          axis.title.y=element_text(size=14),
          axis.text=element_text(size=12)
    )   ,
  
  
  ggplot(Mcowt, aes(x=Mcowt$allyDev, y=..density..)) + 
    geom_histogram(color="black",fill="lightsteelblue",binwidth=.005) + 
    theme_bw() + 
    xlab("Deviation from Mean Alliances") + 
    ylab("Density") + 
    theme(axis.title.x=element_text(size=14,vjust=-.2),
          axis.title.y=element_text(size=14),
          axis.text=element_text(size=12)
    )  , nrow=2
  
)
dev.off()

