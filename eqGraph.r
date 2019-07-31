############################################
############################################
######## Estimating Signaling Games ########
########   Monte Carlo Experiment   ########
########    Graph of Eq. Corresp    ########
############################################
############################################
library(ggplot2)
library(pbivnorm)

source("signalingFunctions_main.r")


##### Multiple Setting #####
xold = -1
xnew = 2

Uold <- list(barWA = -1.9, barWB = -2.9 + 0.1*xold, bara = -1.2, VA = 1, VB = 1, SA = 0, CB = 0, sig=1)
Unew <- list(barWA = -1.9, barWB = -2.9 + 0.1*xnew, bara = -1.2, VA = 1, VB = 1, SA = 0, CB = 0, sig=1)

length.out <- 100
data <- signalEQC(Uold,Unew, length.out)
ggdata = data.frame(x = data[,1], eq = data[,2])
ggdata$select <- pmax(data[,1] %in% which(seq(xold, xnew, length.out=length.out) >=0 & seq(xold, xnew, length.out=length.out) < 1/3) & data[,2] < .39,
                      data[,1] %in% which(seq(xold, xnew, length.out=length.out) >=1/3 & seq(xold, xnew, length.out=length.out) < 2/3) & data[,2] > .39 & data[,2] < .71,
                      data[,1] %in% which(seq(xold, xnew, length.out=length.out) >=2/3 & seq(xold, xnew, length.out=length.out) <= 1) & data[,2] > .71)


ggdata.meq <- ggdata


##### Unique Setting #####

xold = -1
xnew = 2

Uold <- list(barWA = -1.8, barWB = -2.45 + 0.1*xold, bara = -1.2, VA = 1, VB = 1, SA = 0, CB = 0, sig=1)
Unew <- list(barWA = -1.8, barWB = -2.45 + 0.1*xnew, bara = -1.2, VA = 1, VB = 1, SA = 0, CB = 0, sig=1)

length.out <- 100
data <- signalEQC(Uold,Unew, length.out)
ggdata = data.frame(x = data[,1], eq = data[,2])
ggdata$select <- pmax(data[,1] %in% which(seq(xold, xnew, length.out=length.out) >=0 & seq(xold, xnew, length.out=length.out) <= 1))

ggdata$DGP <- "Unique"
ggdata.meq$DGP <- "Multiple"
fulldata <- rbind(ggdata, ggdata.meq)

pfull <- ggplot(fulldata, aes(x=x,y=eq)) +
  geom_point(aes(color=factor(select),shape=factor(select)),size=4) +
  scale_x_continuous(breaks=seq(0, length.out,length.out=10), 
                     labels = round(seq(-1,2,length.out=10),digits=2),
                     limits = c(14,86)) +
  theme_bw(30) + 
  theme(legend.title=element_blank()) +
  scale_colour_manual(values=c("navyblue","orangered1"),labels=c("Not selected", "Selected")) + 
  xlab(expression(paste("Regressor ", italic(x[d])))) + 
  ylab(expression(paste("Equilibria, ", italic(p[dR])))) + 
  scale_shape_discrete(labels=c("Not selected", "Selected")) +
  facet_wrap(~DGP)+
  guides(color = guide_legend(keywidth=2, override.aes = list(size=5)))+
  theme(legend.position="bottom",
        legend.text = element_text(size = 25),
        panel.spacing = unit(2.5, "lines"))
ggsave(filename="figure2.pdf", plot=pfull, width=15, height=6)

warning("End of file. Press enter if the system hangs here.")
