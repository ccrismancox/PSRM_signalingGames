############################################
############################################
######## Estimating Signaling Games ########
########       Install Packages     ########
############################################
############################################
#####R TOOLS 35 REQUIRED#####




install.packages("devtools")
library(devtools)





#####Assuming R tools is in the usual place,  you may need to adjust if it isn't####
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) 


##### Needed for the main packages #####
install_version("mvtnorm","1.0-8",repos="https://cloud.r-project.org/")
install_version("rngtools", "1.3.1",repos="https://cloud.r-project.org/")


##### main packages ####
install.packages("data.table")
install.packages("doParallel")
install.packages("doRNG")
install.packages("e1071")
install.packages("Formula")
install.packages("foreach")
install.packages("foreign")
install.packages("ggplot2")
install.packages("Matrix")
install.packages("matrixStats")
install.packages("maxLik")
install.packages("mc2d")
install.packages("numDeriv")
install.packages("pbivnorm")
install.packages("randomForest")
install.packages("rootSolve")
install.packages("scales")
