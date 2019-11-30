############################################
############################################
######## Estimating Signaling Games ########
########       Install Packages     ########
############################################
############################################
#####R TOOLS 35 REQUIRED#####




install.packages("devtools", repos = "https://cloud.r-project.org/")
library(devtools)





#####Assuming R tools is in the usual place,  you may need to adjust if it isn't####
Sys.setenv(PATH = paste("C:/Rtools/bin", Sys.getenv("PATH"), sep=";")) 


##### Needed for the main packages #####
install_version("mvtnorm","1.0-8",repos="https://cloud.r-project.org/")
install_version("rngtools", "1.3.1",repos="https://cloud.r-project.org/")


##### main packages ####
install.packages("data.table", repos = "https://cloud.r-project.org/")
install.packages("doParallel", repos = "https://cloud.r-project.org/")
install.packages("doRNG", repos = "https://cloud.r-project.org/")
install.packages("e1071", repos = "https://cloud.r-project.org/")
install.packages("Formula", repos = "https://cloud.r-project.org/")
install.packages("foreach", repos = "https://cloud.r-project.org/")
install.packages("foreign", repos = "https://cloud.r-project.org/")
install.packages("ggplot2", repos = "https://cloud.r-project.org/")
install.packages("gridExtra", repos = "https://cloud.r-project.org/")
install.packages("Matrix", repos = "https://cloud.r-project.org/")
install.packages("matrixStats", repos = "https://cloud.r-project.org/")
install.packages("maxLik", repos = "https://cloud.r-project.org/")
install.packages("mc2d", repos = "https://cloud.r-project.org/")
install.packages("numDeriv", repos = "https://cloud.r-project.org/")
install.packages("pbivnorm", repos = "https://cloud.r-project.org/")
install.packages("randomForest", repos = "https://cloud.r-project.org/")
install.packages("rootSolve", repos = "https://cloud.r-project.org/")
install.packages("scales", repos = "https://cloud.r-project.org/")
install.packages("./sigInt_0.0.0.9000.tar.gz", repos=NULL, type="source")