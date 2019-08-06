#!/bin/bash
Rscript eqGraph.r > MonteCarloLog.txt 2>&1   #outputs figure
Rscript MonteCarloMEQ.r >> MonteCarloLog.txt 2>&1 
Rscript MonteCarloUnique.r >> MonteCarloLog.txt 2>&1 
python CMLE_MonteCarlo_meq.py  > /dev/null 2>> MonteCarloLog.txt
python CMLE_MonteCarlo_unique.py  > /dev/null 2>> MonteCarloLog.txt
Rscript analyzeSimulationMEQ.r >> MonteCarloLog.txt 2>&1 #outputs figure
Rscript analyzeSimulationUnique.r >> MonteCarloLog.txt 2>&1 #outputs figure
