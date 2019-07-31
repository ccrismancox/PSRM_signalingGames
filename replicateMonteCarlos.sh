#!/bin/bash
Rscript eqGraph.r #outputs figure
Rscript MonteCarloMEQ.r > r_montecarlo.txt
Rscript MonteCarloUnique.r >> r_montecarlo.txt
python CMLE_MonteCarlo_meq.py > python_montecarlo.txt
python CMLE_MonteCarlo_unique.py >> python_montecarlo.txt
Rscript analyzeSimulationMEQ.r  #outputs figure
Rscript analyzeSimulationUnique.r #outputs figure
