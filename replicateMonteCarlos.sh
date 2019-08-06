#!/bin/bash
Rscript eqGraph.r > /dev/null 2> /dev/null #outputs figure
Rscript MonteCarloMEQ.r > /dev/null 2> /dev/null
Rscript MonteCarloUnique.r > /dev/null 2> /dev/null
python CMLE_MonteCarlo_meq.py  > /dev/null 2> /dev/null
python CMLE_MonteCarlo_unique.py  > /dev/null 2> /dev/null
Rscript analyzeSimulationMEQ.r > /dev/null 2> /dev/null #outputs figure
Rscript analyzeSimulationUnique.r > /dev/null 2> /dev/null #outputs figure
