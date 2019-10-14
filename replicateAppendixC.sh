#!/bin/bash
Rscript appendixC1.r > appendixC.txt 2>&1 #outputs figures
Rscript appendixC2.r >> appendixC.txt 2>&1 #outputs figures
Rscript appendixC3_simulation.r >> appendixC.txt 2>&1
python MonteCarloCMLE_unstable.py > /dev/null 2>> appendixC.txt #suppress ipopt output
Rscript appendixC3.R >> appendixC.txt 2>&1 #outputs figures
Rscript appendixC4_simulation.R >> appendixC.txt 2>&1
python CMLE_sanctionsMC.py > /dev/null 2>> appendixC.txt #suppress ipopt output
Rscript appendixC4.R >> appendixC.txt 2>&1 #outputs table


