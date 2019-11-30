#!/bin/bash
Rscript appendixE_simulation_MEQ.r > appendixE.txt 2>&1 
Rscript appendixE_simulation_Unique.r >> appendixE.txt 2>&1 
Rscript appendixE_MEQ.r >> appendixE.txt 2>&1 #outputs figure
Rscript appendixE_Unique.r >> appendixE.txt 2>&1 #outputs figure

