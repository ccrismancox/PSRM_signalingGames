#!/bin/bash
Rscript appendixF1.r > appendixF.txt 2>&1 #outputs figure
Rscript appendixF2.r >> appendixF.txt 2>&1 
Rscript appendixE_MEQ.r >> appendixE.txt 2>&1 
Rscript appendixE_Unique.r >> appendixE.txt 2>&1 

