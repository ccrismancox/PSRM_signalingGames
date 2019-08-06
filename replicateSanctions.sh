#!/bin/bash
python CMLE_estimation.py > python_estimation.txt 2>&1
Rscript estimation.R > r_estimation.txt 2>&1
Rscript standardErrors.R > final_estimation.txt 2>&1
Rscript comparativeStatics.R > csOut.txt 2>&1 #outputs figure 
