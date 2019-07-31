#!/bin/bash
python CMLE_estimation.py > python_estimation.txt
Rscript estimation.R > r_estimation.txt
Rscript standardErrors.R > final_estimation.txt
Rscript comparativeStatics.R #outputs figure
