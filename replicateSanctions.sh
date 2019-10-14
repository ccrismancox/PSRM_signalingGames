#!/bin/bash
Rscript estimation.R > estimationLog.txt 2>&1
python CMLE_estimation.py >> estimationLog.txt 2>&1
Rscript standardErrors.R >> estimationLog.txt 2>&1
Rscript comparativeStatics.R >> estimationLog.txt 2>&1 #outputs figure 
