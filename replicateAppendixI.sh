#!/bin/bash
Rscript appendixI1_estimation.R > appendixI.txt 2>&1 
python CMLE_sanctionsQuarters.py  > /dev/null 2>> appendixI.txt #suppress ipopt output
Rscript appendixI1_SE.R >> appendixI.txt 2>&1 
Rscript appendixI2_estimation.R >> appendixI.txt 2>&1 
Rscript appendixI3_5years.R >> appendixI.txt 2>&1
python CMLE_sanctions_5years.py  > /dev/null 2>> appendixI.txt #suppress ipopt output
Rscript appendixI3_5yearsSE.R >> appendixI.txt 2>&1  
Rscript appendixI3_1yearsT12.R >> appendixI.txt 2>&1
python CMLE_sanctions_1yearsT12.py  > /dev/null 2>> appendixI.txt #suppress ipopt output
Rscript appendixI3_1yearsT12_SE.R >> appendixI.txt 2>&1  
Rscript appendixI3_1yearsT1.R >> appendixI.txt 2>&1 
Rscript appendixI3_1yearsT1_SE.R >> appendixI.txt 2>&1 
