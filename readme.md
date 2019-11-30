# Replication instructions for "Estimating Crisis Signaling Games in International Relations: Problems and Solutions"
# Casey Crisman-Cox and Michael Gibilisco
# October 15, 2019

## A note for replicators
Conducting constrained maximum likelihood estimation (CMLE)  requires specialized (open source) software that we run using a Ubuntu Linux operating system.
All our PL, NPL, and tML results are based solely on `R` and can be replicated without issue on any operating system.
Here, we provide instructions for reproducing our results using the Windows Subsystem for Linux (WSL) for Windows 10.
All the code should also work on an ordinary Ubuntu system.
We provide detailed setup instructions below.

	
## WSL Setup
All applications of the CMLE to either simulation or actual data was done with Ubuntu 18.04.1 (Bionic Beaver).
The automatic differentiation (AD) software is only tested for Ubuntu operating systems.
To allow for replication we provide these instructions for users using a Windows 10 computer with build 16215 or later (Fall Creators Update and later).
The installation steps are taken from this [link](https://docs.microsoft.com/en-us/windows/wsl/install-win10).


1. Enable WSL by opening the Windows PowerShell as an administrator and running the following command
```bash
Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux
```
2. Restart your computer 
3. Open the Windows store
4. Search for "Ubuntu"
5. Select "Ubuntu 18.04 LTS" and click "Get"
6. Once it is installed, open the Ubuntu app. This may take a few minutes to load the first time.
7. You will be prompted to choose a user name and password.  You will need to use this password throughout.
8. Once you are setup with the Ubuntu app, navigate to the replication folder from within the Ubuntu command line using the `cd` command.  For example, if the replication folder is saved in the Windows Downloads folder you could get there by using the command
```bash
cd /mnt/<c>/Users/<WINDOWS USER>/Downloads/<REPLICATION FOLDER>
```
Where `<c>` refers to your main Windows drive (almost always `c`), `<WINDOWS USER>` refers to your Windows username (**not** the Ubuntu user name you selected in step 7), and `<REPLICATION FOLDER>` is the folder containing all the replication files.

9. Run the file ``CMLE_setup.sh`` using the command
```bash
bash CMLE_setup.sh
```
This step may take up to a few hours depending on network speed and you may be prompted for your Ubuntu password (chosen in step 7), to select "yes",  or to press "Enter" at various points in the process. An Internet connection is required for this step.

10. Once this script has completed we are ready to reproduce all the results.  If Windows is configured for automatic updating, we recommend that you disconnect from the Internet during the Monte Carlo steps as they may take a few hours or days.


## Simulation results
To produce Figures 2-4 run bash script `replicateMonteCarlos.sh`.  Do this by opening the Ubuntu app, navigating to the Replication folder (as in step 8 above) and running the command:
```bash
bash replicateMonteCarlos.sh
```
This command produces four output files (`MonteCarloResults_MEQ.rdata`, `MonteCarloResults_Unique.rdata`, `CMLE_meq.rdata`, and `CMLE_unique.rdata`), a log file (`MonteCarloLog.txt`), and  three figures (`figure2.pdf`, `figure3.pdf`, and `figure4.pdf`).
These three PDF files correspond to Figures 2-4 in the manuscript, respectively.
Individual figures can be reproduced separately by running the commands in `replicateMonteCarlos.sh` one at a time.
Note that this may take several days or longer depending on the computer power available.


## Estimates, standard errors, and  comparative statics
To produce the values from Table 3 and the results in Figure 6 run the bash script `replicateSanctions.sh`
Do this by opening the Ubuntu app, navigating to the Replication folder (as in step 8 above) and running the command:
```bash
bash replicationSanctions.sh
```
This command produces three output files (`CMLE_estimation_output.rdata`, `estimation_output.Rdata`, and `SIGMA.rdata`), a log file (`estimationLog.txt`), and  one figure (`figure6.pdf`).
All the values reported in Table 3 are found in `estimationLog.txt`, and  Figure 6 is reproduced in `figure6.pdf`.
Individual aspects of the analysis can be reproduced separately by running the commands in `replicateSanctions.sh` one at a time.


## Appendices
To replicate the appendices run the commands:
```bash
bash replicateAppendixC.sh
bash replicateAppendixE.sh
bash replicateAppendixF.sh
bash replicateAppendixG.sh
bash replicateAppendixH.sh
bash replicateAppendixI.sh
```
These will produce Figures 8-24 and Tables 5-10. Table 5 is found in `appendixC.txt`, while Tables 6-10 are  in `appendixI.txt`. Additionally, various `rdata` files will be produced along the way.
