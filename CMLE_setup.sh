#!/bin/bash

## DESCRIPTION:
#### In this file we setup all the software
#### needed to replicate our results
#### Note that in setting up Windows 10's Windows Subsystem for Linux (WSL) 
#### you will need to create a user name and password.  That password will 
#### be required several times. For your ease and to ensure full replication,
#### please use the WSL with Ubuntu 18.04.1 (Bionic Beaver). 
#### For with replicating our results on a regular Ubuntu system, contact that authors.
#### Other operating systems are not supported.
#### Detailed instructions on setting up WSL are found in the readme file.


## INSTRUCTIONS: 
#### After setting up WSL. Open the Ubuntu app for Windows 10
#### and navigate to the replication folder.
#### For example, if the replication folder is in your Windows
#### Downloads folder you should use the command
#### cd /mnt/c/Users/USERNAME/Downloads
#### where USERNAME is your Windows user folder name.
#### Once in the replication folder, run this file using the command
#### bash CMLE_setup.sh

## AUTHORS: Casey Crisman-Cox and Michael Gibilisco
##### Created for use Windows 10's Windows Subsystem for Linux (WSL) 
##### running with Ubuntu 18.04.1.
##### Not guaranteed for any other system or version.

## Basic setup and tools
## This list includes 
#### TeX (required for pyadolc and others)
#### Python with packages
#### R with packages
#### basic compilers to build C, C++, and Fotran
#### Git and version control software
#### other dependencies

REPDIR=`pwd`
cd ~
HOMEDIR=`pwd`
sudo apt update
sudo apt -y upgrade


sudo apt -y install build-essential
sudo apt -y install python-dev  git 
sudo apt -y install texlive-full
sudo apt -y install gfortran automake shtool libtool
sudo apt -y install python-matplotlib python-scipy python-pandas python-sympy python-nose spyder
sudo apt -y install subversion swig
sudo apt -y install openmpi-bin openmpi-doc libopenmpi-dev
sudo apt -y  install r-base-dev
mkdir -p ~/R/x86_64-pc-linux-gnu-library/3.4

sudo apt -y install libboost-all-dev
sudo apt -y install libboost-python-dev
sudo apt -y install libboost-system-dev

# Additional python packages
sudo pip install numpy --upgrade
sudo pip install rpy2==2.8.6
sudo pip install mpi4py

# ADOLC, Colpack, and IPOPT, use code stored on my bitbucket to help here
git clone https://ccrismancox@bitbucket.org/ccrismancox/pyopterf_windows.git pyopterf
cd pyopterf

bash setup.sh

libdir=${HOMEDIR}/pyopterf/Ipopt-3.12.3/lib
sudo echo $libdir$'\r' | sudo tee -a  /etc/ld.so.conf
sudo ldconfig
sudo apt -y install libgfortran3 #where does this go?


#Additional R work
sudo apt -y install  libcurl4-gnutls-dev libxml2-dev libssl-dev
cd ~
mkdir -p Documents
Rscript "${REPDIR}/install_Rpackages_noadmin.r"
echo "Ubuntu Setup Complete"







