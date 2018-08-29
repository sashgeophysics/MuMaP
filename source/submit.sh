#!/bin/bash --login
#
#PBS -l select=1
#PBS -o mumap.out
#PBS -e mumap.error

#PBS -l walltime=0:10:00
#PBS -m be
#PBS -M saswata.hier-majumder@rhul.ac.uk
#PBS -A n03-shm
cd /work/n03/n03/saswata/MuMAP/source

aprun -n 24 -d 1 agius_hi.exe
