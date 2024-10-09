#!/bin/bash -l
#$ -S /bin/bash
# Set wallclock time (format hours:minutes:seconds)
#$ -l h_rt=02:00:00
#$ -l mem=1G
#$ -N UFF_FQ_2
#$ -cwd 

source ~/.bashrc
conda activate mdsimul
python CONNECT2GRO.py 
