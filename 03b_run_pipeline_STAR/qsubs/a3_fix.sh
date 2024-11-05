#!/bin/bash

#$ -V -cwd
#$ -l rl9=false
#$ -l h_rt=02:30:00 
#$ -l h_vmem=20G 
#$ -pe sharedmem 10 

. /etc/profile.d/modules.sh
 
module load anaconda/2024.02

conda activate seurat5

Rscript sub_lib3.R
