#!/bin/bash

#$ -V -cwd
#$ -l rl9=false
#$ -l h_rt=20:00:00
#$ -l h_vmem=60G
#$ -pe sharedmem 6

. /etc/profile.d/modules.sh

module load anaconda/2024.02

conda activate seurat5

Rscript Merging_Barcodes_for_STARSOLO_analysis.R
