#!/bin/bash

#$ -V -cwd
#$ -l h_rt=08:00:00
#$ -l h_vmem=40G
#$ -pe sharedmem 1


module load roslin/python/3.8.10

SCRIPTPATH="/exports/eddie/scratch/pdewari/concat/fastq_sep_groups_new2024.py" 

FQ_DIR="/exports/eddie/scratch/pdewari/concat/"

python $SCRIPTPATH \
--chemistry v2 \
--fq1 ${FQ_DIR}A8_EKDL240002480-1A_223M7CLT4_1.fq.gz \
--fq2 ${FQ_DIR}A8_EKDL240002480-1A_223M7CLT4_2.fq.gz \
--opath ${FQ_DIR}split-fq \
--group shrimp A1-C8 \
--group oyster C9-D12
